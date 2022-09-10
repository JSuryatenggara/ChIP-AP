#!/usr/bin/env python3
#pyright: reportUnboundVariable=false


script_version = '2.0'


import tkinter as tk
import tkinter.font as tkFont
from tkinter import filedialog
import os
import multiprocessing
import subprocess
import pandas as pd
from PIL import Image, ImageTk
import sys


chipap_program_name = 'chipap.py'

chipap_icon_full_path = os.path.expanduser('{}/ChIP-AP_icon_GUI.png'.format(sys.path[0])) # Path to ChIP-AP mini icon (png)
chipap_logo_full_path = os.path.expanduser('{}/ChIP-AP_logo_GUI.jpg'.format(sys.path[0])) # Path to ChIP-AP full logo (jpeg)
default_current_dir = os.path.expanduser('~') # Default starting directory for when user choose to browse
genome_folder_full_path = os.path.expanduser('~/genomes') # Path to the default genome folder
default_setting_table_file_full_path = '{}/default_settings_table.tsv'.format(genome_folder_full_path) # Path to the default setting table)

# Create a non-manually-resizable GUI base frame (root)
root = tk.Tk()
root.title(chipap_program_name)
root.resizable(width = False, height = False)


# Set the default fonts to be used in the GUI
default_font = tkFont.nametofont("TkDefaultFont")
default_font.configure(family = 'fixed', size = 20)

text_font = tkFont.nametofont("TkTextFont")
text_font.configure(family = 'fixed', size = 20)

fixed_font = tkFont.nametofont("TkFixedFont")
fixed_font.configure(family = 'fixed', size = 20)

# Set the provided ChIP-AP icon (in the same folder as this script) as the window and taskbar icon for the GUI
root.tk.call('wm', 'iconphoto', root._w, tk.PhotoImage(file = chipap_icon_full_path))


argument_dict = {} # Empty dictionary for programs setting values

# List of accepted extensions of ChIP and control sample files
valid_extension_list = [".fastq",
                        ".fq",
                        ".fastq.gz",
                        ".fq.gz",
                        ".bam"]

# List of program names which setting values are customizable (keys to argument_dict)                        
suite_program_list = ["fastqc1",
                        "clumpify",
                        "bbduk",
                        "trimmomatic",
                        "fastqc2",
                        "bwa_mem",
                        "samtools_view",
                        "plotfingerprint",
                        "fastqc3",
                        "macs2_callpeak",
                        "gem",
                        "sicer2",
                        "homer_findPeaks",
                        "genrich",
                        "homer_mergePeaks",
                        "homer_annotatePeaks",
                        "fold_change_calculator",
                        "homer_findMotifsGenome",
                        "meme_chip"]

# List of supported reference genomes for aligning by bwa mem, etc (options for genome_ref_drop_down menu)
genome_ref_options = ["hg38 (Homo sapiens build 38)", 
                        "hg19 (Homo sapiens build 19)", 
                        "mm9 (Mus musculus build 9)", 
                        "mm10 (Mus musculus build 10)", 
                        "dm6 (Drosophila melanogaster build 6)", 
                        "sacCer3 (Saccharomyces cerevisiae build 3)",
                        "Other [!!!under construction!!!]"]

# List of supported choices of peak list to analyze for motif enrichment by HOMER(options for homer_motif_drop_down menu)
homer_motif_options = ["<Please select peak set(s)>",
                        "Consensus peak set", 
                        "Union peak set", 
                        "Both peak sets"]

# List of supported choices of peak list to analyze for motif enrichment by MEME (options for meme_motif_drop_down menu)
meme_motif_options = ["<Please select peak set(s)>",
                        "Consensus peak set", 
                        "Union peak set", 
                        "Both peak sets"]


########################################################################################################################


# Define tooltip object properties for the GUI
# Copied this python class from some forum then modified it, thus the lack of comments :)
class ToolTip(object):

    def __init__(self, widget):
        self.widget = widget
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0

    def showtip(self, text):
        "Display text in tooltip window"
        self.text = text
        if self.tipwindow or not self.text:
            return
        x, y, cx, cy = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 150
        y = y + cy + self.widget.winfo_rooty() + 45
        self.tipwindow = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        tw.wm_geometry("+%d+%d" % (x, y))
        label = tk.Label(tw, text = self.text, justify = tk.LEFT,
                      background = "#ffffe0", relief = tk.SOLID, borderwidth = 1,
                      font = ("tahoma", "16", "normal"))
        label.pack(ipadx = 1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()

# Tooltip objects are designed to appear at mouse hover or widget focusing with tab key
def CreateToolTip(widget, text):
    toolTip = ToolTip(widget)
    def enter(event):
        toolTip.showtip(text)
    def leave(event):
        toolTip.hidetip()
    widget.bind('<Enter>', enter)
    widget.bind('<Leave>', leave)
    widget.bind('<FocusIn>', enter)
    widget.bind('<FocusOut>', leave)


# Define the function to switch between frames
def change_frame_function(from_frame, to_frame):

    # If sample table loading fails, function is exited prematurely and wizard does not proceed to the ChIP samples frame
    if from_frame == sample_table_frame and to_frame == chip_list_frame and sample_table_question_string_var.get() == 'yes':
        if sample_table_loading_test_function() == False:
            return
        elif sample_table_loading_test_function() == True:
            pass

    # If setting table loading fails (unless setting values modification present),
    #   function is exited prematurely and GUI does not proceed to the customize settings frame
    if from_frame == setting_table_frame and (to_frame == setting_value_frame or to_frame == genome_ref_frame):
        if '[MODIFIED]' in setting_table_string_var.get():
            pass
        else:
            if setting_table_loading_test_function() == False:
                return
            elif setting_table_loading_test_function() == True:
                pass
    
    # Detect setting values modification by comparing values in the entry fields to the values in the loaded settings table
    if from_frame == setting_value_frame and to_frame == setting_table_frame:
        current_setting_table_full_path = setting_table_string_var.get()

        if fastqc1_arg.get() != ' '.join(argument_dict['fastqc1']):
            setting_changed = 1
        elif clumpify_arg.get() != ' '.join(argument_dict['clumpify']):
            setting_changed = 1
        elif bbduk_arg.get() != ' '.join(argument_dict['bbduk']):
            setting_changed = 1
        elif trimmomatic_arg.get() != ' '.join(argument_dict['trimmomatic']):
            setting_changed = 1
        elif fastqc2_arg.get() != ' '.join(argument_dict['fastqc2']):
            setting_changed = 1
        elif bwa_mem_arg.get() != ' '.join(argument_dict['bwa_mem']):
            setting_changed = 1
        elif samtools_view_arg.get() != ' '.join(argument_dict['samtools_view']):
            setting_changed = 1
        elif plotfingerprint_arg.get() != ' '.join(argument_dict['plotfingerprint']):
            setting_changed = 1
        elif fastqc3_arg.get() != ' '.join(argument_dict['fastqc3']):
            setting_changed = 1
        elif macs2_callpeak_arg.get() != ' '.join(argument_dict['macs2_callpeak']):
            setting_changed = 1
        elif gem_arg.get() != ' '.join(argument_dict['gem']):
            setting_changed = 1
        elif sicer2_arg.get() != ' '.join(argument_dict['sicer2']):
            setting_changed = 1
        elif homer_findPeaks_arg.get() != ' '.join(argument_dict['homer_findPeaks']):
            setting_changed = 1
        elif genrich_arg.get() != ' '.join(argument_dict['genrich']):
            setting_changed = 1
        elif homer_mergePeaks_arg.get() != ' '.join(argument_dict['homer_mergePeaks']):
            setting_changed = 1
        elif homer_annotatePeaks_arg.get() != ' '.join(argument_dict['homer_annotatePeaks']):
            setting_changed = 1
        elif fold_change_calculator_arg.get() != ' '.join(argument_dict['fold_change_calculator']):
            setting_changed = 1
        elif homer_findMotifsGenome_arg.get() != ' '.join(argument_dict['homer_findMotifsGenome']):
            setting_changed = 1
        elif meme_chip_arg.get() != ' '.join(argument_dict['meme_chip']):
            setting_changed = 1
        else:
            setting_changed = 0

        # If any of setting values are modified, add a modification marker to the setting table file name and display notification to user
        if setting_changed == 1:
            setting_table_notification_string_var.set('Settings were manually modified by user')
            setting_table_notification_label.config(fg = 'orange2')
            if '[MODIFIED]' not in setting_table_string_var.get(): 
                setting_table_string_var.set('{}{}'.format(current_setting_table_full_path, '[MODIFIED]'))
        elif setting_changed == 0:
            pass

    # Set to skip the HOMER and MEME motif enrichment analysis dialogue frames if the dataset peak type is broad
    if peak_type_string_var.get() == 'broad':
        if from_frame == output_folder_frame and to_frame == homer_motif_frame:
            to_frame = checkbox_frame
        if from_frame == checkbox_frame and to_frame == meme_motif_frame:
            to_frame = output_folder_frame

    # Finally, if the function managed to proceed up to this point, the current GUI frame (from_frame) is going to be colapsed
    #   and replaced with the destination frame (to_frame)
    from_frame.grid_remove()
    to_frame.grid(row = 0, column = 0, sticky = "news")
    root.eval('tk::PlaceWindow . center')
    root.grab_set()
    root.focus_force()

    if to_frame == chip_list_frame: # Trigger to display the right number of ChIP sample entry fields upon entering chip_list_frame
        display_chip_widget_function()

    if to_frame == ctrl_list_frame: # Trigger to display the right number of control sample entry fields upon entering ctrl_list_frame
        display_ctrl_widget_function()

    # Functions to keep the frame window (and its contents) at the center of screen
    root.eval('tk::PlaceWindow . center')
    root.grab_set()
    root.focus_force()
    
    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path


# Define a function to rearrange the entry fields such that they always display the last n characters of entered path,
#   makes it easier for user to check and confirm for the correct files/directories
def rearrange_entry_field_function():
    chip_rep1_r1_entry.xview_moveto(1)
    chip_rep1_r2_entry.xview_moveto(1)
    chip_rep2_r1_entry.xview_moveto(1)
    chip_rep2_r2_entry.xview_moveto(1)
    chip_rep3_r1_entry.xview_moveto(1)
    chip_rep3_r2_entry.xview_moveto(1)
    chip_rep4_r1_entry.xview_moveto(1)
    chip_rep4_r2_entry.xview_moveto(1)
    chip_rep5_r1_entry.xview_moveto(1)
    chip_rep5_r2_entry.xview_moveto(1)

    ctrl_rep1_r1_entry.xview_moveto(1)
    ctrl_rep1_r2_entry.xview_moveto(1)
    ctrl_rep2_r1_entry.xview_moveto(1)
    ctrl_rep2_r2_entry.xview_moveto(1)
    ctrl_rep3_r1_entry.xview_moveto(1)
    ctrl_rep3_r2_entry.xview_moveto(1)
    ctrl_rep4_r1_entry.xview_moveto(1)
    ctrl_rep4_r2_entry.xview_moveto(1)
    ctrl_rep5_r1_entry.xview_moveto(1)
    ctrl_rep5_r2_entry.xview_moveto(1)

    sample_table_entry.xview_moveto(1)
    setting_table_entry.xview_moveto(1)
    genome_folder_entry.xview_moveto(1)
    known_motif_entry.xview_moveto(1)
    output_folder_entry.xview_moveto(1)


########################################################################################################################


# Function to register all the samples stored GUI variables to their respective ChIP or control sample list
#   Subsequently create a python dictionary out of these lists for generation of the new sample table
def register_sample_function(*args):
    
    # Set white as the default background color for the entry field
    chip_rep1_r1_entry.config(bg = 'white')
    chip_rep1_r2_entry.config(bg = 'white')
    chip_rep2_r1_entry.config(bg = 'white')
    chip_rep2_r2_entry.config(bg = 'white')
    chip_rep3_r1_entry.config(bg = 'white')
    chip_rep3_r2_entry.config(bg = 'white')
    chip_rep4_r1_entry.config(bg = 'white')
    chip_rep4_r2_entry.config(bg = 'white')
    chip_rep5_r1_entry.config(bg = 'white')
    chip_rep5_r2_entry.config(bg = 'white')

    ctrl_rep1_r1_entry.config(bg = 'white')
    ctrl_rep1_r2_entry.config(bg = 'white')
    ctrl_rep2_r1_entry.config(bg = 'white')
    ctrl_rep2_r2_entry.config(bg = 'white')
    ctrl_rep3_r1_entry.config(bg = 'white')
    ctrl_rep3_r2_entry.config(bg = 'white')
    ctrl_rep4_r1_entry.config(bg = 'white')
    ctrl_rep4_r2_entry.config(bg = 'white')
    ctrl_rep5_r1_entry.config(bg = 'white')
    ctrl_rep5_r2_entry.config(bg = 'white')

    # Set black as the default font color for the entry field
    chip_rep1_r1_entry.config(fg = 'black')
    chip_rep1_r2_entry.config(fg = 'black')
    chip_rep2_r1_entry.config(fg = 'black')
    chip_rep2_r2_entry.config(fg = 'black')
    chip_rep3_r1_entry.config(fg = 'black')
    chip_rep3_r2_entry.config(fg = 'black')
    chip_rep4_r1_entry.config(fg = 'black')
    chip_rep4_r2_entry.config(fg = 'black')
    chip_rep5_r1_entry.config(fg = 'black')
    chip_rep5_r2_entry.config(fg = 'black')

    ctrl_rep1_r1_entry.config(fg = 'black')
    ctrl_rep1_r2_entry.config(fg = 'black')
    ctrl_rep2_r1_entry.config(fg = 'black')
    ctrl_rep2_r2_entry.config(fg = 'black')
    ctrl_rep3_r1_entry.config(fg = 'black')
    ctrl_rep3_r2_entry.config(fg = 'black')
    ctrl_rep4_r1_entry.config(fg = 'black')
    ctrl_rep4_r2_entry.config(fg = 'black')
    ctrl_rep5_r1_entry.config(fg = 'black')
    ctrl_rep5_r2_entry.config(fg = 'black')

    # Prepare global lists for all ChIP and control samples
    global chip_list_r1
    global chip_list_r2
    global ctrl_list_r1
    global ctrl_list_r2
    
    # Intial reset of all ChIP and control samples lists
    chip_list_r1 = []
    chip_list_r2 = []
    ctrl_list_r1 = []
    ctrl_list_r2 = []

    # Get the ChIP and control samples registered input and try to append it to the correct lists
    # Sequentially, starting from the first replicate to the fifth replicate
    # If a process fails for any reason, proceed and try the subsequent process without returning an error
    try:
        chip_list_r1.append(chip_rep1_r1_string_var.get())
    except:
        pass
    try:
        chip_list_r2.append(chip_rep1_r2_string_var.get())
    except:
        pass
    try:
        chip_list_r1.append(chip_rep2_r1_string_var.get())
    except:
        pass
    try:
        chip_list_r2.append(chip_rep2_r2_string_var.get())
    except:
        pass
    try:
        chip_list_r1.append(chip_rep3_r1_string_var.get())
    except:
        pass
    try:
        chip_list_r2.append(chip_rep3_r2_string_var.get())
    except:
        pass
    try:
        chip_list_r1.append(chip_rep4_r1_string_var.get())
    except:
        pass
    try:
        chip_list_r2.append(chip_rep4_r2_string_var.get())
    except:
        pass
    try:
        chip_list_r1.append(chip_rep5_r1_string_var.get())
    except:
        pass
    try:
        chip_list_r2.append(chip_rep5_r2_string_var.get())
    except:
        pass
    try:
        ctrl_list_r1.append(ctrl_rep1_r1_string_var.get())
    except:
        pass
    try:
        ctrl_list_r2.append(ctrl_rep1_r2_string_var.get())
    except:
        pass
    try:
        ctrl_list_r1.append(ctrl_rep2_r1_string_var.get())
    except:
        pass
    try:
        ctrl_list_r2.append(ctrl_rep2_r2_string_var.get())
    except:
        pass
    try:
        ctrl_list_r1.append(ctrl_rep3_r1_string_var.get())
    except:
        pass
    try:
        ctrl_list_r2.append(ctrl_rep3_r2_string_var.get())
    except:
        pass
    try:
        ctrl_list_r1.append(ctrl_rep4_r1_string_var.get())
    except:
        pass
    try:
        ctrl_list_r2.append(ctrl_rep4_r2_string_var.get())
    except:
        pass
    try:
        ctrl_list_r1.append(ctrl_rep5_r1_string_var.get())
    except:
        pass
    try:
        ctrl_list_r2.append(ctrl_rep5_r2_string_var.get())
    except:
        pass

    # Prepare a global dictionary variable for the sample table to be generated later
    global sample_table_output_dict
    
    # Assign the appended ChIP and control lists above as the values of the dictionary
    sample_table_output_dict = {'chip_read_1' : chip_list_r1,
                                'chip_read_2' : chip_list_r2,
                                'ctrl_read_1' : ctrl_list_r1,
                                'ctrl_read_2' : ctrl_list_r2}

    check_chip_sample_assigned_function() # First step of ChIP samples error checking: Check to see if they are at least assigned
    check_ctrl_sample_assigned_function() # First step of control samples error checking: Check to see if they are at least assigned
    check_required_input_function() # Check if the required inputs other than sample errors are satisifed or contain any errors
    update_command_line_function() # Update the generated ChIP-AP command line based on the latest registered variable values
    

# First step of ChIP samples error checking: Check to see if they are at least assigned
def check_chip_sample_assigned_function():

    if read_mode_string_var.get() == 'single' or read_mode_string_var.get() == 'paired': # Perform checks if dataset sequencing mode is single-end or paired-end
        
        if read_mode_string_var.get() == 'single': # In case of single-end dataset
            if all(chip_r1 == '' for chip_r1 in chip_list_r1): # If r1 list is all empty strings
                chip_sample_notification_string_var.set('Please assign ChIP sample') # Notify the user
                chip_sample_notification_label.config(fg = 'blue')
                chip_list_frame_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
                return # Exit the function

        if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
            if all(chip_r1 == '' for chip_r1 in chip_list_r1): # If r1 list is all empty strings
                chip_sample_notification_string_var.set('Please assign ChIP sample (read 1)') # Notify the user
                chip_sample_notification_label.config(fg = 'blue')
                chip_list_frame_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
                return # Exit the function

            # Check for r2 only if the input files are not aligned
            if all(chip_r2 == '' for chip_r2 in chip_list_r2) and not all('.bam' in sample for sample in (chip_list_r1 + chip_list_r1) if sample != ''):  # If r2 list is all empty strings
                chip_sample_notification_string_var.set('Please assign ChIP sample (read 2)') # Notify the user
                chip_sample_notification_label.config(fg = 'blue')
                chip_list_frame_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
                return # Exit the function

        # If the function reaches this point, it means there are at least one sample assigned to each of the lists necessary

        if read_mode_string_var.get() == 'single': # In case of single-end dataset
            if not bool(chip_rep1_r1_string_var.get()): # If the assigned sample(s) is not the first replicate
                chip_sample_notification_string_var.set('ChIP samples are not assigned from the top (replicate 1)') # Notify the user
                chip_sample_notification_label.config(fg = 'red')
                chip_rep1_r1_entry.config(bg = 'IndianRed1') # Highlight where the problem is
                chip_list_frame_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
                return # Exit the function
        
        if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
            if not bool(chip_rep1_r1_string_var.get()) and not bool(chip_rep1_r2_string_var.get()): # If the assigned sample(s) is not the first replicate of r1 nor r2
                chip_sample_notification_string_var.set('ChIP samples are not assigned from the top (replicate 1)') # Notify the user
                chip_sample_notification_label.config(fg = 'red')
                chip_rep1_r1_entry.config(bg = 'IndianRed1') # Highlight where the problem is
                chip_rep1_r2_entry.config(bg = 'IndianRed1') # Highlight where the problem is
                chip_list_frame_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
                return # Exit the function

            if not bool(chip_rep1_r1_string_var.get()): # If the assigned sample(s) is not the first replicate of r1
                chip_sample_notification_string_var.set('ChIP samples are not assigned from the top (replicate 1)') # Notify the user
                chip_sample_notification_label.config(fg = 'red')
                chip_rep1_r1_entry.config(bg = 'IndianRed1') # Highlight where the problem is
                chip_list_frame_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
                return # Exit the function

            if not bool(chip_rep1_r2_string_var.get()): # If the assigned sample(s) is not the first replicate of r2
                chip_sample_notification_string_var.set('ChIP samples are not assigned from the top (replicate 1)') # Notify the user
                chip_sample_notification_label.config(fg = 'red')
                chip_rep1_r2_entry.config(bg = 'IndianRed1') # Highlight where the problem is
                chip_list_frame_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
                return # Exit the function

        # If no errors found and the function did not exit prematurely because of all the tests above:
        chip_sample_notification_string_var.set('ChIP samples have been assigned') # Notify the user
        chip_sample_notification_label.config(fg = 'green')
        chip_list_frame_continue_button.config(state = tk.NORMAL) # Enable the continue to next frame button
    
        check_chip_sample_file_pair_function() # Second step of ChIP samples error checking: Check to see if paired-end dataset samples are properly paired

    else: # If dataset sequencing mode is not single-end nor paired-end
        return # Exit the function


# First step of control samples error checking: Check to see if they are at least assigned
def check_ctrl_sample_assigned_function():

    if read_mode_string_var.get() == 'single' or read_mode_string_var.get() == 'paired': # Perform checks if dataset sequencing mode is single-end or paired-end
        
        if read_mode_string_var.get() == 'single': # In case of single-end dataset
            if all(ctrl_r1 == '' for ctrl_r1 in ctrl_list_r1): # If r1 list is all empty strings
                ctrl_sample_notification_string_var.set('Please assign control sample') # Notify the user
                ctrl_sample_notification_label.config(fg = 'blue')
                ctrl_list_frame_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
                return

        if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
            if all(ctrl_r1 == '' for ctrl_r1 in ctrl_list_r1): # If r1 list is all empty strings
                ctrl_sample_notification_string_var.set('Please assign control sample (read 1)') # Notify the user
                ctrl_sample_notification_label.config(fg = 'blue')
                ctrl_list_frame_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
                return

            # Check for r2 only if the input files are not aligned
            if all(ctrl_r2 == '' for ctrl_r2 in ctrl_list_r2) and not all('.bam' in sample for sample in (ctrl_list_r1 + ctrl_list_r1) if sample != ''): # If r2 list is all empty strings
                ctrl_sample_notification_string_var.set('Please assign control sample (read 2)') # Notify the user
                ctrl_sample_notification_label.config(fg = 'blue')
                ctrl_list_frame_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
                return

        # If the function reaches this point, it means there are at least one sample assigned to each of the lists necessary

        if read_mode_string_var.get() == 'single': # In case of single-end dataset
            if not bool(ctrl_rep1_r1_string_var.get()): # If the assigned sample(s) is not the first replicate
                ctrl_sample_notification_string_var.set('Control samples are not assigned from the top (replicate 1)') # Notify the user
                ctrl_sample_notification_label.config(fg = 'red')
                ctrl_rep1_r1_entry.config(bg = 'IndianRed1') # Highlight where the problem is
                ctrl_list_frame_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
                return
        
        if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
            if not bool(ctrl_rep1_r1_string_var.get()) and not bool(ctrl_rep1_r2_string_var.get()): # If the assigned sample(s) is not the first replicate of r1 nor r2
                ctrl_sample_notification_string_var.set('Control samples are not assigned from the top (replicate 1)') # Notify the user
                ctrl_sample_notification_label.config(fg = 'red')
                ctrl_rep1_r1_entry.config(bg = 'IndianRed1') # Highlight where the problem is
                ctrl_rep1_r2_entry.config(bg = 'IndianRed1') # Highlight where the problem is
                ctrl_list_frame_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
                return

            if not bool(ctrl_rep1_r1_string_var.get()): # If the assigned sample(s) is not the first replicate of r1
                ctrl_sample_notification_string_var.set('Control samples are not assigned from the top (replicate 1)') # Notify the user
                ctrl_sample_notification_label.config(fg = 'red')
                ctrl_rep1_r1_entry.config(bg = 'IndianRed1') # Highlight where the problem is
                ctrl_list_frame_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
                return

            if not bool(ctrl_rep1_r2_string_var.get()): # If the assigned sample(s) is not the first replicate of r2
                ctrl_sample_notification_string_var.set('Control samples are not assigned from the top (replicate 1)') # Notify the user
                ctrl_sample_notification_label.config(fg = 'red')
                ctrl_rep1_r2_entry.config(bg = 'IndianRed1') # Highlight where the problem is
                ctrl_list_frame_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
                return

        # If no errors found and the function did not exit prematurely because of all the tests above:
        ctrl_sample_notification_string_var.set('Control samples have been assigned') # Notify the user
        ctrl_sample_notification_label.config(fg = 'green')
        ctrl_list_frame_continue_button.config(state = tk.NORMAL) # Enable the continue to next frame button
    
        check_ctrl_sample_file_pair_function() # Second step of control samples error checking: Check to see if paired-end dataset samples are properly paired

    else: # If dataset sequencing mode is not single-end nor paired-end
        return # Exit the function


# Second step of ChIP samples error checking: Check to see if paired-end dataset samples are properly paired
def check_chip_sample_file_pair_function():

    if read_mode_string_var.get() == 'paired': # Perform checks if dataset sequencing mode is paired-end
    
        chip_pair_error_state = 0 # Zero the ChIP pairing error detector variable

        if bool(chip_rep1_r1_string_var.get()) != bool(chip_rep1_r2_string_var.get()): # If only one of either the ChIP samples R1 or R2 is assigned to this replicate
            chip_pair_error_state = 1 # Set the ChIP pairing error to 1
            
            if not bool(chip_rep1_r1_string_var.get()): # Check is R1 is still unassigned
                chip_rep1_r1_entry.config(bg = 'IndianRed1') # Highlight where the problem is
            
            if not bool(chip_rep1_r2_string_var.get()): # Check is R2 is still unassigned
                chip_rep1_r2_entry.config(bg = 'IndianRed1') # Highlight where the problem is

        if bool(chip_rep2_r1_string_var.get()) != bool(chip_rep2_r2_string_var.get()): # If only one of either the ChIP samples R1 or R2 is assigned to this replicate
            chip_pair_error_state = 1 # Set the ChIP pairing error to 1

            if not bool(chip_rep2_r1_string_var.get()): # Check is R1 is still unassigned
                chip_rep2_r1_entry.config(bg = 'IndianRed1') # Highlight where the problem is
            
            if not bool(chip_rep2_r2_string_var.get()): # Check is R2 is still unassigned
                chip_rep2_r2_entry.config(bg = 'IndianRed1') # Highlight where the problem is

        if bool(chip_rep3_r1_string_var.get()) != bool(chip_rep3_r2_string_var.get()): # If only one of either the ChIP samples R1 or R2 is assigned to this replicate
            chip_pair_error_state = 1 # Set the ChIP pairing error to 1

            if not bool(chip_rep3_r1_string_var.get()): # Check is R1 is still unassigned
                chip_rep3_r1_entry.config(bg = 'IndianRed1') # Highlight where the problem is
            
            if not bool(chip_rep3_r2_string_var.get()): # Check is R2 is still unassigned
                chip_rep3_r2_entry.config(bg = 'IndianRed1') # Highlight where the problem is

        if bool(chip_rep4_r1_string_var.get()) != bool(chip_rep4_r2_string_var.get()): # If only one of either the ChIP samples R1 or R2 is assigned to this replicate
            chip_pair_error_state = 1 # Set the ChIP pairing error to 1

            if not bool(chip_rep4_r1_string_var.get()): # Check is R1 is still unassigned
                chip_rep4_r1_entry.config(bg = 'IndianRed1') # Highlight where the problem is
            
            if not bool(chip_rep4_r2_string_var.get()): # Check is R2 is still unassigned
                chip_rep4_r2_entry.config(bg = 'IndianRed1') # Highlight where the problem is

        if bool(chip_rep5_r1_string_var.get()) != bool(chip_rep5_r2_string_var.get()): # If only one of either the ChIP samples R1 or R2 is assigned to this replicate
            chip_pair_error_state = 1 # Set the ChIP pairing error to 1

            if not bool(chip_rep5_r1_string_var.get()): # Check is R1 is still unassigned
                chip_rep5_r1_entry.config(bg = 'IndianRed1') # Highlight where the problem is
            
            if not bool(chip_rep5_r2_string_var.get()): # Check is R2 is still unassigned
                chip_rep5_r2_entry.config(bg = 'IndianRed1') # Highlight where the problem is
    
        if chip_pair_error_state == 1: # If any sample pairing error is detected above 
            chip_sample_notification_string_var.set('One or more of ChIP samples have missing pair') # Notify the user
            chip_sample_notification_label.config(fg = 'red')
            chip_list_frame_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
            return # Exit the function

        elif chip_pair_error_state == 0: # If no sample pairing error is detected above 
            chip_sample_notification_string_var.set('All ChIP samples are properly paired') # Notify the user
            chip_sample_notification_label.config(fg = 'green')
            # Continue to next frame button is still enabled here as the result of previous function
        
        check_chip_sample_validity_function() # Third step of ChIP samples error checking: Check to see if samples have valid file extensions

    else: # If dataset sequencing mode is single-end, skip the proper-pairing checks above
        check_chip_sample_validity_function() # Third step of ChIP samples error checking: Check to see if samples have valid file extensions


# Second step of control samples error checking: Check to see if paired-end dataset samples are properly paired
def check_ctrl_sample_file_pair_function():

    if read_mode_string_var.get() == 'paired': # Perform checks if dataset sequencing mode is paired-end
    
        ctrl_pair_error_state = 0  # Zero the control pairing error detector variable

        if bool(ctrl_rep1_r1_string_var.get()) != bool(ctrl_rep1_r2_string_var.get()): # If only one of either the ChIP samples R1 or R2 is assigned to this replicate
            ctrl_pair_error_state = 1  # Set the ChIP pairing error to 1
            
            if not bool(ctrl_rep1_r1_string_var.get()): # Check is R1 is still unassigned
                ctrl_rep1_r1_entry.config(bg = 'IndianRed1') # Highlight where the problem is
            
            if not bool(ctrl_rep1_r2_string_var.get()): # Check is R2 is still unassigned
                ctrl_rep1_r2_entry.config(bg = 'IndianRed1') # Highlight where the problem is

        if bool(ctrl_rep2_r1_string_var.get()) != bool(ctrl_rep2_r2_string_var.get()): # If only one of either the ChIP samples R1 or R2 is assigned to this replicate
            ctrl_pair_error_state = 1 # Set the ChIP pairing error to 1

            if not bool(ctrl_rep2_r1_string_var.get()): # Check is R1 is still unassigned
                ctrl_rep2_r1_entry.config(bg = 'IndianRed1') # Highlight where the problem is
            
            if not bool(ctrl_rep2_r2_string_var.get()): # Check is R2 is still unassigned
                ctrl_rep2_r2_entry.config(bg = 'IndianRed1') # Highlight where the problem is

        if bool(ctrl_rep3_r1_string_var.get()) != bool(ctrl_rep3_r2_string_var.get()): # If only one of either the ChIP samples R1 or R2 is assigned to this replicate
            ctrl_pair_error_state = 1 # Set the ChIP pairing error to 1

            if not bool(ctrl_rep3_r1_string_var.get()): # Check is R1 is still unassigned
                ctrl_rep3_r1_entry.config(bg = 'IndianRed1') # Highlight where the problem is
            
            if not bool(ctrl_rep3_r2_string_var.get()): # Check is R2 is still unassigned
                ctrl_rep3_r2_entry.config(bg = 'IndianRed1') # Highlight where the problem is

        if bool(ctrl_rep4_r1_string_var.get()) != bool(ctrl_rep4_r2_string_var.get()): # If only one of either the ChIP samples R1 or R2 is assigned to this replicate
            ctrl_pair_error_state = 1 # Set the ChIP pairing error to 1

            if not bool(ctrl_rep4_r1_string_var.get()): # Check is R1 is still unassigned
                ctrl_rep4_r1_entry.config(bg = 'IndianRed1') # Highlight where the problem is
            
            if not bool(ctrl_rep4_r2_string_var.get()): # Check is R2 is still unassigned
                ctrl_rep4_r2_entry.config(bg = 'IndianRed1') # Highlight where the problem is

        if bool(ctrl_rep5_r1_string_var.get()) != bool(ctrl_rep5_r2_string_var.get()): # If only one of either the ChIP samples R1 or R2 is assigned to this replicate
            ctrl_pair_error_state = 1 # Set the ChIP pairing error to 1

            if not bool(ctrl_rep5_r1_string_var.get()): # Check is R1 is still unassigned
                ctrl_rep5_r1_entry.config(bg = 'IndianRed1') # Highlight where the problem is
            
            if not bool(ctrl_rep5_r2_string_var.get()): # Check is R2 is still unassigned
                ctrl_rep5_r2_entry.config(bg = 'IndianRed1') # Highlight where the problem is
    
        if ctrl_pair_error_state == 1: # If any sample pairing error is detected above 
            ctrl_sample_notification_string_var.set('One or more of control samples have missing pair') # Notify the user
            ctrl_sample_notification_label.config(fg = 'red')
            ctrl_list_frame_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
            return # Exit the function

        elif ctrl_pair_error_state == 0: # If no sample pairing error is detected above 
            ctrl_sample_notification_string_var.set('All control samples are properly paired') # Notify the user
            ctrl_sample_notification_label.config(fg = 'green')
            # Continue to next frame button is still enabled here as the result of previous function

        check_ctrl_sample_validity_function() # Third step of control samples error checking: Check to see if samples have valid file extensions

    else: # If dataset sequencing mode is single-end, skip the proper-pairing checks above
        check_ctrl_sample_validity_function() # Third step of control samples error checking: Check to see if samples have valid file extensions


# Third step of ChIP samples error checking: Check to see if samples have valid file extensions
def check_chip_sample_validity_function():

    chip_sample_format_error_state = 0 # Zero the ChIP file extension error detector variable

    if read_mode_string_var.get() == 'single': # In case of single-end dataset
        
        # If the sample extension is not among the accepted extensions list
        if bool(chip_rep1_r1_string_var.get()) and not any(chip_rep1_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep1_r1_entry.config(fg = 'red') # Highlight where the problem is
            chip_sample_format_error_state = 1 # Set the ChIP file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(chip_rep2_r1_string_var.get()) and not any(chip_rep2_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep2_r1_entry.config(fg = 'red') # Highlight where the problem is
            chip_sample_format_error_state = 1 # Set the ChIP file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(chip_rep3_r1_string_var.get()) and not any(chip_rep3_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep3_r1_entry.config(fg = 'red') # Highlight where the problem is
            chip_sample_format_error_state = 1 # Set the ChIP file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(chip_rep4_r1_string_var.get()) and not any(chip_rep4_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep4_r1_entry.config(fg = 'red') # Highlight where the problem is
            chip_sample_format_error_state = 1 # Set the ChIP file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(chip_rep5_r1_string_var.get()) and not any(chip_rep5_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep5_r1_entry.config(fg = 'red') # Highlight where the problem is
            chip_sample_format_error_state = 1 # Set the ChIP file extension error to 1


    if read_mode_string_var.get() == 'paired': # In case of paired-end dataset

        # If the sample extension is not among the accepted extensions list
        if bool(chip_rep1_r1_string_var.get()) and not any(chip_rep1_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep1_r1_entry.config(fg = 'red') # Highlight where the problem is
            chip_sample_format_error_state = 1 # Set the ChIP file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(chip_rep1_r2_string_var.get()) and not any(chip_rep1_r2_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep1_r2_entry.config(fg = 'red') # Highlight where the problem is
            chip_sample_format_error_state = 1 # Set the ChIP file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(chip_rep2_r1_string_var.get()) and not any(chip_rep2_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep2_r1_entry.config(fg = 'red') # Highlight where the problem is
            chip_sample_format_error_state = 1 # Set the ChIP file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(chip_rep2_r2_string_var.get()) and not any(chip_rep2_r2_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep2_r2_entry.config(fg = 'red') # Highlight where the problem is
            chip_sample_format_error_state = 1 # Set the ChIP file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(chip_rep3_r1_string_var.get()) and not any(chip_rep3_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep3_r1_entry.config(fg = 'red') # Highlight where the problem is
            chip_sample_format_error_state = 1 # Set the ChIP file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(chip_rep3_r2_string_var.get()) and not any(chip_rep3_r2_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep3_r2_entry.config(fg = 'red') # Highlight where the problem is
            chip_sample_format_error_state = 1 # Set the ChIP file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(chip_rep4_r1_string_var.get()) and not any(chip_rep4_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep4_r1_entry.config(fg = 'red') # Highlight where the problem is
            chip_sample_format_error_state = 1 # Set the ChIP file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(chip_rep4_r2_string_var.get()) and not any(chip_rep4_r2_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep4_r2_entry.config(fg = 'red') # Highlight where the problem is
            chip_sample_format_error_state = 1 # Set the ChIP file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(chip_rep5_r1_string_var.get()) and not any(chip_rep5_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep5_r1_entry.config(fg = 'red') # Highlight where the problem is
            chip_sample_format_error_state = 1 # Set the ChIP file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(chip_rep5_r2_string_var.get()) and not any(chip_rep5_r2_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep5_r2_entry.config(fg = 'red') # Highlight where the problem is
            chip_sample_format_error_state = 1 # Set the ChIP file extension error to 1


    if chip_sample_format_error_state == 1: # If any sample file extension error is detected above
        chip_sample_notification_string_var.set('One or more ChIP samples do not have a valid file extension') # Notify the user
        chip_sample_notification_label.config(fg = 'red')
        chip_list_frame_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
        return # Exit the function

    elif chip_sample_format_error_state == 0: # If no sample file extension error is detected above
        chip_sample_notification_string_var.set('All ChIP samples have valid file extension') # Notify the user
        chip_sample_notification_label.config(fg = 'green')
        # Continue to next frame button is still enabled here as the result of previous function

    check_chip_sample_file_exist_function() # Fourth step of control samples error checking: Check to see if assigned sample files actually exist


# Third step of control samples error checking: Check to see if samples have valid file extensions
def check_ctrl_sample_validity_function():

    ctrl_sample_format_error_state = 0 # Zero the control file extension error detector variable

    if read_mode_string_var.get() == 'single': # In case of single-end dataset
        
        # If the sample extension is not among the accepted extensions list
        if bool(ctrl_rep1_r1_string_var.get()) and not any(ctrl_rep1_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep1_r1_entry.config(fg = 'red') # Highlight where the problem is
            ctrl_sample_format_error_state = 1 # Set the control file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(ctrl_rep2_r1_string_var.get()) and not any(ctrl_rep2_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep2_r1_entry.config(fg = 'red') # Highlight where the problem is
            ctrl_sample_format_error_state = 1 # Set the control file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(ctrl_rep3_r1_string_var.get()) and not any(ctrl_rep3_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep3_r1_entry.config(fg = 'red') # Highlight where the problem is
            ctrl_sample_format_error_state = 1 # Set the control file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(ctrl_rep4_r1_string_var.get()) and not any(ctrl_rep4_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep4_r1_entry.config(fg = 'red') # Highlight where the problem is
            ctrl_sample_format_error_state = 1 # Set the control file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(ctrl_rep5_r1_string_var.get()) and not any(ctrl_rep5_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep5_r1_entry.config(fg = 'red') # Highlight where the problem is
            ctrl_sample_format_error_state = 1 # Set the control file extension error to 1


    if read_mode_string_var.get() == 'paired': # In case of paired-end dataset

        # If the sample extension is not among the accepted extensions list
        if bool(ctrl_rep1_r1_string_var.get()) and not any(ctrl_rep1_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep1_r1_entry.config(fg = 'red') # Highlight where the problem is
            ctrl_sample_format_error_state = 1 # Set the control file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(ctrl_rep1_r2_string_var.get()) and not any(ctrl_rep1_r2_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep1_r2_entry.config(fg = 'red') # Highlight where the problem is
            ctrl_sample_format_error_state = 1 # Set the control file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(ctrl_rep2_r1_string_var.get()) and not any(ctrl_rep2_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep2_r1_entry.config(fg = 'red') # Highlight where the problem is
            ctrl_sample_format_error_state = 1 # Set the control file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(ctrl_rep2_r2_string_var.get()) and not any(ctrl_rep2_r2_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep2_r2_entry.config(fg = 'red') # Highlight where the problem is
            ctrl_sample_format_error_state = 1 # Set the control file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(ctrl_rep3_r1_string_var.get()) and not any(ctrl_rep3_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep3_r1_entry.config(fg = 'red') # Highlight where the problem is
            ctrl_sample_format_error_state = 1 # Set the control file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(ctrl_rep3_r2_string_var.get()) and not any(ctrl_rep3_r2_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep3_r2_entry.config(fg = 'red') # Highlight where the problem is
            ctrl_sample_format_error_state = 1 # Set the control file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(ctrl_rep4_r1_string_var.get()) and not any(ctrl_rep4_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep4_r1_entry.config(fg = 'red') # Highlight where the problem is
            ctrl_sample_format_error_state = 1 # Set the control file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(ctrl_rep4_r2_string_var.get()) and not any(ctrl_rep4_r2_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep4_r2_entry.config(fg = 'red') # Highlight where the problem is
            ctrl_sample_format_error_state = 1 # Set the control file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(ctrl_rep5_r1_string_var.get()) and not any(ctrl_rep5_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep5_r1_entry.config(fg = 'red') # Highlight where the problem is
            ctrl_sample_format_error_state = 1 # Set the control file extension error to 1

        # If the sample extension is not among the accepted extensions list
        if bool(ctrl_rep5_r2_string_var.get()) and not any(ctrl_rep5_r2_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep5_r2_entry.config(fg = 'red') # Highlight where the problem is
            ctrl_sample_format_error_state = 1 # Set the control file extension error to 1


    if ctrl_sample_format_error_state == 1: # If any sample file extension error is detected above
        ctrl_sample_notification_string_var.set('One or more control samples do not have a valid file extension') # Notify the user
        ctrl_sample_notification_label.config(fg = 'red')
        ctrl_list_frame_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
        return # Exit the function

    elif ctrl_sample_format_error_state == 0: # If no sample file extension error is detected above
        ctrl_sample_notification_string_var.set('All control samples have valid file extension') # Notify the user
        ctrl_sample_notification_label.config(fg = 'green')
        # Continue to next frame button is still enabled here as the result of previous function

    check_ctrl_sample_file_exist_function() # Fourth step of control samples error checking: Check to see if assigned sample files actually exist


# Fourth step of ChIP samples error checking: Check to see if assigned sample files actually exist
def check_chip_sample_file_exist_function():

    chip_file_exist_error_state = 0 # Zero the ChIP file not found error detector variable

    if chip_rep1_r1_string_var.get(): # Perform the check only if the sample variable value is not empty
        if not os.path.isfile(chip_rep1_r1_string_var.get()): # If the entered path leads to a non-existent file
            chip_file_exist_error_state = 1 # Set ChIP file not found error to 1
            chip_rep1_r1_entry.config(fg = 'red') # Highlight where the problem is

    if chip_rep1_r2_string_var.get(): # Perform the check only if the sample variable value is not empty
        if not os.path.isfile(chip_rep1_r2_string_var.get()): # If the entered path leads to a non-existent file
            chip_file_exist_error_state = 1 # Set ChIP file not found error to 1
            chip_rep1_r2_entry.config(fg = 'red') # Highlight where the problem is

    if chip_rep2_r1_string_var.get(): # Perform the check only if the sample variable value is not empty
        if not os.path.isfile(chip_rep2_r1_string_var.get()): # If the entered path leads to a non-existent file
            chip_file_exist_error_state = 1 # Set ChIP file not found error to 1
            chip_rep2_r1_entry.config(fg = 'red') # Highlight where the problem is

    if chip_rep2_r2_string_var.get(): # Perform the check only if the sample variable value is not empty
        if not os.path.isfile(chip_rep2_r2_string_var.get()): # If the entered path leads to a non-existent file
            chip_file_exist_error_state = 1 # Set ChIP file not found error to 1
            chip_rep2_r2_entry.config(fg = 'red') # Highlight where the problem is

    if chip_rep3_r1_string_var.get(): # Perform the check only if the sample variable value is not empty
        if not os.path.isfile(chip_rep3_r1_string_var.get()): # If the entered path leads to a non-existent file
            chip_file_exist_error_state = 1 # Set ChIP file not found error to 1
            chip_rep3_r1_entry.config(fg = 'red') # Highlight where the problem is

    if chip_rep3_r2_string_var.get(): # Perform the check only if the sample variable value is not empty
        if not os.path.isfile(chip_rep3_r2_string_var.get()): # If the entered path leads to a non-existent file
            chip_file_exist_error_state = 1 # Set ChIP file not found error to 1
            chip_rep3_r2_entry.config(fg = 'red') # Highlight where the problem is

    if chip_rep4_r1_string_var.get(): # Perform the check only if the sample variable value is not empty
        if not os.path.isfile(chip_rep4_r1_string_var.get()): # If the entered path leads to a non-existent file
            chip_file_exist_error_state = 1 # Set ChIP file not found error to 1
            chip_rep4_r1_entry.config(fg = 'red') # Highlight where the problem is

    if chip_rep4_r2_string_var.get(): # Perform the check only if the sample variable value is not empty
        if not os.path.isfile(chip_rep4_r2_string_var.get()): # If the entered path leads to a non-existent file
            chip_file_exist_error_state = 1 # Set ChIP file not found error to 1
            chip_rep4_r2_entry.config(fg = 'red') # Highlight where the problem is

    if chip_rep5_r1_string_var.get(): # Perform the check only if the sample variable value is not empty
        if not os.path.isfile(chip_rep5_r1_string_var.get()): # If the entered path leads to a non-existent file
            chip_file_exist_error_state = 1 # Set ChIP file not found error to 1
            chip_rep5_r1_entry.config(fg = 'red') # Highlight where the problem is

    if chip_rep5_r2_string_var.get(): # Perform the check only if the sample variable value is not empty
        if not os.path.isfile(chip_rep5_r2_string_var.get()): # If the entered path leads to a non-existent file
            chip_file_exist_error_state = 1 # Set ChIP file not found error to 1
            chip_rep5_r2_entry.config(fg = 'red') # Highlight where the problem is

    if chip_file_exist_error_state == 1: # If any sample file not found error is detected above
        chip_sample_notification_string_var.set('One or more ChIP sample files do not exist') # Notify the user
        chip_sample_notification_label.config(fg = 'red')
        chip_list_frame_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button

    elif chip_file_exist_error_state == 0: # If no sample file not found error is detected above
        chip_sample_notification_string_var.set('No problem found in ChIP samples') # Notify the user
        chip_sample_notification_label.config(fg = 'green')
        chip_list_frame_continue_button.config(state = tk.NORMAL) # Enable the continue to next frame button


# Fourth step of control samples error checking: Check to see if assigned sample files actually exist
def check_ctrl_sample_file_exist_function():
    
    ctrl_file_exist_error_state = 0 # Zero the control file not found error detector variable
    
    if ctrl_rep1_r1_string_var.get(): # Perform the check only if the sample variable value is not empty
        if not os.path.isfile(ctrl_rep1_r1_string_var.get()): # If the entered path leads to a non-existent file
            ctrl_file_exist_error_state = 1 # Set control file not found error to 1
            ctrl_rep1_r1_entry.config(fg = 'red') # Highlight where the problem is

    if ctrl_rep1_r2_string_var.get(): # Perform the check only if the sample variable value is not empty
        if not os.path.isfile(ctrl_rep1_r2_string_var.get()): # If the entered path leads to a non-existent file
            ctrl_file_exist_error_state = 1 # Set control file not found error to 1
            ctrl_rep1_r2_entry.config(fg = 'red') # Highlight where the problem is

    if ctrl_rep2_r1_string_var.get(): # Perform the check only if the sample variable value is not empty
        if not os.path.isfile(ctrl_rep2_r1_string_var.get()): # If the entered path leads to a non-existent file
            ctrl_file_exist_error_state = 1 # Set control file not found error to 1
            ctrl_rep2_r1_entry.config(fg = 'red') # Highlight where the problem is

    if ctrl_rep2_r2_string_var.get(): # Perform the check only if the sample variable value is not empty
        if not os.path.isfile(ctrl_rep2_r2_string_var.get()): # If the entered path leads to a non-existent file
            ctrl_file_exist_error_state = 1 # Set control file not found error to 1
            ctrl_rep2_r2_entry.config(fg = 'red') # Highlight where the problem is

    if ctrl_rep3_r1_string_var.get(): # Perform the check only if the sample variable value is not empty
        if not os.path.isfile(ctrl_rep3_r1_string_var.get()): # If the entered path leads to a non-existent file
            ctrl_file_exist_error_state = 1 # Set control file not found error to 1
            ctrl_rep3_r1_entry.config(fg = 'red') # Highlight where the problem is

    if ctrl_rep3_r2_string_var.get(): # Perform the check only if the sample variable value is not empty
        if not os.path.isfile(ctrl_rep3_r2_string_var.get()): # If the entered path leads to a non-existent file
            ctrl_file_exist_error_state = 1 # Set control file not found error to 1
            ctrl_rep3_r2_entry.config(fg = 'red') # Highlight where the problem is

    if ctrl_rep4_r1_string_var.get(): # Perform the check only if the sample variable value is not empty
        if not os.path.isfile(ctrl_rep4_r1_string_var.get()): # If the entered path leads to a non-existent file
            ctrl_file_exist_error_state = 1 # Set control file not found error to 1
            ctrl_rep4_r1_entry.config(fg = 'red') # Highlight where the problem is

    if ctrl_rep4_r2_string_var.get(): # Perform the check only if the sample variable value is not empty
        if not os.path.isfile(ctrl_rep4_r2_string_var.get()): # If the entered path leads to a non-existent file
            ctrl_file_exist_error_state = 1 # Set control file not found error to 1
            ctrl_rep4_r2_entry.config(fg = 'red') # Highlight where the problem is

    if ctrl_rep5_r1_string_var.get(): # Perform the check only if the sample variable value is not empty
        if not os.path.isfile(ctrl_rep5_r1_string_var.get()): # If the entered path leads to a non-existent file
            ctrl_file_exist_error_state = 1 # Set control file not found error to 1
            ctrl_rep5_r1_entry.config(fg = 'red') # Highlight where the problem is

    if ctrl_rep5_r2_string_var.get(): # Perform the check only if the sample variable value is not empty
        if not os.path.isfile(ctrl_rep5_r2_string_var.get()): # If the entered path leads to a non-existent file
            ctrl_file_exist_error_state = 1 # Set control file not found error to 1
            ctrl_rep5_r2_entry.config(fg = 'red') # Highlight where the problem is

    if ctrl_file_exist_error_state == 1: # If any sample file not found error is detected above
        ctrl_sample_notification_string_var.set('One or more control sample files do not exist') # Notify the user
        ctrl_sample_notification_label.config(fg = 'red')
        ctrl_list_frame_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button

    elif ctrl_file_exist_error_state == 0: # If no sample file not found error is detected above
        ctrl_sample_notification_string_var.set('No problem found in control samples') # Notify the user
        ctrl_sample_notification_label.config(fg = 'green')
        ctrl_list_frame_continue_button.config(state = tk.NORMAL) # Enable the continue to next frame button


# Check if the required inputs other than sample errors are satisifed or contain any errors, before giving user access to proceed
def check_required_input_function(*args):
    
    if read_mode_string_var.get() != 'single' and read_mode_string_var.get() != 'paired': # If inputted sequencing mode is none of the accepted choices
        read_mode_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
    else: # If inputted sequencing mode is within the accepted choices
        read_mode_continue_button.config(state = tk.NORMAL) # Enable the continue to next frame button

    if peak_type_string_var.get() != 'narrow' and peak_type_string_var.get() != 'broad' and peak_type_string_var.get() != 'unsure': # If inputted peak type is none of the accepted choices
        peak_type_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
    else: # If inputted peak type is within the accepted choices
        peak_type_continue_button.config(state = tk.NORMAL) # Enable the continue to next frame button


    if not bool(sample_table_question_string_var.get()): # If user has not decided whether to use sample table
        sample_table_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
    elif bool(sample_table_question_string_var.get()): # If user has decided whether to use sample table
        if sample_table_question_string_var.get() == 'yes': # If user has not decided to use sample table
            if not bool(sample_table_string_var.get()): # If the sample table variable value is empty
                sample_table_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
            elif not os.path.isfile(sample_table_string_var.get()): # If the entered path leads to a non-existent file
                sample_table_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
            elif os.path.isfile(sample_table_string_var.get()): # If the entered path leads to an existing file
                sample_table_continue_button.config(state = tk.NORMAL) # Enable the continue to next frame button
        if sample_table_question_string_var.get() == 'no': # If user has decided to not use sample table 
            sample_table_continue_button.config(state = tk.NORMAL) # Enable the continue to next frame button


    if not bool(setting_table_question_string_var.get()): # If user has not decided whether to use setting table
        setting_table_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
    else:
        if setting_table_question_string_var.get() == 'yes': # If user has decided whether to use sample table
            if not bool(setting_table_string_var.get()): # If the setting table variable value is empty
                setting_table_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
                setting_table_customize_button.config(state = tk.DISABLED) # Disable the go to customize settings frame button
            if '[MODIFIED]' in setting_table_string_var.get(): # If the setting table variable value contains an indicator of setting values manual modification
                pass # Do nothing
            elif not os.path.isfile(setting_table_string_var.get()): # If the entered path leads to a non-existent directory
                setting_table_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
                setting_table_customize_button.config(state = tk.DISABLED) # Disable the go to customize settings frame button
            elif os.path.isfile(setting_table_string_var.get()): # If the entered path leads to an existing directory
                setting_table_continue_button.config(state = tk.NORMAL) # Enable the continue to next frame button
                setting_table_customize_button.config(state = tk.NORMAL) # Enable the go to customize settings frame button
        if setting_table_question_string_var.get() == 'no': # If user has decided to not use setting table     
            setting_table_continue_button.config(state = tk.NORMAL) # Enable the continue to next frame button
            setting_table_customize_button.config(state = tk.NORMAL) # Enable the go to customize settings frame button


    if not bool(genome_ref_string_var.get()): # If the reference genome build variable value is empty
        genome_ref_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button

    # If the reference genome build variable value is not empty
    elif not bool(genome_folder_string_var.get()): # If the reference genome directory variable value is empty
        genome_ref_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
    elif not os.path.isdir(genome_folder_string_var.get()): # If the entered path leads to a non-existent directory
        genome_ref_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
        genome_folder_entry.config(fg = 'red') # Highlight where the problem is
    elif os.path.isdir(genome_folder_string_var.get()): # If the entered path leads to an existing directory
        genome_ref_continue_button.config(state = tk.NORMAL) # Enable the continue to next frame button


    if not bool(known_motif_question_string_var.get()): # If user has not decided whether to use a known motif file
        known_motif_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
        known_motif_entry.config(fg = 'black') # Remove the highlight, if previously highlighted
    elif known_motif_question_string_var.get() == 'no': # If user has decided to not use a known motif file
        known_motif_continue_button.config(state = tk.NORMAL) # Enable the continue to next frame button
        known_motif_entry.config(fg = 'black') # Remove the highlight, if previously highlighted
    elif known_motif_question_string_var.get() == 'yes': # If user has decided to use a known motif file 
        if known_motif_string_var.get() == '': # If the known motif variable value is empty
            known_motif_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
            known_motif_entry.config(fg = 'red') # Highlight where the problem is
        elif known_motif_string_var.get() != '': # If the known motif variable value is not empty
            if not os.path.isfile(known_motif_string_var.get()): # If the entered path leads to a non-existent directory
                known_motif_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
                known_motif_entry.config(fg = 'red') # Highlight where the problem is
            if os.path.isfile(known_motif_string_var.get()): # If the entered path leads to an existing directory
                known_motif_continue_button.config(state = tk.NORMAL) # Enable the continue to next frame button
                known_motif_entry.config(fg = 'black') # Remove the highlight, if previously highlighted


    if not bool(setname_string_var.get()): # If the dataset name variable value is empty
        output_folder_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
    elif not bool(output_folder_string_var.get()): # If the output directory variable value is empty
        output_folder_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button     
    elif bool(output_folder_string_var.get()): # If the output directory variable value is not empty
        output_folder_continue_button.config(state = tk.NORMAL) # Enable the continue to next frame button


    if not bool(homer_motif_question_string_var.get()): # If user has not decided whether to perform HOMER motif enrichment analysis
        homer_motif_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
    elif homer_motif_question_string_var.get() == 'no': # If user has decided to not perform HOMER motif enrichment analysis
        homer_motif_continue_button.config(state = tk.NORMAL) # Enable the continue to next frame button
    elif homer_motif_question_string_var.get() == 'yes': # If user has decided to perform HOMER motif enrichment analysis
        if homer_motif_string_var.get() == "<Please select peak set(s)>": # If user has not decided which peak set to perform HOMER motif enrichment analysis on
            homer_motif_drop_down.config(fg = 'blue')
            homer_motif_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
        elif homer_motif_string_var.get() != "<Please select peak set(s)>": # If user has decided which peak set to perform HOMER motif enrichment analysis on
            homer_motif_drop_down.config(fg = 'black') 
            homer_motif_continue_button.config(state = tk.NORMAL) # Enable the continue to next frame button


    if not bool(meme_motif_question_string_var.get()): # If user has not decided whether to perform MEME motif enrichment analysis
        meme_motif_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
    elif meme_motif_question_string_var.get() == 'no': # If user has decided to not perform MEME motif enrichment analysis
        meme_motif_continue_button.config(state = tk.NORMAL) # Enable the continue to next frame button
    elif meme_motif_question_string_var.get() == 'yes': # If user has decided to perform MEME motif enrichment analysis
        if meme_motif_string_var.get() == "<Please select peak set(s)>": # If user has not decided which peak set to perform MEME motif enrichment analysis on
            meme_motif_drop_down.config(fg = 'blue')
            meme_motif_continue_button.config(state = tk.DISABLED) # Disable the continue to next frame button
        elif meme_motif_string_var.get() != "<Please select peak set(s)>": # If user has decided which peak set to perform MEME motif enrichment analysis on
            meme_motif_drop_down.config(fg = 'black')
            meme_motif_continue_button.config(state = tk.NORMAL) # Enable the continue to next frame button

 
    if not bool(cpu_count_string_var.get()): # If the CPU count variable value is empty
        generate_scripts_button.config(state = tk.DISABLED)
        generate_and_run_scripts_button.config(state = tk.DISABLED)
        cpu_count_entry.config(fg = 'black') # Set font color for CPU count entry field prior to any input
        cpu_count_notification_string_var.set("Maximum number of CPU cores available: {}".format(max_cpu)) # Inform the user of the acceptable range of values
        cpu_count_notification_label.config(fg = 'blue')
    
    elif bool(cpu_count_string_var.get()): # If the CPU count variable value is not empty
        try: # If the entered value can be converted into python integer
            input_cpu_count = int(cpu_count_string_var.get())

            if int(cpu_count_string_var.get()) > max_cpu: # If the CPU count variable value is larger than maximum capacity
                generate_scripts_button.config(state = tk.DISABLED) # Disable the big red button
                generate_and_run_scripts_button.config(state = tk.DISABLED) # Disable the big red button
                cpu_count_notification_string_var.set("Entered number exceeds available CPU cores ({})".format(max_cpu)) # Notify the user
                cpu_count_notification_label.config(fg = 'red')

            elif int(cpu_count_string_var.get()) < 1: # If the CPU count variable value is lower than minimum needed
                generate_scripts_button.config(state = tk.DISABLED) # Disable the big red button
                generate_and_run_scripts_button.config(state = tk.DISABLED) # Disable the big red button
                cpu_count_notification_string_var.set("Need at least one CPU core to run the pipeline") # Notify the user
                cpu_count_notification_label.config(fg = 'red')

            elif int(cpu_count_string_var.get()) <= max_cpu: # If the CPU count variable value is within acceptable range of values
                generate_scripts_button.config(state = tk.NORMAL) # Enable the big red button
                generate_and_run_scripts_button.config(state = tk.NORMAL) # Enable the big red button
                cpu_count_notification_string_var.set("ChIP-AP ready for action!") # Notify the user
                cpu_count_notification_label.config(fg = 'green')
        
        except: # If the entered value cannot be converted into python integer
            generate_scripts_button.config(state = tk.DISABLED) # Disable the big red button
            generate_and_run_scripts_button.config(state = tk.DISABLED) # Disable the big red button            
            cpu_count_notification_string_var.set("Entered value is not an integer") # Notify the user
            cpu_count_notification_label.config(fg = 'red')
    

# Update the generated ChIP-AP command line based on the latest registered variable values
def update_command_line_function(*args):

    if read_mode_string_var.get() == 'single' or read_mode_string_var.get() == 'paired': # If inputted sequencing mode is within the accepted choices
        read_mode_arg.set(' --mode {}'.format(read_mode_string_var.get())) # Assign the sequencing mode argument for the command line behind the corresponding flag
    else: # If inputted sequencing mode is none of the accepted choices
        read_mode_arg.set('') # Empty the sequencing mode argument for the command line

    if peak_type_string_var.get() == 'narrow' or peak_type_string_var.get() == 'broad' or peak_type_string_var.get() == 'unsure': # If inputted peak type is within the accepted choices
        peak_type_arg.set(' --peak {}'.format(peak_type_string_var.get())) # Assign the peak type argument for the command line behind the corresponding flag
    else: # If inputted peak type is none of the accepted choices
        peak_type_arg.set('') # Empty the peak type argument for the command line

    if bool(output_folder_string_var.get()): # If the output directory variable value is not empty
        output_folder_arg.set(' --output {}'.format(output_folder_string_var.get())) # Assign the output_directory argument for the command line behind the corresponding flag
    if not bool(output_folder_string_var.get()): # If the output directory variable value is empty
        output_folder_arg.set('') # Empty the output directory argument for the command line

    if bool(setname_string_var.get()): # If the dataset name variable value is not empty
        setname_arg.set(' --setname {}'.format(setname_string_var.get())) # Assign the dataset name argument for the command line behind the corresponding flag
    if not bool(setname_string_var.get()): # If the dataset name variable value is empty
        setname_arg.set('') # Empty the dataset name argument for the command line

    output_dir_arg.set('{}/{}'.format(os.path.abspath(output_folder_string_var.get()), setname_string_var.get())) # All results are actually saved in [output folder]/[dataset name]

    if bool(genome_ref_string_var.get()): # If the reference genome build variable value is not empty
        genome_ref_arg.set(' --ref {}'.format(genome_ref_string_var.get().split(' ')[0])) # Assign the reference genome build argument for the command line behind the corresponding flag
    if not bool(genome_ref_string_var.get()): # If the reference genome build variable value is empty
        genome_ref_arg.set('') # Empty the reference genome build argument for the command line

    if bool(genome_folder_string_var.get()): # If the genome directory variable value is not empty
        genome_folder_arg.set(' --genome {}'.format(genome_folder_string_var.get())) # Assign the genome directory argument for the command line behind the corresponding flag
    if not bool(genome_folder_string_var.get()): # If the genome directory variable value is empty
        genome_folder_arg.set('') # Empty the genome directory argument for the command line

    if bool(output_folder_string_var.get()) and bool(setname_string_var.get()): # If neither the output directory nor dataset name variable values are empty
        # Assign the sample table argument for the command line behind the corresponding flag
        sample_table_arg.set(' --sample_table {}/{}_sample_table.tsv'.format(output_dir_arg.get(), setname_string_var.get()))
    else: # If either the output directory or dataset name variable values are empty
        sample_table_arg.set('') # Empty the sample table argument for the command line

    if bool(output_folder_string_var.get()) and bool(setname_string_var.get()): # If neither the output directory nor dataset name variable values are empty
        # Assign the setting table argument for the command line behind the corresponding flag
        setting_table_arg.set(' --custom_setting_table {}/{}_setting_table.tsv'.format(output_dir_arg.get(), setname_string_var.get()))
    else: # If either the output directory or dataset name variable values are empty
        setting_table_arg.set('') # Empty the setting table argument for the command line

    if bool(known_motif_string_var.get()): # If the known motif variable value is not empty
        known_motif_arg.set(' --motif {}'.format(known_motif_string_var.get())) # Assign the known motif argument for the command line behind the corresponding flag
    if not bool(known_motif_string_var.get()): # If the known motif variable value is empty
        known_motif_arg.set('') # Empty the known motif argument for the command line

    # This argument does not apply and will be absent from dataset with broad peak type
    if peak_type_string_var.get() != 'broad': # If dataset peak type is not broad
        if homer_motif_string_var.get() != '<Please select peak set(s)>': # If user has decided which peak set to perform HOMER motif enrichment analysis on
            # Assign the HOMER motif enrichment analysis argument (consensus, union, or both) for the command line behind the corresponding flag
            homer_motif_arg.set(' --homer_motif {}'.format(homer_motif_string_var.get().split(' ')[0].lower()))
        elif homer_motif_string_var.get() == '<Please select peak set(s)>': # If user has not decided which peak set to perform HOMER motif enrichment analysis on
            homer_motif_arg.set('') # Empty the HOMER motif enrichment analysis argument for the command line
    elif peak_type_string_var.get() == 'broad': # If dataset peak type is broad
            homer_motif_arg.set('') # Empty the HOMER motif enrichment analysis argument for the command line

    # This argument does not apply and will be absent from dataset with broad peak type        
    if peak_type_string_var.get() != 'broad': # If dataset peak type is not broad
        if meme_motif_string_var.get() != '<Please select peak set(s)>': # If user has decided which peak set to perform MEME motif enrichment analysis on
            # Assign the MEME motif enrichment analysis argument (consensus, union, or both) for the command line behind the corresponding flag
            meme_motif_arg.set(' --meme_motif {}'.format(meme_motif_string_var.get().split(' ')[0].lower()))
        elif meme_motif_string_var.get() == '<Please select peak set(s)>': # If user has not decided which peak set to perform MEME motif enrichment analysis on
            meme_motif_arg.set('') # Empty the MEME motif enrichment analysis argument for the command line
    elif peak_type_string_var.get() == 'broad': # If dataset peak type is broad
            meme_motif_arg.set('') # Empty the MEME motif enrichment analysis argument for the command line

    if fcmerge_var.get() == 1: # If the force merge variable value is set to 1
        fcmerge_arg.set(' --fcmerge') # Assign the force merge flag for the command line
    if fcmerge_var.get() == 0: # If the force merge variable value is set to 0
        fcmerge_arg.set('') # Empty the force merge argument for the command line

    if goann_var.get() == 1: # If the GO annotation variable value is set to 1
        goann_arg.set(' --goann') # Assign the GO annotation flag for the command line
    if goann_var.get() == 0: # If the GO annotation variable value is set to 0
        goann_arg.set('') # Empty the GO annotation argument for the command line

    if pathann_var.get() == 1: # If the pathway annotation variable value is set to 1
        pathann_arg.set(' --pathann') # Assign the pathway annotation flag for the command line
    if pathann_var.get() == 0: # If the pathway annotation variable value is set to 0
        pathann_arg.set('') # Empty the pathway annotation argument for the command line

    if deltemp_var.get() == 1: # If the delete temporary files variable value is set to 1
        deltemp_arg.set(' --deltemp') # Assign the delete temporary files flag for the command line
    if deltemp_var.get() == 0: # If the delete temporary files variable value is set to 0
        deltemp_arg.set('')  # Empty the delete temporary files argument for the command line

    if bool(cpu_count_string_var.get()): # If the CPU count variable value is not empty
        cpu_count_arg.set(' --thread {}'.format(cpu_count_string_var.get())) # Assign the CPU count argument for the command line behind the corresponding flag
    if not bool(cpu_count_string_var.get()): # If the CPU count variable value is empty
        cpu_count_arg.set('') # Empty the CPU count argument for the command line

    # If neither the standard output/error, output directory, nor dataset name variable values are empty
    if stdout_var.get() == 1 and bool(output_folder_string_var.get()) and bool(setname_string_var.get()):
        stdout_arg.set(' 1> {}/{}.out'.format(output_dir_arg.get(), setname_string_var.get())) # Assign the standard output argument for the command line behind the corresponding flag
        stderr_arg.set(' 2> {}/{}.err'.format(output_dir_arg.get(), setname_string_var.get())) # Assign the standard error argument for the command line behind the corresponding flag
    else: # If either the standard output/error, output directory, or dataset name variable values are empty
        stdout_arg.set('') # Empty the standard output argument for the command line
        stderr_arg.set('') # Empty the standard error argument for the command line

    # Combine all arguments above into a full command line to call for whole pipeline, which begins with the ChIP-AP main script name (chipap.py)
    # Empty arguments will simply be absent from the command line without disrupting the command
    command_line_output_string_var.set('{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}'.format(chipap_program_name,
                                                                                        read_mode_arg.get(),
                                                                                        peak_type_arg.get(),
                                                                                        output_folder_arg.get(),
                                                                                        setname_arg.get(),
                                                                                        genome_ref_arg.get(),
                                                                                        genome_folder_arg.get(),
                                                                                        sample_table_arg.get(),
                                                                                        setting_table_arg.get(),
                                                                                        known_motif_arg.get(),
                                                                                        homer_motif_arg.get(),
                                                                                        meme_motif_arg.get(),
                                                                                        fcmerge_arg.get(),
                                                                                        goann_arg.get(),
                                                                                        pathann_arg.get(),
                                                                                        deltemp_arg.get(),
                                                                                        cpu_count_arg.get(),
                                                                                        stdout_arg.get(),
                                                                                        stderr_arg.get()))

# Function to be run on non-sample variables which are traced for any value change
def check_traced_input_function(*args):
    check_required_input_function() # Check if the required inputs other than sample errors are satisifed or contain any errors, before giving user access to proceed
    update_command_line_function() # Update the generated ChIP-AP command line based on the latest registered variable values

# Function to be run on particularly the traced sample table variable, enabling the values in the file to be loaded properly upon typing inside the sample table entry field
def sample_table_entry_trace_load_function(*args):
    if bool(sample_table_string_var.get()): # If the sample table variable value is not empty
        if not os.path.isfile(sample_table_string_var.get()): # If the entered path leads to a non-existent file
            sample_table_notification_string_var.set('Sample table not found') # Notify the user
            sample_table_notification_label.config(fg = 'red')
            clear_sample_function('valuesonly') # Clear all ChIP and control samples variable values, but not the sample table variable value
        elif os.path.isfile(sample_table_string_var.get()): # If the entered path leads to an existing file
            sample_table_loading_test_function() # Load the sample table normally, as in when using the file browser via button click

# Function to be run on particularly the traced setting table variable, enabling the values in the file to be loaded properly upon typing inside the setting table entry field
def setting_table_entry_trace_load_function(*args):
    if bool(setting_table_string_var.get()):
        if '[MODIFIED]' in setting_table_string_var.get(): # If the setting table variable value contains an indicator of setting values manual modification
            pass # Do nothing
        elif not os.path.isfile(setting_table_string_var.get()): # If the entered path leads to a non-existent file
            setting_table_notification_string_var.set('Setting table not found') # Notify the user
            setting_table_notification_label.config(fg = 'red')
            clear_setting_function('valuesonly') # Clear all program settings variable values, but not the setting table variable value
            read_setting_table_function(default_setting_table_file_full_path) # Load the default setting values from the provided default setting table
        elif os.path.isfile(setting_table_string_var.get()): # If the entered path leads to an existing file
            setting_table_loading_test_function() # Load the setting table normally, as in when using the file browser via button click


########################################################################################################################


max_cpu = multiprocessing.cpu_count() # Read available CPU cores
current_dir_string_var = tk.StringVar(value = default_current_dir) # Path the the last browsed folder

read_mode_string_var = tk.StringVar() # Argument for the --mode flag
peak_type_string_var = tk.StringVar() # Argument for the --peak flag
sample_table_question_string_var = tk.StringVar()
sample_table_string_var = tk.StringVar() # Path to setting table file which values are to be loaded
chip_rep_number_int_var = tk.IntVar(value = 1)
ctrl_rep_number_int_var = tk.IntVar(value = 1)
sample_table_notification_string_var = tk.StringVar() # GUI text notification for the sample table loading status
chip_sample_notification_string_var = tk.StringVar() # GUI text notification for the ChIP samples check status
ctrl_sample_notification_string_var = tk.StringVar() # GUI text notification for the control samples check status
setting_table_question_string_var = tk.StringVar()
setting_table_string_var = tk.StringVar(value = default_setting_table_file_full_path) # Path to setting table file which values are to be loaded
setting_table_notification_string_var = tk.StringVar(value = "Currently using default settings table") # GUI text notification for the setting table loading status
genome_ref_string_var = tk.StringVar(value = genome_ref_options[0]) # Argument for the --ref flag
genome_folder_string_var = tk.StringVar(value = genome_folder_full_path) # Argument for --genome flag
known_motif_question_string_var = tk.StringVar()
known_motif_string_var = tk.StringVar() # Argument for the --motif flag
setname_string_var = tk.StringVar() # Argument for the --setname flag
output_folder_string_var = tk.StringVar() # Argument for the --output flag
fcmerge_var = tk.IntVar() # Binary value for the --fcmerge flag
goann_var = tk.IntVar() # Binary value for the --goann flag
pathann_var = tk.IntVar() # Binary value for the --pathann flag
deltemp_var = tk.IntVar(value = 1) # Binary value for the --deltemp flag.
stdout_var = tk.IntVar() # Binary value for the 1> and 2> channel flags
homer_motif_question_string_var = tk.StringVar()
homer_motif_string_var = tk.StringVar(value = homer_motif_options[0]) # Argument for --homer_motif flag
meme_motif_question_string_var = tk.StringVar()
meme_motif_string_var = tk.StringVar(value = meme_motif_options[0]) # Argument for --meme_motif flag
cpu_count_string_var = tk.StringVar() # Argument for --thread flag
cpu_count_notification_string_var = tk.StringVar() # GUI text notification for inputted number of cores validity, and overall readiness of the command line

# Argument values for ChIP and control samples, obtained from loaded sample table or assigned manually
chip_rep1_r1_string_var = tk.StringVar() # 1st argument for the --chipR1 flag
chip_rep1_r2_string_var = tk.StringVar() # 1st argument for the --chipR2 flag
chip_rep2_r1_string_var = tk.StringVar() # 2nd argument for the --chipR1 flag
chip_rep2_r2_string_var = tk.StringVar() # 2nd argument for the --chipR2 flag
chip_rep3_r1_string_var = tk.StringVar() # 3rd argument for the --chipR1 flag
chip_rep3_r2_string_var = tk.StringVar() # 3rd argument for the --chipR2 flag
chip_rep4_r1_string_var = tk.StringVar() # 4th argument for the --chipR1 flag
chip_rep4_r2_string_var = tk.StringVar() # 4th argument for the --chipR2 flag
chip_rep5_r1_string_var = tk.StringVar() # 5th argument for the --chipR1 flag
chip_rep5_r2_string_var = tk.StringVar() # 5th argument for the --chipR2 flag
ctrl_rep1_r1_string_var = tk.StringVar() # 1st argument for the --ctrlR1 flag
ctrl_rep1_r2_string_var = tk.StringVar() # 1st argument for the --ctrlR2 flag
ctrl_rep2_r1_string_var = tk.StringVar() # 2nd argument for the --ctrlR1 flag
ctrl_rep2_r2_string_var = tk.StringVar() # 2nd argument for the --ctrlR2 flag
ctrl_rep3_r1_string_var = tk.StringVar() # 3rd argument for the --ctrlR1 flag
ctrl_rep3_r2_string_var = tk.StringVar() # 3rd argument for the --ctrlR2 flag
ctrl_rep4_r1_string_var = tk.StringVar() # 4th argument for the --ctrlR1 flag
ctrl_rep4_r2_string_var = tk.StringVar() # 4th argument for the --ctrlR2 flag
ctrl_rep5_r1_string_var = tk.StringVar() # 5th argument for the --ctrlR1 flag
ctrl_rep5_r2_string_var = tk.StringVar() # 5th argument for the --ctrlR2 flag

# Flags and argument values for every program in the pipeline, obtained from loaded setting table
fastqc1_arg = tk.StringVar()
clumpify_arg = tk.StringVar()
bbduk_arg = tk.StringVar()
trimmomatic_arg = tk.StringVar()
fastqc2_arg = tk.StringVar()
bwa_mem_arg = tk.StringVar()
samtools_view_arg = tk.StringVar()
plotfingerprint_arg = tk.StringVar()
fastqc3_arg = tk.StringVar()
macs2_callpeak_arg = tk.StringVar()
gem_arg = tk.StringVar()
sicer2_arg = tk.StringVar()
homer_findPeaks_arg = tk.StringVar()
genrich_arg = tk.StringVar()
homer_mergePeaks_arg = tk.StringVar()
homer_annotatePeaks_arg = tk.StringVar()
fold_change_calculator_arg = tk.StringVar()
homer_findMotifsGenome_arg = tk.StringVar()
meme_chip_arg = tk.StringVar()

# Flags and argument values for ChIP-AP command line call, obtained from this GUI
command_line_output_string_var = tk.StringVar() # The full ChIP-AP command line. A string of flags and arguments for all the variables below.
read_mode_arg = tk.StringVar()
peak_type_arg = tk.StringVar()
output_folder_arg = tk.StringVar()
setname_arg = tk.StringVar()
output_dir_arg = tk.StringVar()
genome_ref_arg = tk.StringVar()
genome_folder_arg = tk.StringVar()
sample_table_arg = tk.StringVar()
setting_table_arg = tk.StringVar()
known_motif_arg = tk.StringVar()
fcmerge_arg = tk.StringVar()
goann_arg = tk.StringVar()
pathann_arg = tk.StringVar()
deltemp_arg = tk.StringVar()
stdout_arg = tk.StringVar()
stderr_arg = tk.StringVar()
cpu_count_arg = tk.StringVar()
homer_motif_arg = tk.StringVar()
meme_motif_arg = tk.StringVar()


########################################################################################################################


# Traced variables. The assigned function will be triggered every time the variable's value changes.
# For variables which value only affects the resulting ChIP-AP command line and nothing else, update_command_line_function is used.
# For variables with more complex rules and greater influence over the values of other variables, other functions are used.
# See the comments of all the tracing functions (above) to learn what things are done when any of these traced variables change in value.
read_mode_string_var.trace('w', check_traced_input_function)
peak_type_string_var.trace('w', check_traced_input_function)
chip_rep1_r1_string_var.trace('w', register_sample_function)
chip_rep1_r2_string_var.trace('w', register_sample_function)
chip_rep2_r1_string_var.trace('w', register_sample_function)
chip_rep2_r2_string_var.trace('w', register_sample_function)
chip_rep3_r1_string_var.trace('w', register_sample_function)
chip_rep3_r2_string_var.trace('w', register_sample_function)
chip_rep4_r1_string_var.trace('w', register_sample_function)
chip_rep4_r2_string_var.trace('w', register_sample_function)
chip_rep5_r1_string_var.trace('w', register_sample_function)
chip_rep5_r2_string_var.trace('w', register_sample_function)
ctrl_rep1_r1_string_var.trace('w', register_sample_function)
ctrl_rep1_r2_string_var.trace('w', register_sample_function)
ctrl_rep2_r1_string_var.trace('w', register_sample_function)
ctrl_rep2_r2_string_var.trace('w', register_sample_function)
ctrl_rep3_r1_string_var.trace('w', register_sample_function)
ctrl_rep3_r2_string_var.trace('w', register_sample_function)
ctrl_rep4_r1_string_var.trace('w', register_sample_function)
ctrl_rep4_r2_string_var.trace('w', register_sample_function)
ctrl_rep5_r1_string_var.trace('w', register_sample_function)
ctrl_rep5_r2_string_var.trace('w', register_sample_function)
sample_table_question_string_var.trace('w', register_sample_function)
sample_table_string_var.trace('w', check_traced_input_function)
sample_table_string_var.trace('w', sample_table_entry_trace_load_function)
setting_table_question_string_var.trace('w', check_required_input_function)
setting_table_string_var.trace('w', check_traced_input_function)
setting_table_string_var.trace('w', setting_table_entry_trace_load_function)
genome_ref_string_var.trace('w', check_traced_input_function)
genome_folder_string_var.trace('w', check_traced_input_function)
known_motif_question_string_var.trace('w', check_required_input_function)
known_motif_string_var.trace('w', check_traced_input_function)
setname_string_var.trace('w', check_traced_input_function)
output_folder_string_var.trace('w', check_traced_input_function)
fcmerge_var.trace('w', update_command_line_function)
goann_var.trace('w', update_command_line_function)
pathann_var.trace('w', update_command_line_function)
deltemp_var.trace('w', update_command_line_function)
stdout_var.trace('w', update_command_line_function)
homer_motif_question_string_var.trace('w', check_traced_input_function)
homer_motif_string_var.trace('w', check_traced_input_function)
meme_motif_question_string_var.trace('w', check_traced_input_function)
meme_motif_string_var.trace('w', check_traced_input_function)
cpu_count_string_var.trace('w', check_traced_input_function)


########################################################################################################################


opening_frame = tk.Frame(root) # Make a new frame

# Insert full ChIP-AP logo
image = Image.open(chipap_logo_full_path)
image = image.resize((800, 295), Image.ANTIALIAS)
photo = ImageTk.PhotoImage(image)

# Place the logo in the center of the frame
chipap_logo = tk.Label(opening_frame, image = photo, anchor = "center")
chipap_logo.image = photo
chipap_logo.grid(row = 2, column = 1, padx = 10, pady = 10, columnspan = 2)

# Write few sentences about ChIP-AP below the logo
chipap_about = tk.Text(opening_frame, width = 70, height = 7, relief = tk.FLAT, bg = 'gray85',font = (None, 15))
chipap_about_text = 'ChIP-AP – Integrated Analysis Pipeline for Unbiased ChIP-seq Analysis.\n\nComplete guides and walkthroughs can be found on our github page\n(https://github.com/JSuryatenggara/ChIP-AP).\n\nIf you use ChIP-AP please cite our bioRxiv pre-print article\n(https://www.biorxiv.org/content/10.1101/2021.04.18.440382v1).'
chipap_about.insert(tk.END, chipap_about_text)
chipap_about.config(state = tk.DISABLED)
chipap_about.tag_configure("center", justify = tk.CENTER)
chipap_about.tag_add("center", 1.0, tk.END)
chipap_about.grid(row = 3, column = 1, columnspan = 2, sticky = "we", padx = 10, pady = (5,10))

# Create navigation buttons
opening_frame_exit_button = tk.Button(opening_frame, text = "Exit wizard", command = lambda : exit(), width = 20)
opening_frame_continue_button = tk.Button(opening_frame, text = "Continue >>", command = lambda : change_frame_function(opening_frame, read_mode_frame), width = 20)
opening_frame_exit_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
opening_frame_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))

opening_frame.grid(row = 0, column = 0, sticky = "news") # Display this frame immediately after finished loading, as the front page of the wizard


########################################################################################################################


read_mode_frame = tk.Frame(root) # Make a new frame

# Simple label with static text
read_mode_label = tk.Label(read_mode_frame, text = "Dataset sequencing mode:", justify = tk.LEFT, width = 30)
read_mode_label.grid(row = 1, column = 1, padx = 10, pady = 5, columnspan = 2)

# Radio button of the first choice that triggers a function when activated
single_end_radio = tk.Radiobutton(read_mode_frame, text = "Single end", padx = 5, variable = read_mode_string_var, value = 'single', width = 30)
CreateToolTip(single_end_radio, text = 'Select this if there is one\noutput file per sample')
single_end_radio.grid(row = 2, column = 1, padx = 10, pady = 2, columnspan = 2)

# Radio button of the second choice that triggers a function when activated
paired_end_radio = tk.Radiobutton(read_mode_frame, text = "Paired end", padx = 5, variable = read_mode_string_var, value = 'paired', width = 30)
CreateToolTip(paired_end_radio, text = 'Select this if there are two\noutput files per sample.\nThey are typically in pairs:\nR1 and R2 for every sample')
paired_end_radio.grid(row = 3, column = 1, padx = 10, pady = 2, columnspan = 2)

# Create navigation buttons
read_mode_back_button = tk.Button(read_mode_frame, text = "<< Back", command = lambda : change_frame_function(read_mode_frame, opening_frame), width = 20)
read_mode_continue_button = tk.Button(read_mode_frame, text = "Continue >>", command = lambda : change_frame_function(read_mode_frame, peak_type_frame), width = 20)
read_mode_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
read_mode_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################


peak_type_frame = tk.Frame(root) # Make a new frame

# Simple label with static text
peak_type_label = tk.Label(peak_type_frame, text = "Dataset peak type:", justify = tk.LEFT, width = 30)
peak_type_label.grid(row = 1, column = 1, padx = 10, pady = 5, columnspan = 2)

# Radio button of the first choice that triggers a function when activated
narrow_peak_radio = tk.Radiobutton(peak_type_frame, text = "Narrow peaks", padx = 5, variable = peak_type_string_var, value = 'narrow', width = 30)
CreateToolTip(narrow_peak_radio, text = 'Select this for ChIP-seq experiment\nusing transcription factor protein')
narrow_peak_radio.grid(row = 2, column = 1, padx = 10, pady = 2, columnspan = 2)

# Radio button of the second choice that triggers a function when activated
broad_peak_radio = tk.Radiobutton(peak_type_frame, text = "Broad peaks", padx = 5, variable = peak_type_string_var, value = 'broad', width = 30)
CreateToolTip(broad_peak_radio, text = 'Select this for ChIP-seq experiment\nusing chromatin modifier protein')
broad_peak_radio.grid(row = 3, column = 1, padx = 10, pady = 2, columnspan = 2)

# Radio button of the third choice that triggers a function when activated
unsure_peak_radio = tk.Radiobutton(peak_type_frame, text = "Unsure", padx = 5, variable = peak_type_string_var, value = 'unsure', width = 30)
CreateToolTip(unsure_peak_radio, text = 'Select this for ChIP-seq experiment\nusing protein with both possibilities')
unsure_peak_radio.grid(row = 4, column = 1, padx = 10, pady = 2, columnspan = 2)

# Create navigation buttons
peak_type_back_button = tk.Button(peak_type_frame, text = "<< Back", command = lambda : change_frame_function(peak_type_frame, read_mode_frame), width = 20)
peak_type_continue_button = tk.Button(peak_type_frame, text = "Continue >>", command = lambda : change_frame_function(peak_type_frame, sample_table_frame), width = 20)
peak_type_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
peak_type_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################


# Function triggered by the "Load sample table" button, to open a browser window for choosing a local file
def sample_table_button_function():
    sample_table_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a sample table file')
    if sample_table_button_input: # If the user choose a file from the browser
        sample_table_string_var.set(sample_table_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(sample_table_button_input.name)) # Set the last visited directory to the user-chosen file directory

    else: # If the user closed the browser without choosing any file from the browser
        clear_sample_function('withfile') # Clear all sample-related variables
        sample_table_notification_string_var.set('No sample table was selected. No samples were loaded') # Notify the user
        sample_table_notification_label.config(fg = 'red')
        rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path
        return # Exit the function


# Function to test whether the sample table is valid and its contents can be loaded properly without any error
def sample_table_loading_test_function():
    try: # Try to load the sample table with pd.read_csv
        sample_table_absolute_path = sample_table_string_var.get()
        sample_table_df = pd.read_csv(sample_table_absolute_path, delimiter='\t')
        sample_table_df.fillna('', inplace = True) # Replace with empty strings all 'NaN' fields resulting from empty fields in the sample table

    except: # If sample table loading fails
        # Clear all sample-related lists
        chip_r1_sample_list = []
        ctrl_r1_sample_list = []
        chip_r2_sample_list = []
        ctrl_r2_sample_list = []
        clear_sample_function('withfile') # Clear all sample-related variables
        sample_table_notification_string_var.set('Sample table loading error! No samples were loaded') # Notify the user
        sample_table_notification_label.config(fg = 'red')
        return False # Exit the function and return False as a sign of failed test
    
    # If no loading error occured up to this point
    try: # Try to register the contents into the prepared lists
        chip_r1_sample_list = [chip_r1_sample if str(chip_r1_sample) != 'nan' else '' for chip_r1_sample in sample_table_df['chip_read_1']]
        ctrl_r1_sample_list = [ctrl_r1_sample if str(ctrl_r1_sample) != 'nan' else '' for ctrl_r1_sample in sample_table_df['ctrl_read_1']]
        chip_r2_sample_list = [chip_r2_sample if str(chip_r2_sample) != 'nan' else '' for chip_r2_sample in sample_table_df['chip_read_2']]
        ctrl_r2_sample_list = [ctrl_r2_sample if str(ctrl_r2_sample) != 'nan' else '' for ctrl_r2_sample in sample_table_df['ctrl_read_2']]

    except: # If the read sample table contents failed to be loaded into the prepared lists
        # Clear all sample-related lists
        chip_r1_sample_list = []
        ctrl_r1_sample_list = []
        chip_r2_sample_list = []
        ctrl_r2_sample_list = []
        clear_sample_function('withfile') # Clear all sample-related variables
        sample_table_notification_string_var.set('Sample table reading error! No samples were loaded') # Notify the user
        sample_table_notification_label.config(fg = 'red')
        return False # Exit the function and return False as a sign of failed test

    if read_mode_string_var.get() == 'single': # In case of single-end dataset

        if not all(chip_r2 == '' for chip_r2 in chip_r2_sample_list) or not all(ctrl_r2 == '' for ctrl_r2 in ctrl_r2_sample_list): # If there is any sample registered into the R2 lists
            clear_sample_function('withfile') # Clear all sample-related variables
            sample_table_notification_string_var.set('Error! Paired samples detected. No samples were loaded') # Notify the user
            sample_table_notification_label.config(fg = 'red')
            return False # Exit the function and return False as a sign of failed test
        
        else: # If there is no sample registered into the R2 lists
            try: # Set the GUI sample-related variables accordingly with respect to the order and number of samples in the ChIP R1 list contents
                chip_rep1_r1_string_var.set(chip_r1_sample_list[0])
                chip_rep2_r1_string_var.set(chip_r1_sample_list[1])
                chip_rep3_r1_string_var.set(chip_r1_sample_list[2])
                chip_rep4_r1_string_var.set(chip_r1_sample_list[3])
                chip_rep5_r1_string_var.set(chip_r1_sample_list[4])

            except: # When the list is out of contents to register to the GUI variable
                pass # Do nothing and return no error. Only as many variables as the number of samples in the list needs to be set

            try: # Set the GUI sample-related variables accordingly with respect to the order and number of samples in the control R1 list contents
                ctrl_rep1_r1_string_var.set(ctrl_r1_sample_list[0])
                ctrl_rep2_r1_string_var.set(ctrl_r1_sample_list[1])
                ctrl_rep3_r1_string_var.set(ctrl_r1_sample_list[2])
                ctrl_rep4_r1_string_var.set(ctrl_r1_sample_list[3])
                ctrl_rep5_r1_string_var.set(ctrl_r1_sample_list[4])

            except: # When the list is out of contents to register to the GUI variable
                pass # Do nothing and return no error. Only as many variables as the number of samples in the list needs to be set
            
            # If all tests pass up to this point with no error
            sample_table_notification_string_var.set('Sample table loading successful') # Notify the user
            sample_table_notification_label.config(fg = 'green')
            rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path
            return True # Exit the function and return False as a sign of failed test

    if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
        try: # Set the GUI sample-related variables accordingly with respect to the order and number of samples in the ChIP R1 and R2 list contents
            chip_rep1_r1_string_var.set(chip_r1_sample_list[0])
            chip_rep1_r2_string_var.set(chip_r2_sample_list[0])
            chip_rep2_r1_string_var.set(chip_r1_sample_list[1])
            chip_rep2_r2_string_var.set(chip_r2_sample_list[1])
            chip_rep3_r1_string_var.set(chip_r1_sample_list[2])
            chip_rep3_r2_string_var.set(chip_r2_sample_list[2])
            chip_rep4_r1_string_var.set(chip_r1_sample_list[3])
            chip_rep4_r2_string_var.set(chip_r2_sample_list[3])
            chip_rep5_r1_string_var.set(chip_r1_sample_list[4])
            chip_rep5_r2_string_var.set(chip_r2_sample_list[4])

        except: # When the list is out of contents to register to the GUI variable
            pass # Do nothing and return no error. Only as many variables as the number of samples in the list needs to be set

        try: # Set the GUI sample-related variables accordingly with respect to the order and number of samples in the control R1 and R2 list contents
            ctrl_rep1_r1_string_var.set(ctrl_r1_sample_list[0])
            ctrl_rep1_r2_string_var.set(ctrl_r2_sample_list[0])
            ctrl_rep2_r1_string_var.set(ctrl_r1_sample_list[1])
            ctrl_rep2_r2_string_var.set(ctrl_r2_sample_list[1])
            ctrl_rep3_r1_string_var.set(ctrl_r1_sample_list[2])
            ctrl_rep3_r2_string_var.set(ctrl_r2_sample_list[2])
            ctrl_rep4_r1_string_var.set(ctrl_r1_sample_list[3])
            ctrl_rep4_r2_string_var.set(ctrl_r2_sample_list[3])
            ctrl_rep5_r1_string_var.set(ctrl_r1_sample_list[4])
            ctrl_rep5_r2_string_var.set(ctrl_r2_sample_list[4])

        except: # When the list is out of contents to register to the GUI variable
            pass # Do nothing and return no error. Only as many variables as the number of samples in the list needs to be set

        # If all tests pass up to this point with no error
        sample_table_notification_string_var.set('Sample table loading successful') # Notify the user
        sample_table_notification_label.config(fg = 'green')
        rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path
        return True # Exit the function and return True as a sign of passed test


# Function that manages the widgets in the GUI's sample table loading frame
# Toggled by the user choice between loading or not loading a sample table
def sample_table_entry_popup_function(sample_table_entry_popup_arg):
    if sample_table_entry_popup_arg == 'yes': # If the user choose to load a sample table
        # Display the buttons and entry field for loading the sample table
        sample_table_entry.grid(row = 3, column = 1, padx = (10,5), ipady = 3)
        sample_table_button.grid(row = 3, column = 2, padx = (5,10), pady = 2)
        # Remove the drop-down options widget for number of replicates of ChIP and control samples
        chip_rep_number_label.grid_remove()
        chip_rep_number_drop_down.grid_remove()
        ctrl_rep_number_label.grid_remove()
        ctrl_rep_number_drop_down.grid_remove()
        clear_sample_function('withfile') # Clear all sample-related variables
        sample_table_notification_string_var.set('Please load your sample table file') # Prompt the user
        sample_table_notification_label.config(fg = 'blue')
        # Functions to keep the frame window (and its contents) at the center of screen
        root.eval('tk::PlaceWindow . center')
        root.grab_set()
        root.focus_force()

    if sample_table_entry_popup_arg == 'no': # If the user choose to load the samples manually
        # Remove the buttons and entry field for loading the sample table
        sample_table_entry.grid_remove()
        sample_table_button.grid_remove()
        # Display the drop-down options widget for number of replicates of ChIP and control samples
        chip_rep_number_label.grid(row = 5, column = 1, padx = (10,5), sticky = "e")
        chip_rep_number_drop_down.grid(row = 5, column = 2, padx = (5,10), sticky = "w")
        ctrl_rep_number_label.grid(row = 6, column = 1, padx = (10,5), sticky = "e")
        ctrl_rep_number_drop_down.grid(row = 6, column = 2, padx = (5,10), sticky = "w")
        clear_sample_function('withfile') # Clear all sample-related variables
        sample_table_notification_string_var.set('Please choose the number of sample replicates') # Prompt the user
        sample_table_notification_label.config(fg = 'blue')
        # Functions to keep the frame window (and its contents) at the center of screen
        root.eval('tk::PlaceWindow . center')
        root.grab_set()
        root.focus_force()

    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path


# Function to clear all sample-related variables in the GUI
def clear_sample_function(clear_sample_function_arg):
    
    if clear_sample_function_arg == 'withfile': # When the user wants to also clear the registered path to the sample table file
        sample_table_string_var.set('') # Register empty strings to the variable that stores the path to the sample table file

    # Register empty strings to all ChIP samples GUI variables
    chip_rep1_r1_string_var.set('')
    chip_rep1_r2_string_var.set('')
    chip_rep2_r1_string_var.set('')
    chip_rep2_r2_string_var.set('')
    chip_rep3_r1_string_var.set('')
    chip_rep3_r2_string_var.set('')
    chip_rep4_r1_string_var.set('')
    chip_rep4_r2_string_var.set('')
    chip_rep5_r1_string_var.set('')
    chip_rep5_r2_string_var.set('')

    # Register empty strings to all control samples GUI variables
    ctrl_rep1_r1_string_var.set('')
    ctrl_rep1_r2_string_var.set('')
    ctrl_rep2_r1_string_var.set('')
    ctrl_rep2_r2_string_var.set('')
    ctrl_rep3_r1_string_var.set('')
    ctrl_rep3_r2_string_var.set('')
    ctrl_rep4_r1_string_var.set('')
    ctrl_rep4_r2_string_var.set('')
    ctrl_rep5_r1_string_var.set('')
    ctrl_rep5_r2_string_var.set('')

    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path


sample_table_frame = tk.Frame(root) # Make a new frame

# Simple label with static text
sample_table_label = tk.Label(sample_table_frame, text = "Do you want to use your sample table?", width = 70)
sample_table_label.grid(row = 1, column = 1, padx = 10, pady = 5, columnspan = 2)

# Radio button of the first choice that triggers a function when activated
sample_table_yes_radio = tk.Radiobutton(sample_table_frame, text = "Yes, I would like use my sample table", padx = 5, variable = sample_table_question_string_var, value = 'yes', command = lambda : sample_table_entry_popup_function('yes'), width = 50)
CreateToolTip(sample_table_yes_radio, text = 'Select this to load a pre-populated\nsample table in ChIP-AP format\n(see GitHub documentation).\nThe ChIP and control samples list\nwill be filled automatically\nbased on the loaded sample table')
sample_table_yes_radio.grid(row = 2, column = 1, padx = 10, pady = 2, columnspan = 2)

# Entry field where user can manually key in their input to be stored as a variable
sample_table_entry = tk.Entry(sample_table_frame, textvariable = sample_table_string_var, width = 50, justify = tk.RIGHT)
sample_table_entry.xview_moveto(1)

# Clickable button that triggers a function when activated
sample_table_button = tk.Button(sample_table_frame, text = 'Browse', width = 20, command = lambda : sample_table_button_function())
CreateToolTip(sample_table_button, text = 'Click here to browse and\nselect your sample table file')

# Radio button of the first choice that triggers a function when activated
sample_table_no_radio = tk.Radiobutton(sample_table_frame, text = "No, I would like to assign my samples manually", padx = 5, variable = sample_table_question_string_var, value = 'no', command = lambda : sample_table_entry_popup_function('no'), width = 50)
CreateToolTip(sample_table_no_radio, text = 'Select this if you want to\ntype in, or browse and select\nall your sample files yourself')
sample_table_no_radio.grid(row = 4, column = 1, padx = 10, pady = 2, columnspan = 2)

replicate_number_choice = [1, 2, 3, 4, 5] # Choices for the number of replicates drop-down menu

# Simple label with static text
chip_rep_number_label = tk.Label(sample_table_frame, text = "Number of ChIP replicate(s):")

# Drop-down options widget for number of replicates of ChIP samples
chip_rep_number_drop_down = tk.OptionMenu(sample_table_frame, chip_rep_number_int_var, *replicate_number_choice)
chip_rep_number_drop_down.config(width = 3, takefocus = 1)
chip_rep_number_drop_down_menu = root.nametowidget(chip_rep_number_drop_down.menuname) # Get the drop-down menu options to follow the default fonts
chip_rep_number_drop_down_menu.config(font = default_font) # Get the drop-down menu options to follow the default fonts

# Simple label with static text
ctrl_rep_number_label = tk.Label(sample_table_frame, text = "Number of control replicate(s):")

# Drop-down options widget for number of replicates of control samples
ctrl_rep_number_drop_down = tk.OptionMenu(sample_table_frame, ctrl_rep_number_int_var, *replicate_number_choice)
ctrl_rep_number_drop_down.config(width = 3, takefocus = 1)
ctrl_rep_number_drop_down_menu = root.nametowidget(ctrl_rep_number_drop_down.menuname) # Get the drop-down menu options to follow the default fonts
ctrl_rep_number_drop_down_menu.config(font = default_font) # Get the drop-down menu options to follow the default fonts

# Dynamic label with changeable text
sample_table_notification_label = tk.Label(sample_table_frame, textvariable = sample_table_notification_string_var, width = 70, padx = 5, pady = 5)
sample_table_notification_label.grid(row = 7, column = 1, padx = 10, pady = 5, columnspan = 2)

# Create navigation buttons
sample_table_back_button = tk.Button(sample_table_frame, text = "<< Back", command = lambda : change_frame_function(sample_table_frame, peak_type_frame), width = 20)
sample_table_continue_button = tk.Button(sample_table_frame, text = "Continue >>", command = lambda : change_frame_function(sample_table_frame, chip_list_frame), width = 20)
sample_table_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
sample_table_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################


# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def chip_rep1_r1_button_function():
    chip_rep1_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = chip_list_frame, mode = 'rb', title = 'Choose a file')
    if chip_rep1_r1_button_input: # If the user choose a file from the browser
        chip_rep1_r1_string_var.set(chip_rep1_r1_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(chip_rep1_r1_button_input.name)) # Set the last visited directory to the user-chosen file directory
    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def chip_rep1_r2_button_function():
    chip_rep1_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = chip_list_frame, mode = 'rb', title = 'Choose a file')
    if chip_rep1_r2_button_input: # If the user choose a file from the browser
        chip_rep1_r2_string_var.set(chip_rep1_r2_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(chip_rep1_r2_button_input.name)) # Set the last visited directory to the user-chosen file directory
    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def chip_rep2_r1_button_function():
    chip_rep2_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = chip_list_frame, mode = 'rb', title = 'Choose a file')
    if chip_rep2_r1_button_input: # If the user choose a file from the browser
        chip_rep2_r1_string_var.set(chip_rep2_r1_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(chip_rep2_r1_button_input.name)) # Set the last visited directory to the user-chosen file directory
    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def chip_rep2_r2_button_function():
    chip_rep2_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = chip_list_frame, mode = 'rb', title = 'Choose a file')
    if chip_rep2_r2_button_input: # If the user choose a file from the browser
        chip_rep2_r2_string_var.set(chip_rep2_r2_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(chip_rep2_r2_button_input.name)) # Set the last visited directory to the user-chosen file directory
    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def chip_rep3_r1_button_function():
    chip_rep3_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = chip_list_frame, mode = 'rb', title = 'Choose a file')
    if chip_rep3_r1_button_input: # If the user choose a file from the browser
        chip_rep3_r1_string_var.set(chip_rep3_r1_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(chip_rep3_r1_button_input.name)) # Set the last visited directory to the user-chosen file directory
    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def chip_rep3_r2_button_function():
    chip_rep3_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = chip_list_frame, mode = 'rb', title = 'Choose a file')
    if chip_rep3_r2_button_input: # If the user choose a file from the browser
        chip_rep3_r2_string_var.set(chip_rep3_r2_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(chip_rep3_r2_button_input.name)) # Set the last visited directory to the user-chosen file directory
    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def chip_rep4_r1_button_function():
    chip_rep4_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = chip_list_frame, mode = 'rb', title = 'Choose a file')
    if chip_rep4_r1_button_input: # If the user choose a file from the browser
        chip_rep4_r1_string_var.set(chip_rep4_r1_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(chip_rep4_r1_button_input.name)) # Set the last visited directory to the user-chosen file directory
    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def chip_rep4_r2_button_function():
    chip_rep4_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = chip_list_frame, mode = 'rb', title = 'Choose a file')
    if chip_rep4_r2_button_input: # If the user choose a file from the browser
        chip_rep4_r2_string_var.set(chip_rep4_r2_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(chip_rep4_r2_button_input.name)) # Set the last visited directory to the user-chosen file directory
    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def chip_rep5_r1_button_function():
    chip_rep5_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = chip_list_frame, mode = 'rb', title = 'Choose a file')
    if chip_rep5_r1_button_input: # If the user choose a file from the browser
        chip_rep5_r1_string_var.set(chip_rep5_r1_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(chip_rep5_r1_button_input.name)) # Set the last visited directory to the user-chosen file directory
    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def chip_rep5_r2_button_function():
    chip_rep5_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = chip_list_frame, mode = 'rb', title = 'Choose a file')
    if chip_rep5_r2_button_input: # If the user choose a file from the browser
        chip_rep5_r2_string_var.set(chip_rep5_r2_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(chip_rep5_r2_button_input.name)) # Set the last visited directory to the user-chosen file directory
    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path


# Function that manages the widgets in the GUI's ChIP samples loading frame
# Toggled by the user choice between loading and not loading a sample table, and the number of ChIP sample replicates
def display_chip_widget_function():
    for chip_widget in chip_widget_list: # Reset the frame contents by first removing every ChIP samples widgets (button and entry field)
        chip_widget.grid_remove()

    if sample_table_question_string_var.get() == 'yes': # In case of user choosing to load a sample table
        if bool(chip_rep1_r1_string_var.get()): # If this sample replicate exists
            # Display the sample loader button and entry field for the first read of this replicate 
            chip_rep1_r1_button.grid(column = 1, row = 1, padx = (10,5))
            chip_rep1_r1_entry.grid(column = 2, row = 1, padx = (5,10), ipady = 3)
            
            if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
                # Also display the sample loader button and entry field for the second read of this replicate
                chip_rep1_r2_button.grid(column = 1, row = 2, padx = (10,5))
                chip_rep1_r2_entry.grid(column = 2, row = 2, padx = (5,10), ipady = 3)

        if bool(chip_rep2_r1_string_var.get()): # If this sample replicate exists
            # Display the sample loader button and entry field for the first read of this replicate 
            chip_rep2_r1_button.grid(column = 1, row = 3, padx = (10,5), pady = (10,0))
            chip_rep2_r1_entry.grid(column = 2, row = 3, padx = (5,10), ipady = 3, pady = (10,0))
            
            if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
                # Also display the sample loader button and entry field for the second read of this replicate
                chip_rep2_r2_button.grid(column = 1, row = 4, padx = (10,5))
                chip_rep2_r2_entry.grid(column = 2, row = 4, padx = (5,10), ipady = 3)
    
        if bool(chip_rep3_r1_string_var.get()): # If this sample replicate exists
            # Display the sample loader button and entry field for the first read of this replicate
            chip_rep3_r1_button.grid(column = 1, row = 5, padx = (10,5), pady = (10,0))
            chip_rep3_r1_entry.grid(column = 2, row = 5, padx = (5,10), ipady = 3, pady = (10,0))     
            
            if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
                # Also display the sample loader button and entry field for the second read of this replicate
                chip_rep3_r2_button.grid(column = 1, row = 6, padx = (10,5))
                chip_rep3_r2_entry.grid(column = 2, row = 6, padx = (5,10), ipady = 3)

        if bool(chip_rep4_r1_string_var.get()): # If this sample replicate exists
            # Display the sample loader button and entry field for the first read of this replicate
            chip_rep4_r1_button.grid(column = 1, row = 7, padx = (10,5), pady = (10,0))
            chip_rep4_r1_entry.grid(column = 2, row = 7, padx = (5,10), ipady = 3, pady = (10,0))
            
            if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
                # Also display the sample loader button and entry field for the second read of this replicate
                chip_rep4_r2_button.grid(column = 1, row = 8, padx = (10,5))
                chip_rep4_r2_entry.grid(column = 2, row = 8, padx = (5,10), ipady = 3)

        if bool(chip_rep5_r1_string_var.get()): # If this sample replicate exists
            # Display the sample loader button and entry field for the first read of this replicate
            chip_rep5_r1_button.grid(column = 1, row = 9, padx = (10,5), pady = (10,0))
            chip_rep5_r1_entry.grid(column = 2, row = 9, padx = (5,10), ipady = 3, pady = (10,0))
            
            if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
                # Also display the sample loader button and entry field for the second read of this replicate
                chip_rep5_r2_button.grid(column = 1, row = 10, padx = (10,5))
                chip_rep5_r2_entry.grid(column = 2, row = 10, padx = (5,10), ipady = 3)

    if sample_table_question_string_var.get() == 'no': # In case of user choosing to load the samples manually
        if chip_rep_number_int_var.get() >= 1: # If this sample replicate exists
            # Display the sample loader button and entry field for the first read of this replicate
            chip_rep1_r1_button.grid(column = 1, row = 1, padx = (10,5))
            chip_rep1_r1_entry.grid(column = 2, row = 1, padx = (5,10), ipady = 3)
            
            if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
                # Also display the sample loader button and entry field for the second read of this replicate
                chip_rep1_r2_button.grid(column = 1, row = 2, padx = (10,5))
                chip_rep1_r2_entry.grid(column = 2, row = 2, padx = (5,10), ipady = 3)

        if chip_rep_number_int_var.get() >= 2: # If this sample replicate exists
            # Display the sample loader button and entry field for the first read of this replicate
            chip_rep2_r1_button.grid(column = 1, row = 3, padx = (10,5), pady = (10,0))
            chip_rep2_r1_entry.grid(column = 2, row = 3, padx = (5,10), ipady = 3, pady = (10,0))
            
            if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
                # Also display the sample loader button and entry field for the second read of this replicate
                chip_rep2_r2_button.grid(column = 1, row = 4, padx = (10,5))
                chip_rep2_r2_entry.grid(column = 2, row = 4, padx = (5,10), ipady = 3)
    
        if chip_rep_number_int_var.get() >= 3: # If this sample replicate exists
            # Display the sample loader button and entry field for the first read of this replicate
            chip_rep3_r1_button.grid(column = 1, row = 5, padx = (10,5), pady = (10,0))
            chip_rep3_r1_entry.grid(column = 2, row = 5, padx = (5,10), ipady = 3, pady = (10,0))
            
            if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
                # Also display the sample loader button and entry field for the second read of this replicate
                chip_rep3_r2_button.grid(column = 1, row = 6, padx = (10,5))
                chip_rep3_r2_entry.grid(column = 2, row = 6, padx = (5,10), ipady = 3)

        if chip_rep_number_int_var.get() >= 4: # If this sample replicate exists
            # Display the sample loader button and entry field for the first read of this replicate
            chip_rep4_r1_button.grid(column = 1, row = 7, padx = (10,5), pady = (10,0))
            chip_rep4_r1_entry.grid(column = 2, row = 7, padx = (5,10), ipady = 3, pady = (10,0))

            if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
                # Also display the sample loader button and entry field for the second read of this replicate
                chip_rep4_r2_button.grid(column = 1, row = 8, padx = (10,5))
                chip_rep4_r2_entry.grid(column = 2, row = 8, padx = (5,10), ipady = 3)

        if chip_rep_number_int_var.get() >= 5: # If this sample replicate exists
            # Display the sample loader button and entry field for the first read of this replicate
            chip_rep5_r1_button.grid(column = 1, row = 9, padx = (10,5), pady = (10,0))
            chip_rep5_r1_entry.grid(column = 2, row = 9, padx = (5,10), ipady = 3, pady = (10,0))

            if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
                # Also display the sample loader button and entry field for the second read of this replicate
                chip_rep5_r2_button.grid(column = 1, row = 10, padx = (10,5))
                chip_rep5_r2_entry.grid(column = 2, row = 10, padx = (5,10), ipady = 3)

    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path


chip_list_frame = tk.Frame(root) # Make a new frame

# Simple label with static text
chip_list_label = tk.Label(chip_list_frame, text = "Please type in, or browse and select your ChIP sample files below", justify = tk.LEFT, width = 70)
chip_list_label.grid(row = 0, column = 1, padx = 10, pady = (5,10), columnspan = 2)

# Clickable button that triggers a function when activated and an entry field where user can manually key in their input to be stored as a variable
# For each of the five (maximum number that can be handled by ChIP-AP) first reads of the ChIP sample replicates, and each of the five second reads 
chip_rep1_r1_button = tk.Button(chip_list_frame, text ='ChIP rep 1 read 1', command = lambda : chip_rep1_r1_button_function(), width = 20)
chip_rep1_r1_entry = tk.Entry(chip_list_frame, textvariable = chip_rep1_r1_string_var, width = 50, justify = tk.RIGHT)
chip_rep1_r2_button = tk.Button(chip_list_frame, text ='ChIP rep 1 read 2', command = lambda : chip_rep1_r2_button_function(), width = 20)
chip_rep1_r2_entry = tk.Entry(chip_list_frame, textvariable = chip_rep1_r2_string_var, width = 50, justify = tk.RIGHT)
chip_rep2_r1_button = tk.Button(chip_list_frame, text ='ChIP rep 2 read 1', command = lambda : chip_rep2_r1_button_function(), width = 20)
chip_rep2_r1_entry = tk.Entry(chip_list_frame, textvariable = chip_rep2_r1_string_var, width = 50, justify = tk.RIGHT)
chip_rep2_r2_button = tk.Button(chip_list_frame, text ='ChIP rep 2 read 2', command = lambda : chip_rep2_r2_button_function(), width = 20)
chip_rep2_r2_entry = tk.Entry(chip_list_frame, textvariable = chip_rep2_r2_string_var, width = 50, justify = tk.RIGHT)
chip_rep3_r1_button = tk.Button(chip_list_frame, text ='ChIP rep 3 read 1', command = lambda : chip_rep3_r1_button_function(), width = 20)
chip_rep3_r1_entry = tk.Entry(chip_list_frame, textvariable = chip_rep3_r1_string_var, width = 50, justify = tk.RIGHT)
chip_rep3_r2_button = tk.Button(chip_list_frame, text ='ChIP rep 3 read 2', command = lambda : chip_rep3_r2_button_function(), width = 20)
chip_rep3_r2_entry = tk.Entry(chip_list_frame, textvariable = chip_rep3_r2_string_var, width = 50, justify = tk.RIGHT)
chip_rep4_r1_button = tk.Button(chip_list_frame, text ='ChIP rep 4 read 1', command = lambda : chip_rep4_r1_button_function(), width = 20)
chip_rep4_r1_entry = tk.Entry(chip_list_frame, textvariable = chip_rep4_r1_string_var, width = 50, justify = tk.RIGHT)
chip_rep4_r2_button = tk.Button(chip_list_frame, text ='ChIP rep 4 read 2', command = lambda : chip_rep4_r2_button_function(), width = 20)
chip_rep4_r2_entry = tk.Entry(chip_list_frame, textvariable = chip_rep4_r2_string_var, width = 50, justify = tk.RIGHT)
chip_rep5_r1_button = tk.Button(chip_list_frame, text ='ChIP rep 5 read 1', command = lambda : chip_rep5_r1_button_function(), width = 20)
chip_rep5_r1_entry = tk.Entry(chip_list_frame, textvariable = chip_rep5_r1_string_var, width = 50, justify = tk.RIGHT)
chip_rep5_r2_button = tk.Button(chip_list_frame, text ='ChIP rep 5 read 2', command = lambda : chip_rep5_r2_button_function(), width = 20)
chip_rep5_r2_entry = tk.Entry(chip_list_frame, textvariable = chip_rep5_r2_string_var, width = 50, justify = tk.RIGHT)

# The help tooltips for each of the ten buttons above
CreateToolTip(chip_rep1_r1_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam\nfile extension accepted)')
CreateToolTip(chip_rep1_r2_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam\nfile extension accepted)')
CreateToolTip(chip_rep2_r1_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam\nfile extension accepted)')
CreateToolTip(chip_rep2_r2_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam\nfile extension accepted)')
CreateToolTip(chip_rep3_r1_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam\nfile extension accepted)')
CreateToolTip(chip_rep3_r2_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam\nfile extension accepted)')
CreateToolTip(chip_rep4_r1_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam\nfile extension accepted)')
CreateToolTip(chip_rep4_r2_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam\nfile extension accepted)')
CreateToolTip(chip_rep5_r1_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam\nfile extension accepted)')
CreateToolTip(chip_rep5_r2_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam\nfile extension accepted)')

# Listing of all the ChIP samples-related widgets
chip_widget_list = [chip_rep1_r1_button, 
                    chip_rep1_r2_button, 
                    chip_rep2_r1_button, 
                    chip_rep2_r2_button, 
                    chip_rep3_r1_button, 
                    chip_rep3_r2_button, 
                    chip_rep4_r1_button, 
                    chip_rep4_r2_button, 
                    chip_rep5_r1_button, 
                    chip_rep5_r2_button, 
                    chip_rep1_r1_entry, 
                    chip_rep1_r2_entry, 
                    chip_rep2_r1_entry, 
                    chip_rep2_r2_entry, 
                    chip_rep3_r1_entry, 
                    chip_rep3_r2_entry, 
                    chip_rep4_r1_entry, 
                    chip_rep4_r2_entry, 
                    chip_rep5_r1_entry, 
                    chip_rep5_r2_entry]

# Dynamic label with changeable text
chip_sample_notification_label = tk.Label(chip_list_frame, textvariable = chip_sample_notification_string_var, width = 70, padx = 5, pady = 5)
chip_sample_notification_label.grid(row = 15, column = 1, padx = 10, pady = 5, columnspan = 2)

# Create navigation buttons
chip_list_frame_back_button = tk.Button(chip_list_frame, text = "<< Back", command = lambda : change_frame_function(chip_list_frame, sample_table_frame), width = 20)
chip_list_frame_continue_button = tk.Button(chip_list_frame, text = "Continue >>", command = lambda : change_frame_function(chip_list_frame, ctrl_list_frame), width = 20)
chip_list_frame_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
chip_list_frame_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################


# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def ctrl_rep1_r1_button_function():
    ctrl_rep1_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = ctrl_list_frame, mode = 'rb', title = 'Choose a file')
    if ctrl_rep1_r1_button_input: # If the user choose a file from the browser
        ctrl_rep1_r1_string_var.set(ctrl_rep1_r1_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(ctrl_rep1_r1_button_input.name)) # Set the last visited directory to the user-chosen file directory
    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def ctrl_rep1_r2_button_function():
    ctrl_rep1_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = ctrl_list_frame, mode = 'rb', title = 'Choose a file')
    if ctrl_rep1_r2_button_input: # If the user choose a file from the browser
        ctrl_rep1_r2_string_var.set(ctrl_rep1_r2_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(ctrl_rep1_r2_button_input.name)) # Set the last visited directory to the user-chosen file directory
    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def ctrl_rep2_r1_button_function():
    ctrl_rep2_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = ctrl_list_frame, mode = 'rb', title = 'Choose a file')
    if ctrl_rep2_r1_button_input: # If the user choose a file from the browser
        ctrl_rep2_r1_string_var.set(ctrl_rep2_r1_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(ctrl_rep2_r1_button_input.name)) # Set the last visited directory to the user-chosen file directory
    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def ctrl_rep2_r2_button_function():
    ctrl_rep2_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = ctrl_list_frame, mode = 'rb', title = 'Choose a file')
    if ctrl_rep2_r2_button_input: # If the user choose a file from the browser
        ctrl_rep2_r2_string_var.set(ctrl_rep2_r2_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(ctrl_rep2_r2_button_input.name)) # Set the last visited directory to the user-chosen file directory
    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def ctrl_rep3_r1_button_function():
    ctrl_rep3_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = ctrl_list_frame, mode = 'rb', title = 'Choose a file')
    if ctrl_rep3_r1_button_input: # If the user choose a file from the browser
        ctrl_rep3_r1_string_var.set(ctrl_rep3_r1_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(ctrl_rep3_r1_button_input.name)) # Set the last visited directory to the user-chosen file directory
    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def ctrl_rep3_r2_button_function():
    ctrl_rep3_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = ctrl_list_frame, mode = 'rb', title = 'Choose a file')
    if ctrl_rep3_r2_button_input: # If the user choose a file from the browser
        ctrl_rep3_r2_string_var.set(ctrl_rep3_r2_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(ctrl_rep3_r2_button_input.name)) # Set the last visited directory to the user-chosen file directory
    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def ctrl_rep4_r1_button_function():
    ctrl_rep4_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = ctrl_list_frame, mode = 'rb', title = 'Choose a file')
    if ctrl_rep4_r1_button_input: # If the user choose a file from the browser
        ctrl_rep4_r1_string_var.set(ctrl_rep4_r1_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(ctrl_rep4_r1_button_input.name)) # Set the last visited directory to the user-chosen file directory
    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def ctrl_rep4_r2_button_function():
    ctrl_rep4_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = ctrl_list_frame, mode = 'rb', title = 'Choose a file')
    if ctrl_rep4_r2_button_input: # If the user choose a file from the browser
        ctrl_rep4_r2_string_var.set(ctrl_rep4_r2_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(ctrl_rep4_r2_button_input.name)) # Set the last visited directory to the user-chosen file directory
    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def ctrl_rep5_r1_button_function():
    ctrl_rep5_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = ctrl_list_frame, mode = 'rb', title = 'Choose a file')
    if ctrl_rep5_r1_button_input: # If the user choose a file from the browser
        ctrl_rep5_r1_string_var.set(ctrl_rep5_r1_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(ctrl_rep5_r1_button_input.name)) # Set the last visited directory to the user-chosen file directory
    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def ctrl_rep5_r2_button_function():
    ctrl_rep5_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = ctrl_list_frame, mode = 'rb', title = 'Choose a file')
    if ctrl_rep5_r2_button_input: # If the user choose a file from the browser
        ctrl_rep5_r2_string_var.set(ctrl_rep5_r2_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(ctrl_rep5_r2_button_input.name)) # Set the last visited directory to the user-chosen file directory
    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path


# Function that manages the widgets in the GUI's control samples loading frame
# Toggled by the user choice between loading and not loading a sample table, and the number of control sample replicates
def display_ctrl_widget_function():
    for ctrl_widget in ctrl_widget_list: # Reset the frame contents by first removing every control samples widgets (button and entry field)
        ctrl_widget.grid_remove()

    if sample_table_question_string_var.get() == 'yes': # In case of user choosing to load a sample table
        if bool(ctrl_rep1_r1_string_var.get()): # If this sample replicate exists
            # Display the sample loader button and entry field for the first read of this replicate
            ctrl_rep1_r1_button.grid(column = 1, row = 1, padx = (10,5))
            ctrl_rep1_r1_entry.grid(column = 2, row = 1, padx = (5,10), ipady = 3)
        
            if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
                # Also display the sample loader button and entry field for the second read of this replicate
                ctrl_rep1_r2_button.grid(column = 1, row = 2, padx = (10,5))
                ctrl_rep1_r2_entry.grid(column = 2, row = 2, padx = (5,10), ipady = 3)

        if bool(ctrl_rep2_r1_string_var.get()): # If this sample replicate exists
            # Display the sample loader button and entry field for the first read of this replicate
            ctrl_rep2_r1_button.grid(column = 1, row = 3, padx = (10,5), pady = (10,0))
            ctrl_rep2_r1_entry.grid(column = 2, row = 3, padx = (5,10), ipady = 3, pady = (10,0))
            
            if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
                # Also display the sample loader button and entry field for the second read of this replicate
                ctrl_rep2_r2_button.grid(column = 1, row = 4, padx = (10,5))
                ctrl_rep2_r2_entry.grid(column = 2, row = 4, padx = (5,10), ipady = 3)
    
        if bool(ctrl_rep3_r1_string_var.get()): # If this sample replicate exists
            # Display the sample loader button and entry field for the first read of this replicate
            ctrl_rep3_r1_button.grid(column = 1, row = 5, padx = (10,5), pady = (10,0))
            ctrl_rep3_r1_entry.grid(column = 2, row = 5, padx = (5,10), ipady = 3, pady = (10,0))
            
            if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
                # Also display the sample loader button and entry field for the second read of this replicate
                ctrl_rep3_r2_button.grid(column = 1, row = 6, padx = (10,5))
                ctrl_rep3_r2_entry.grid(column = 2, row = 6, padx = (5,10), ipady = 3)

        if bool(ctrl_rep4_r1_string_var.get()): # If this sample replicate exists
            # Display the sample loader button and entry field for the first read of this replicate
            ctrl_rep4_r1_button.grid(column = 1, row = 7, padx = (10,5), pady = (10,0))
            ctrl_rep4_r1_entry.grid(column = 2, row = 7, padx = (5,10), ipady = 3, pady = (10,0))

            if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
                # Also display the sample loader button and entry field for the second read of this replicate
                ctrl_rep4_r2_button.grid(column = 1, row = 8, padx = (10,5))
                ctrl_rep4_r2_entry.grid(column = 2, row = 8, padx = (5,10), ipady = 3)

        if bool(ctrl_rep5_r1_string_var.get()): # If this sample replicate exists
            # Display the sample loader button and entry field for the first read of this replicate
            ctrl_rep5_r1_button.grid(column = 1, row = 9, padx = (10,5), pady = (10,0))
            ctrl_rep5_r1_entry.grid(column = 2, row = 9, padx = (5,10), ipady = 3, pady = (10,0))

            if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
                # Also display the sample loader button and entry field for the second read of this replicate
                ctrl_rep5_r2_button.grid(column = 1, row = 10, padx = (10,5))
                ctrl_rep5_r2_entry.grid(column = 2, row = 10, padx = (5,10), ipady = 3)

    if sample_table_question_string_var.get() == 'no': # In case of user choosing to load the samples manually
        if ctrl_rep_number_int_var.get() >= 1: # If this sample replicate exists
            # Display the sample loader button and entry field for the first read of this replicate
            ctrl_rep1_r1_button.grid(column = 1, row = 1, padx = (10,5))
            ctrl_rep1_r1_entry.grid(column = 2, row = 1, padx = (5,10), ipady = 3)
            
            if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
                # Also display the sample loader button and entry field for the second read of this replicate
                ctrl_rep1_r2_button.grid(column = 1, row = 2, padx = (10,5))
                ctrl_rep1_r2_entry.grid(column = 2, row = 2, padx = (5,10), ipady = 3)

        if ctrl_rep_number_int_var.get() >= 2: # If this sample replicate exists
            # Display the sample loader button and entry field for the first read of this replicate
            ctrl_rep2_r1_button.grid(column = 1, row = 3, padx = (10,5), pady = (10,0))
            ctrl_rep2_r1_entry.grid(column = 2, row = 3, padx = (5,10), ipady = 3, pady = (10,0))
            
            if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
                # Also display the sample loader button and entry field for the second read of this replicate
                ctrl_rep2_r2_button.grid(column = 1, row = 4, padx = (10,5))
                ctrl_rep2_r2_entry.grid(column = 2, row = 4, padx = (5,10), ipady = 3)
    
        if ctrl_rep_number_int_var.get() >= 3: # If this sample replicate exists
            # Display the sample loader button and entry field for the first read of this replicate
            ctrl_rep3_r1_button.grid(column = 1, row = 5, padx = (10,5), pady = (10,0))
            ctrl_rep3_r1_entry.grid(column = 2, row = 5, padx = (5,10), ipady = 3, pady = (10,0))
            
            if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
                # Also display the sample loader button and entry field for the second read of this replicate
                ctrl_rep3_r2_button.grid(column = 1, row = 6, padx = (10,5))
                ctrl_rep3_r2_entry.grid(column = 2, row = 6, padx = (5,10), ipady = 3)

        if ctrl_rep_number_int_var.get() >= 4: # If this sample replicate exists
            # Display the sample loader button and entry field for the first read of this replicate
            ctrl_rep4_r1_button.grid(column = 1, row = 7, padx = (10,5), pady = (10,0))
            ctrl_rep4_r1_entry.grid(column = 2, row = 7, padx = (5,10), ipady = 3, pady = (10,0))

            if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
                # Also display the sample loader button and entry field for the second read of this replicate
                ctrl_rep4_r2_button.grid(column = 1, row = 8, padx = (10,5))
                ctrl_rep4_r2_entry.grid(column = 2, row = 8, padx = (5,10), ipady = 3)

        if ctrl_rep_number_int_var.get() >= 5: # If this sample replicate exists
            # Display the sample loader button and entry field for the first read of this replicate
            ctrl_rep5_r1_button.grid(column = 1, row = 9, padx = (10,5), pady = (10,0))
            ctrl_rep5_r1_entry.grid(column = 2, row = 9, padx = (5,10), ipady = 3, pady = (10,0))

            if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
                # Also display the sample loader button and entry field for the second read of this replicate
                ctrl_rep5_r2_button.grid(column = 1, row = 10, padx = (10,5))
                ctrl_rep5_r2_entry.grid(column = 2, row = 10, padx = (5,10), ipady = 3)

    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path


ctrl_list_frame = tk.Frame(root) # Make a new frame

# Simple label with static text
ctrl_list_label = tk.Label(ctrl_list_frame, text = "Please type in, or browse and select your control sample files below", justify = tk.LEFT, width = 70)
ctrl_list_label.grid(row = 0, column = 1, padx = 10, pady = (5,10), columnspan = 2)

# Clickable button that triggers a function when activated and an entry field where user can manually key in their input to be stored as a variable
# For each of the five (maximum number that can be handled by ChIP-AP) first reads of the control sample replicates, and each of the five second reads 
ctrl_rep1_r1_button = tk.Button(ctrl_list_frame, text ='ctrl rep 1 read 1', command = lambda : ctrl_rep1_r1_button_function(), width = 20)
ctrl_rep1_r1_entry = tk.Entry(ctrl_list_frame, textvariable = ctrl_rep1_r1_string_var, width = 50, justify = tk.RIGHT)
ctrl_rep1_r2_button = tk.Button(ctrl_list_frame, text ='ctrl rep 1 read 2', command = lambda : ctrl_rep1_r2_button_function(), width = 20)
ctrl_rep1_r2_entry = tk.Entry(ctrl_list_frame, textvariable = ctrl_rep1_r2_string_var, width = 50, justify = tk.RIGHT)
ctrl_rep2_r1_button = tk.Button(ctrl_list_frame, text ='ctrl rep 2 read 1', command = lambda : ctrl_rep2_r1_button_function(), width = 20)
ctrl_rep2_r1_entry = tk.Entry(ctrl_list_frame, textvariable = ctrl_rep2_r1_string_var, width = 50, justify = tk.RIGHT)
ctrl_rep2_r2_button = tk.Button(ctrl_list_frame, text ='ctrl rep 2 read 2', command = lambda : ctrl_rep2_r2_button_function(), width = 20)
ctrl_rep2_r2_entry = tk.Entry(ctrl_list_frame, textvariable = ctrl_rep2_r2_string_var, width = 50, justify = tk.RIGHT)
ctrl_rep3_r1_button = tk.Button(ctrl_list_frame, text ='ctrl rep 3 read 1', command = lambda : ctrl_rep3_r1_button_function(), width = 20)
ctrl_rep3_r1_entry = tk.Entry(ctrl_list_frame, textvariable = ctrl_rep3_r1_string_var, width = 50, justify = tk.RIGHT)
ctrl_rep3_r2_button = tk.Button(ctrl_list_frame, text ='ctrl rep 3 read 2', command = lambda : ctrl_rep3_r2_button_function(), width = 20)
ctrl_rep3_r2_entry = tk.Entry(ctrl_list_frame, textvariable = ctrl_rep3_r2_string_var, width = 50, justify = tk.RIGHT)
ctrl_rep4_r1_button = tk.Button(ctrl_list_frame, text ='ctrl rep 4 read 1', command = lambda : ctrl_rep4_r1_button_function(), width = 20)
ctrl_rep4_r1_entry = tk.Entry(ctrl_list_frame, textvariable = ctrl_rep4_r1_string_var, width = 50, justify = tk.RIGHT)
ctrl_rep4_r2_button = tk.Button(ctrl_list_frame, text ='ctrl rep 4 read 2', command = lambda : ctrl_rep4_r2_button_function(), width = 20)
ctrl_rep4_r2_entry = tk.Entry(ctrl_list_frame, textvariable = ctrl_rep4_r2_string_var, width = 50, justify = tk.RIGHT)
ctrl_rep5_r1_button = tk.Button(ctrl_list_frame, text ='ctrl rep 5 read 1', command = lambda : ctrl_rep5_r1_button_function(), width = 20)
ctrl_rep5_r1_entry = tk.Entry(ctrl_list_frame, textvariable = ctrl_rep5_r1_string_var, width = 50, justify = tk.RIGHT)
ctrl_rep5_r2_button = tk.Button(ctrl_list_frame, text ='ctrl rep 5 read 2', command = lambda : ctrl_rep5_r2_button_function(), width = 20)
ctrl_rep5_r2_entry = tk.Entry(ctrl_list_frame, textvariable = ctrl_rep5_r2_string_var, width = 50, justify = tk.RIGHT)

# The help tooltips for each of the ten buttons above
CreateToolTip(ctrl_rep1_r1_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam\nfile extension accepted)')
CreateToolTip(ctrl_rep1_r2_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam\nfile extension accepted)')
CreateToolTip(ctrl_rep2_r1_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam\nfile extension accepted)')
CreateToolTip(ctrl_rep2_r2_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam\nfile extension accepted)')
CreateToolTip(ctrl_rep3_r1_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam\nfile extension accepted)')
CreateToolTip(ctrl_rep3_r2_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam\nfile extension accepted)')
CreateToolTip(ctrl_rep4_r1_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam\nfile extension accepted)')
CreateToolTip(ctrl_rep4_r2_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam\nfile extension accepted)')
CreateToolTip(ctrl_rep5_r1_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam\nfile extension accepted)')
CreateToolTip(ctrl_rep5_r2_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam\nfile extension accepted)')

# Listing of all the control samples-related widgets
ctrl_widget_list = [ctrl_rep1_r1_button, 
                    ctrl_rep1_r2_button, 
                    ctrl_rep2_r1_button, 
                    ctrl_rep2_r2_button, 
                    ctrl_rep3_r1_button, 
                    ctrl_rep3_r2_button, 
                    ctrl_rep4_r1_button, 
                    ctrl_rep4_r2_button, 
                    ctrl_rep5_r1_button, 
                    ctrl_rep5_r2_button, 
                    ctrl_rep1_r1_entry, 
                    ctrl_rep1_r2_entry, 
                    ctrl_rep2_r1_entry, 
                    ctrl_rep2_r2_entry, 
                    ctrl_rep3_r1_entry, 
                    ctrl_rep3_r2_entry, 
                    ctrl_rep4_r1_entry, 
                    ctrl_rep4_r2_entry, 
                    ctrl_rep5_r1_entry, 
                    ctrl_rep5_r2_entry]

# Dynamic label with changeable text
ctrl_sample_notification_label = tk.Label(ctrl_list_frame, textvariable = ctrl_sample_notification_string_var, width = 70, padx = 5, pady = 5)
ctrl_sample_notification_label.grid(row = 15, column = 1, padx = 10, pady = 5, columnspan = 2)

# Create navigation buttons
ctrl_list_frame_back_button = tk.Button(ctrl_list_frame, text = "<< Back", command = lambda : change_frame_function(ctrl_list_frame, chip_list_frame), width = 20)
ctrl_list_frame_continue_button = tk.Button(ctrl_list_frame, text = "Continue >>", command = lambda : change_frame_function(ctrl_list_frame, setting_table_frame), width = 20)
ctrl_list_frame_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
ctrl_list_frame_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################


# Function triggered by the "Load setting table" button, to open a browser window for choosing a local file
def setting_table_button_function():
    setting_table_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), mode = 'rb', title = 'Choose a setting table file')
    if setting_table_button_input: # If the user choose a file from the browser
        setting_table_string_var.set(setting_table_button_input.name) # Set the setting table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(setting_table_button_input.name)) # Set the last visited directory to the user-chosen file directory

    else: # If the user closed the browser without choosing any file from the browser
        clear_setting_function('withfile') # Clear all setting-related variables
        setting_table_question_string_var.set('no') # Return the chosen option into using a setting table
        setting_table_entry_popup_function('no') # Trigger the function where the chosen option is using a setting table i.e., restore all settings to default values
        setting_table_notification_string_var.set('No custom settings table was selected, reverted to default') # Notify the user
        setting_table_notification_label.config(fg = 'red')
        rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path
        return


# Function to test whether the setting table is valid and its contents can be loaded properly without any error
def setting_table_loading_test_function():
    try: # Try to load the setting table with pd.read_csv
        setting_table_absolute_path = setting_table_string_var.get()
        setting_table_df = pd.read_csv(setting_table_absolute_path, delimiter='\t')
        setting_table_df.fillna('', inplace = True) # Replace with empty strings all 'NaN' fields resulting from empty fields in the setting table
        setting_table_header = setting_table_df.columns.values.tolist() # Get the headers name to be parsed below

        # Check the formatting of the custom settings table, to ensure correct program-argument readings.
        if len(setting_table_header) != 2: # If the table does not consist of two columns 
            clear_setting_function('withfile') # Clear all setting-related variables
            setting_table_question_string_var.set('no') # Return the chosen option into using a setting table
            setting_table_entry_popup_function('no') # Trigger the function where the chosen option is using a setting table i.e., restore all settings to default values
            setting_table_notification_string_var.set('Columns error, reverted to default') # Notify the user
            setting_table_notification_label.config(fg = 'red')
            rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path
            return False

        # Check first if they are both strings to avoid TypeError, then check if the table headers are 'program' and 'argument'
        if isinstance(setting_table_header[0], str) and isinstance(setting_table_header[1], str): # If both header values are string type
            if setting_table_header[0].strip().lower() != 'program' or setting_table_header[1].strip().lower() != 'argument': # If the header values are not 'program' and 'argument'
                clear_setting_function('withfile') # Clear all setting-related variables
                setting_table_question_string_var.set('no') # Return the chosen option into using a setting table
                setting_table_entry_popup_function('no') # Trigger the function where the chosen option is using a setting table i.e., restore all settings to default values
                setting_table_notification_string_var.set('Header error, reverted to default') # Notify the user
                setting_table_notification_label.config(fg = 'red')
                rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path
                return False

            else: # If the header values are 'program' and 'argument'
                if setting_table_string_var.get() == default_setting_table_file_full_path: # If the currently registered path in the GUI variable is the path to the default settings table
                    setting_table_notification_string_var.set('Currently using default settings table') # Notify the user
                    setting_table_notification_label.config(fg = 'blue')
                    read_setting_table_function(setting_table_string_var.get()) # Read the registered settings table and load the contained values to their respective GUI variables
                    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path
                    return True

                else: # If the currently registered path in the GUI variable is the path to user-provided settings table
                    setting_table_notification_string_var.set('Custom settings table loading successful') # Notify the user
                    setting_table_notification_label.config(fg = 'green')
                    read_setting_table_function(setting_table_string_var.get()) # Read the registered settings table and load the contained values to their respective GUI variables
                    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path
                    return True

        else: # If both header values are not even string type
            clear_setting_function('withfile') # Clear all setting-related variables
            setting_table_question_string_var.set('no') # Return the chosen option into using a setting table
            setting_table_entry_popup_function('no') # Trigger the function where the chosen option is using a setting table i.e., restore all settings to default values
            setting_table_notification_string_var.set('Header error, reverted to default') # Notify the user
            setting_table_notification_label.config(fg = 'red')
            rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path
            return False

    except: # If setting table loading fails
        clear_setting_function('withfile') # Clear all setting-related variables
        setting_table_question_string_var.set('no') # Return the chosen option into using a setting table
        setting_table_entry_popup_function('no') # Trigger the function where the chosen option is using a setting table i.e., restore all settings to default values
        setting_table_notification_string_var.set('Custom settings table loading error, reverted to default') # Notify the user
        setting_table_notification_label.config(fg = 'red')
        rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path
        return False


# Function to read the registered settings table and load the contained values to their respective GUI variables
def read_setting_table_function(read_setting_table_arg):
    setting_table_absolute_path = read_setting_table_arg
    setting_table_df = pd.read_csv(setting_table_absolute_path, delimiter='\t') # Try to load the setting table with pd.read_csv
    setting_table_df.fillna('', inplace = True) # Replace with empty strings all 'NaN' fields resulting from empty fields in the setting table
    setting_table_header = setting_table_df.columns.values.tolist() # Get the headers name to be parsed below

    # Parse the location of 'program' and 'argument' columns
    setting_table_program_colnum    = setting_table_header.index('program')
    setting_table_argument_colnum   = setting_table_header.index('argument')
    setting_table_array             = setting_table_df.values.tolist()

    # Prepare a global dictionary to store the program names with their respective registered arguments
    global argument_dict
    argument_dict = {}

    # Register all the program names first into the dictionary, with empty strings as the argument values
    for suite_program in suite_program_list:
        argument_dict[suite_program] = []

    # For each entry in the 'program' column that matched with any of the program name in suite_program_list,
    #   assign the argument as the value and program as the key in dictionary argument_dict
    for setting_table_array_row in setting_table_array:
        if setting_table_array_row[setting_table_program_colnum] in suite_program_list:
            current_setting_table_program = setting_table_array_row[setting_table_program_colnum]
            current_setting_table_argument = setting_table_array_row[setting_table_argument_colnum]
            argument_dict[current_setting_table_program].append(current_setting_table_argument)

    # Join all arguments value within each program key with a single space, 
    #   and assign the joined string into their own variable for easier calling later downstream
    fastqc1_arg.set(' '.join(argument_dict['fastqc1']))
    clumpify_arg.set(' '.join(argument_dict['clumpify']))
    bbduk_arg.set(' '.join(argument_dict['bbduk']))
    trimmomatic_arg.set(' '.join(argument_dict['trimmomatic']))
    fastqc2_arg.set(' '.join(argument_dict['fastqc2']))
    bwa_mem_arg.set(' '.join(argument_dict['bwa_mem']))
    samtools_view_arg.set(' '.join(argument_dict['samtools_view']))
    plotfingerprint_arg.set(' '.join(argument_dict['plotfingerprint']))
    fastqc3_arg.set(' '.join(argument_dict['fastqc3']))
    macs2_callpeak_arg.set(' '.join(argument_dict['macs2_callpeak']))
    gem_arg.set(' '.join(argument_dict['gem']))
    sicer2_arg.set(' '.join(argument_dict['sicer2']))
    homer_findPeaks_arg.set(' '.join(argument_dict['homer_findPeaks']))
    genrich_arg.set(' '.join(argument_dict['genrich']))
    homer_mergePeaks_arg.set(' '.join(argument_dict['homer_mergePeaks']))
    homer_annotatePeaks_arg.set(' '.join(argument_dict['homer_annotatePeaks']))
    fold_change_calculator_arg.set(' '.join(argument_dict['fold_change_calculator']))
    homer_findMotifsGenome_arg.set(' '.join(argument_dict['homer_findMotifsGenome']))
    meme_chip_arg.set(' '.join(argument_dict['meme_chip']))

# Function that manages the widgets in the GUI's setting table loading frame
# Toggled by the user choice between loading a setting table and using default settings
def setting_table_entry_popup_function(setting_table_entry_popup_arg):
    if setting_table_entry_popup_arg == 'yes': # If the user choose to load a setting table
        # Display the buttons and entry field for loading the setting table
        setting_table_entry.grid(row = 3, column = 1, padx = (5,10), ipady = 3)
        setting_table_button.grid(row = 3, column = 2, padx = (10,5), pady = 2)
        clear_setting_function('withfile') # Clear all sample-related variables
        setting_table_notification_string_var.set('Please load your custom settings table file') # Prompt the user
        setting_table_notification_label.config(fg = 'blue')
        # Functions to keep the frame window (and its contents) at the center of screen
        root.eval('tk::PlaceWindow . center')
        root.grab_set()
        root.focus_force()

    if setting_table_entry_popup_arg == 'no': # If the user choose to use the default settings
        # Remove the buttons and entry field for loading the setting table
        setting_table_entry.grid_remove()
        setting_table_button.grid_remove()
        setting_table_question_string_var.set('no') # Return the chosen option into using a setting table
        setting_table_string_var.set(default_setting_table_file_full_path) # Set the setting table GUI variable with the path to the default settings table file
        read_setting_table_function(default_setting_table_file_full_path) # Read the default settings table and load the contained values to their respective GUI variables
        setting_table_notification_string_var.set('Currently using default settings table') # Prompt the user
        setting_table_notification_label.config(fg = 'blue')
        # Functions to keep the frame window (and its contents) at the center of screen
        root.eval('tk::PlaceWindow . center')
        root.grab_set()
        root.focus_force()

    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path


# Function to clear all setting-related variables in the GUI
def clear_setting_function(clear_setting_function_arg):

    if clear_setting_function_arg == 'withfile': # When the user wants to also clear the registered path to the setting table file
        setting_table_string_var.set('') # Register empty strings to the variable that stores the path to the setting table file

    # Register empty strings to all setting-related GUI variables
    fastqc1_arg.set('')
    clumpify_arg.set('')
    bbduk_arg.set('')
    trimmomatic_arg.set('')
    fastqc2_arg.set('')
    bwa_mem_arg.set('')
    samtools_view_arg.set('')
    plotfingerprint_arg.set('')
    fastqc3_arg.set('')
    macs2_callpeak_arg.set('')
    gem_arg.set('')
    sicer2_arg.set('')
    homer_findPeaks_arg.set('')
    genrich_arg.set('')
    homer_mergePeaks_arg.set('')
    homer_annotatePeaks_arg.set('')
    fold_change_calculator_arg.set('')
    homer_findMotifsGenome_arg.set('')
    meme_chip_arg.set('')

    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path


setting_table_frame = tk.Frame(root) # Make a new frame

# Simple label with static text
setting_table_label = tk.Label(setting_table_frame, text = "Do you want to use your custom settings table?", width = 70)
setting_table_label.grid(row = 1, column = 1, padx = 10, pady = 5, columnspan = 2)

# Radio button of the first choice that triggers a function when activated
setting_table_yes_radio = tk.Radiobutton(setting_table_frame, text = "Yes, I would like use my custom settings table", padx = 5, variable = setting_table_question_string_var, value = 'yes', command = lambda : setting_table_entry_popup_function('yes'), width = 45)
CreateToolTip(setting_table_yes_radio, text = 'Select this to load a custom\nsetting table in ChIP-AP format\n(see GitHub documentation).\nThe ChIP-AP pipeline will then be\nrun based on the custom settings\nfrom loaded setting table')
setting_table_yes_radio.grid(row = 2, column = 1, padx = 10, pady = 2, columnspan = 2)

# Entry field where user can manually key in their input to be stored as a variable
setting_table_entry = tk.Entry(setting_table_frame, textvariable = setting_table_string_var, width = 50, justify = tk.RIGHT)
setting_table_button = tk.Button(setting_table_frame, text = 'Browse', command = lambda : setting_table_button_function(), width = 20)
CreateToolTip(setting_table_button, text = 'Click here to browse and select\nyour setting table file')

# Radio button of the second choice that triggers a function when activated
setting_table_no_radio = tk.Radiobutton(setting_table_frame, text = "No, I would like to use the default settings", padx = 5, variable = setting_table_question_string_var, value = 'no', command = lambda : setting_table_entry_popup_function('no'), width = 45)
CreateToolTip(setting_table_no_radio, text = 'Select this if you want to use\nthe ChIP-AP pipeline default settings')
setting_table_no_radio.grid(row = 4, column = 1, padx = 10, pady = 2, columnspan = 2)

# Dynamic label with changeable text
setting_table_notification_label = tk.Label(setting_table_frame, textvariable = setting_table_notification_string_var, width = 70, padx = 5, pady = 5)
setting_table_notification_label.grid(row = 5, column = 1, padx = 10, pady = 5, columnspan = 2)

# Create navigation buttons
setting_table_back_button = tk.Button(setting_table_frame, text = "<< Back", command = lambda : change_frame_function(setting_table_frame, ctrl_list_frame), width = 20)
setting_table_customize_button = tk.Button(setting_table_frame, text = "Customize settings", command = lambda : change_frame_function(setting_table_frame, setting_value_frame), width = 20)
setting_table_continue_button = tk.Button(setting_table_frame, text = "Continue >>", command = lambda : change_frame_function(setting_table_frame, genome_ref_frame), width = 20)
setting_table_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
setting_table_customize_button.grid(row = 21, column = 1, columnspan = 2, pady = (5,10))
setting_table_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################


setting_value_frame = tk.Frame(root) # Make a new frame

# Simple label with static text
setting_value_label = tk.Label(setting_value_frame, text = "Please modify the setting values below as needed", width = 70)
CreateToolTip(setting_value_label, text = "Warning: proceed only when\nyou know what you are doing.\nOtherwise, leave all settings\nat their default values.\nInvalid values will\nlikely cause the whole\npipeline run to break.\nCheck GitHub documentation\nto learn more about\ncustom settings.")
setting_value_label.grid(row = 0, column = 1, pady = 5, columnspan = 2)

# Simple labels with static texts
fastqc1_arg_label = tk.Label(setting_value_frame, text = 'fastqc1', width = 22, pady = 5, anchor = 'e')
clumpify_arg_label = tk.Label(setting_value_frame, text = 'clumpify', width = 22, pady = 5, anchor = 'e')
bbduk_arg_label = tk.Label(setting_value_frame, text = 'bbduk', width = 22, pady = 5, anchor = 'e')
trimmomatic_arg_label = tk.Label(setting_value_frame, text = 'trimmomatic', width = 22, pady = 5, anchor = 'e')
fastqc2_arg_label = tk.Label(setting_value_frame, text = 'fastqc2', width = 22, pady = 5, anchor = 'e')
bwa_mem_arg_label = tk.Label(setting_value_frame, text = 'bwa_mem', width = 22, pady = 5, anchor = 'e')
samtools_view_arg_label = tk.Label(setting_value_frame, text = 'samtools_view', width = 22, pady = 5, anchor = 'e')
plotfingerprint_arg_label = tk.Label(setting_value_frame, text = 'plotfingerprint', width = 22, pady = 5, anchor = 'e')
fastqc3_arg_label = tk.Label(setting_value_frame, text = 'fastqc3', width = 22, pady = 5, anchor = 'e')
macs2_callpeak_arg_label = tk.Label(setting_value_frame, text = 'macs2_callpeak', width = 22, pady = 5, anchor = 'e')
gem_arg_label = tk.Label(setting_value_frame, text = 'gem', width = 22, pady = 5, anchor = 'e')
sicer2_arg_label = tk.Label(setting_value_frame, text = 'sicer2', width = 22, pady = 5, anchor = 'e')
homer_findPeaks_arg_label = tk.Label(setting_value_frame, text = 'homer_findPeaks', width = 22, pady = 5, anchor = 'e')
genrich_arg_label = tk.Label(setting_value_frame, text = 'genrich', width = 22, pady = 5, anchor = 'e')
homer_mergePeaks_arg_label = tk.Label(setting_value_frame, text = 'homer_mergePeaks', width = 22, pady = 5, anchor = 'e')
homer_annotatePeaks_arg_label = tk.Label(setting_value_frame, text = 'homer_annotatePeaks', width = 22, pady = 5, anchor = 'e')
fold_change_calculator_arg_label = tk.Label(setting_value_frame, text = 'fold_change_calculator', width = 22, pady = 5, anchor = 'e')
homer_findMotifsGenome_label = tk.Label(setting_value_frame, text = 'homer_findMotifsGenome', width = 22, pady = 5, anchor = 'e')
meme_chip_label = tk.Label(setting_value_frame, text = 'meme_chip', width = 22, pady = 5, anchor = 'e')

# Display all the labels created above
fastqc1_arg_label.grid(row = 1, column = 1, padx = 5)
clumpify_arg_label.grid(row = 2, column = 1, padx = 5)
bbduk_arg_label.grid(row = 3, column = 1, padx = 5)
trimmomatic_arg_label.grid(row = 4, column = 1, padx = 5)
fastqc2_arg_label.grid(row = 5, column = 1, padx = 5)
bwa_mem_arg_label.grid(row = 6, column = 1, padx = 5)
samtools_view_arg_label.grid(row = 7, column = 1, padx = 5)
plotfingerprint_arg_label.grid(row = 8, column = 1, padx = 5)
fastqc3_arg_label.grid(row = 9, column = 1, padx = 5)
macs2_callpeak_arg_label.grid(row = 10, column = 1, padx = 5)
gem_arg_label.grid(row = 11, column = 1, padx = 5)
sicer2_arg_label.grid(row = 12, column = 1, padx = 5)
homer_findPeaks_arg_label.grid(row = 13, column = 1, padx = 5)
genrich_arg_label.grid(row = 14, column = 1, padx = 5)
homer_mergePeaks_arg_label.grid(row = 15, column = 1, padx = 5)
homer_annotatePeaks_arg_label.grid(row = 16, column = 1, padx = 5)
fold_change_calculator_arg_label.grid(row = 17, column = 1, padx = 5)
homer_findMotifsGenome_label.grid(row = 18, column = 1, padx = 5)
meme_chip_label.grid(row = 19, column = 1, padx = 5)

# Entry fields where user can manually key in their input to be stored as a variable
fastqc1_arg_entry = tk.Entry(setting_value_frame, textvariable = fastqc1_arg, width = 50, justify = tk.LEFT)
clumpify_arg_entry = tk.Entry(setting_value_frame, textvariable = clumpify_arg, width = 50, justify = tk.LEFT)
bbduk_arg_entry = tk.Entry(setting_value_frame, textvariable = bbduk_arg, width = 50, justify = tk.LEFT)
trimmomatic_arg_entry = tk.Entry(setting_value_frame, textvariable = trimmomatic_arg, width = 50, justify = tk.LEFT)
fastqc2_arg_entry = tk.Entry(setting_value_frame, textvariable = fastqc2_arg, width = 50, justify = tk.LEFT)
bwa_mem_arg_entry = tk.Entry(setting_value_frame, textvariable = bwa_mem_arg, width = 50, justify = tk.LEFT)
samtools_view_arg_entry = tk.Entry(setting_value_frame, textvariable = samtools_view_arg, width = 50, justify = tk.LEFT)
plotfingerprint_arg_entry = tk.Entry(setting_value_frame, textvariable = plotfingerprint_arg, width = 50, justify = tk.LEFT)
fastqc3_arg_entry = tk.Entry(setting_value_frame, textvariable = fastqc3_arg, width = 50, justify = tk.LEFT)
macs2_callpeak_arg_entry = tk.Entry(setting_value_frame, textvariable = macs2_callpeak_arg, width = 50, justify = tk.LEFT)
gem_arg_entry = tk.Entry(setting_value_frame, textvariable = gem_arg, width = 50, justify = tk.LEFT)
sicer2_arg_entry = tk.Entry(setting_value_frame, textvariable = sicer2_arg, width = 50, justify = tk.LEFT)
homer_findPeaks_arg_entry = tk.Entry(setting_value_frame, textvariable = homer_findPeaks_arg, width = 50, justify = tk.LEFT)
genrich_arg_entry = tk.Entry(setting_value_frame, textvariable = genrich_arg, width = 50, justify = tk.LEFT)
homer_mergePeaks_arg_entry = tk.Entry(setting_value_frame, textvariable = homer_mergePeaks_arg, width = 50, justify = tk.LEFT)
homer_annotatePeaks_arg_entry = tk.Entry(setting_value_frame, textvariable = homer_annotatePeaks_arg, width = 50, justify = tk.LEFT)
fold_change_calculator_arg_entry = tk.Entry(setting_value_frame, textvariable = fold_change_calculator_arg, width = 50, justify = tk.LEFT)
homer_findMotifsGenome_arg_entry = tk.Entry(setting_value_frame, textvariable = homer_findMotifsGenome_arg, width = 50, justify = tk.LEFT)
meme_chip_arg_entry = tk.Entry(setting_value_frame, textvariable = meme_chip_arg, width = 50, justify = tk.LEFT)

# Display all the entry fields created above
fastqc1_arg_entry.grid(row = 1, column = 2, padx = 5, ipady = 3)
clumpify_arg_entry.grid(row = 2, column = 2, padx = 5, ipady = 3)
bbduk_arg_entry.grid(row = 3, column = 2, padx = 5, ipady = 3)
trimmomatic_arg_entry.grid(row = 4, column = 2, padx = 5, ipady = 3)
fastqc2_arg_entry.grid(row = 5, column = 2, padx = 5, ipady = 3)
bwa_mem_arg_entry.grid(row = 6, column = 2, padx = 5, ipady = 3)
samtools_view_arg_entry.grid(row = 7, column = 2, padx = 5, ipady = 3)
plotfingerprint_arg_entry.grid(row = 8, column = 2, padx = 5, ipady = 3)
fastqc3_arg_entry.grid(row = 9, column = 2, padx = 5, ipady = 3)
macs2_callpeak_arg_entry.grid(row = 10, column = 2, padx = 5, ipady = 3)
gem_arg_entry.grid(row = 11, column = 2, padx = 5, ipady = 3)
sicer2_arg_entry.grid(row = 12, column = 2, padx = 5, ipady = 3)
homer_findPeaks_arg_entry.grid(row = 13, column = 2, padx = 5, ipady = 3)
genrich_arg_entry.grid(row = 14, column = 2, padx = 5, ipady = 3)
homer_mergePeaks_arg_entry.grid(row = 15, column = 2, padx = 5, ipady = 3)
homer_annotatePeaks_arg_entry.grid(row = 16, column = 2, padx = 5, ipady = 3)
fold_change_calculator_arg_entry.grid(row = 17, column = 2, padx = 5, ipady = 3)
homer_findMotifsGenome_arg_entry.grid(row = 18, column = 2, padx = 5, ipady = 3)
meme_chip_arg_entry.grid(row = 19, column = 2, padx = 5, ipady = 3)

# Clickable button that triggers a function when activated
default_setting_button = tk.Button(setting_value_frame, text = "Restore defaults", command = lambda : setting_table_entry_popup_function('no'), width = 22)
CreateToolTip(default_setting_button, text = 'Click here if you want to restore\nto ChIP-AP pipeline default settings')
default_setting_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))

# Create navigation buttons
setting_value_accept_button = tk.Button(setting_value_frame, text = "Accept and close", command = lambda : change_frame_function(setting_value_frame, setting_table_frame), width = 22)
setting_value_accept_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################


# Function to auto-assign the genome directory to the ChIP-AP's default when the user-chosen reference genome is one of the supported builds
# Opens up the button and entry field widgets for inputting by user only when user chooses "others" reference genome, which genome files needs to be self-provided by user  
def auto_assign_genome_folder_function():
    if genome_ref_string_var.get() == "hg38 (Homo sapiens build 38)":
        genome_folder_string_var.set(genome_folder_full_path)
        genome_folder_button.grid_remove()
        genome_folder_entry.grid_remove()
    elif genome_ref_string_var.get() == "hg19 (Homo sapiens build 19)":
        genome_folder_string_var.set(genome_folder_full_path)
        genome_folder_button.grid_remove()
        genome_folder_entry.grid_remove()
    elif genome_ref_string_var.get() == "mm9 (Mus musculus build 9)":
        genome_folder_string_var.set(genome_folder_full_path)
        genome_folder_button.grid_remove()
        genome_folder_entry.grid_remove()
    elif genome_ref_string_var.get() == "mm10 (Mus musculus build 10)":
        genome_folder_string_var.set(genome_folder_full_path)
        genome_folder_button.grid_remove()
        genome_folder_entry.grid_remove()
    elif genome_ref_string_var.get() == "dm6 (Drosophila melanogaster build 6)":
        genome_folder_string_var.set(genome_folder_full_path)
        genome_folder_button.grid_remove()
        genome_folder_entry.grid_remove()
    elif genome_ref_string_var.get() == "sacCer3 (Saccharomyces cerevisiae build 3)":
        genome_folder_string_var.set(genome_folder_full_path)
        genome_folder_button.grid_remove()
        genome_folder_entry.grid_remove()
    else:
        genome_folder_button.grid(row = 2, column = 1, padx = (10,5), pady = 5)
        genome_folder_entry.grid(row = 2, column = 2, padx = (5,10), pady = 5, ipady = 3)

    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path


# Function triggered by the "Genome folder" button, to open a browser window for choosing a local directory
def genome_folder_button_function():
    genome_folder_button_input = filedialog.askdirectory(initialdir = current_dir_string_var.get(), parent = root, title = 'Choose the Genome Directory')
    if genome_folder_button_input: # If the user choose a directory from the browser
        genome_folder_string_var.set(genome_folder_button_input) # Set the setting table GUI variable with the path to the user-chosen directory
        current_dir_string_var.set(genome_folder_button_input) # Set the last visited directory to the user-chosen file directory

    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path


genome_ref_frame = tk.Frame(root) # Make a new frame

# Simple label with static text
genome_ref_label = tk.Label(genome_ref_frame, text = "Reference genome:")
genome_ref_label.grid(row = 1, column = 1, padx = (10,5), pady = (10,5))

# Drop-down options widget that triggers a function when an option is chosen
genome_ref_drop_down = tk.OptionMenu(genome_ref_frame, genome_ref_string_var, *genome_ref_options, command = lambda x = None : auto_assign_genome_folder_function())
CreateToolTip(genome_ref_drop_down, text = 'Select the genome assembly you want\nyour ChIP-seq reads to be aligned to.\nCurrently, ChIP-AP supports\nsix genome assemblies.\nOutside those supported by\nChIP-AP, you will need to\ngenerate the files yourself')
genome_ref_drop_down.config(width = 45, takefocus = 1)
genome_ref_drop_down_menu = root.nametowidget(genome_ref_drop_down.menuname) # Get the drop-down menu options to follow the default fonts
genome_ref_drop_down_menu.config(font = default_font) # Get the drop-down menu options to follow the default fonts 
genome_ref_drop_down.grid(row = 1, column = 2, padx = (5,10), pady = (10,5))

# Clickable button that triggers a function when activated
genome_folder_button = tk.Button(genome_ref_frame, text = 'Genome folder', command = lambda : genome_folder_button_function(), state = tk.NORMAL, width = 20)
CreateToolTip(genome_folder_button, text = 'Click here to browse and select\nthe directory containing your\ncustom genome reference')
genome_folder_entry = tk.Entry(genome_ref_frame, textvariable = genome_folder_string_var, width = 50, justify = tk.RIGHT)

# Create navigation buttons
genome_ref_back_button = tk.Button(genome_ref_frame, text = "<< Back", command = lambda : change_frame_function(genome_ref_frame, setting_table_frame), width = 20)
genome_ref_continue_button = tk.Button(genome_ref_frame, text = "Continue >>", command = lambda : change_frame_function(genome_ref_frame, known_motif_frame), width = 20)
genome_ref_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
genome_ref_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################


# Function that manages the widgets in the GUI's known motif file loading frame
# Toggled by the user choice between loading and not loading a known motif file
def known_motif_options_popup_function(known_motif_options_popup_arg):
    if known_motif_options_popup_arg == 'yes': # If the user choose to load a known motif file
        # Display the buttons and entry field for loading the known motif file
        known_motif_button.grid(sticky = "e", row = 3, column = 1, pady = (0,10), padx = (10,5))
        known_motif_entry.grid(sticky = "w", row = 3, column = 2, pady = (0,10), padx = (5,10), ipady = 3)
        # Functions to keep the frame window (and its contents) at the center of screen
        root.eval('tk::PlaceWindow . center')
        root.grab_set()
        root.focus_force()

    if known_motif_options_popup_arg == 'no': # If the user choose to not load a known motif file
        # Remove the buttons and entry field for loading the known motif file
        known_motif_button.grid_remove()
        known_motif_entry.grid_remove()
        known_motif_string_var.set('') # Clear the known motif file GUI variable
        # Functions to keep the frame window (and its contents) at the center of screen
        root.eval('tk::PlaceWindow . center')
        root.grab_set()
        root.focus_force()

    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path


# Function triggered by the "Known motif file" button, to open a browser window for choosing a local file
def known_motif_button_function():
    known_motif_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a .motif file')
    if known_motif_button_input: # If the user choose a file from the browser
        known_motif_string_var.set(known_motif_button_input.name) # Set the setting table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(known_motif_button_input.name)) # Set the last visited directory to the user-chosen file directory

    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path


known_motif_frame = tk.Frame(root) # Make a new frame

# Simple label with static text
known_motif_label = tk.Label(known_motif_frame, text = "Do you want to load a known DNA-binding motif file?", width = 70)
known_motif_label.grid(row = 1, column = 1, padx = 10, pady = 5, columnspan = 2)

# Radio button of the first choice that triggers a function when activated
known_motif_yes_radio = tk.Radiobutton(known_motif_frame, text = "Yes, I would like to use it in my analysis", padx = 5, variable = known_motif_question_string_var, value = 'yes', command = lambda : known_motif_options_popup_function('yes'), width = 70)
CreateToolTip(known_motif_yes_radio, text = 'Select this if you wish to load\na known DNA-binding motif file\n(in HOMER matrix format)')
known_motif_yes_radio.grid(row = 2, column = 1, padx = 10, pady = 5, columnspan = 2)

# Clickable button that triggers a function when activated
known_motif_button = tk.Button(known_motif_frame, text = 'Known motif file', command = lambda : known_motif_button_function(), state = tk.NORMAL, width = 20)
CreateToolTip(known_motif_button, text = 'Click here to browse and select your\n.motif file (in HOMER matrix format)')

# Entry field where user can manually key in their input to be stored as a variable
known_motif_entry = tk.Entry(known_motif_frame, textvariable = known_motif_string_var, width = 50, justify = tk.RIGHT)

# Radio button of the second choice that triggers a function when activated
known_motif_no_radio = tk.Radiobutton(known_motif_frame, text = "No, that would not be necessary", padx = 5, variable = known_motif_question_string_var, value = 'no', command = lambda : known_motif_options_popup_function('no'), width = 70)
CreateToolTip(known_motif_no_radio, text = 'Select this to do not wish\nto load, or do not have, any\nknown DNA-binding motif file')
known_motif_no_radio.grid(row = 4, column = 1, padx = 10, pady = 5, columnspan = 2)

# Create navigation buttons
known_motif_back_button = tk.Button(known_motif_frame, text = "<< Back", command = lambda : change_frame_function(known_motif_frame, genome_ref_frame), width = 20)
known_motif_continue_button = tk.Button(known_motif_frame, text = "Continue >>", command = lambda : change_frame_function(known_motif_frame, output_folder_frame), width = 20)
known_motif_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
known_motif_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################

# Function triggered by the "Output save folder" button, to open a browser window for choosing a local directory
def output_folder_button_function():
    output_folder_button_input = filedialog.askdirectory(initialdir = current_dir_string_var.get(), parent = root, title = 'Choose the Output Directory')
    if output_folder_button_input: # If the user choose a directory from the browser
        output_folder_string_var.set(output_folder_button_input) # Set the setting table GUI variable with the path to the user-chosen directory
        current_dir_string_var.set(output_folder_button_input) # Set the last visited directory to the user-chosen file directory

    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path


output_folder_frame = tk.Frame(root) # Make a new frame

# Clickable button that triggers a function when activated
output_folder_button = tk.Button(output_folder_frame, text = 'Output save folder', command = lambda : output_folder_button_function(), state = tk.NORMAL, width = 20)
CreateToolTip(output_folder_button, text = 'Click here to browse and select\nthe directory to save the results\nof your ChIP-AP pipeline run.\nTo save in a new folder,\ntype in the entry field\nyour new folder name after\nthe output save directory:\nfull_path_to_output_save_directory/\nnew_folder_name')
output_folder_button.grid(row = 1, column = 1, padx = (10,5), pady = (10,5))

# Entry field where user can manually key in their input to be stored as a variable
output_folder_entry = tk.Entry(output_folder_frame, textvariable = output_folder_string_var, width = 50, justify = tk.RIGHT)
output_folder_entry.grid(row = 1, column = 2, padx = (5,10), pady = (10,5), ipady = 3)

# Simple label with static text
setname_label = tk.Label(output_folder_frame, text = "Dataset prefix:", width = 20)
setname_label.grid(row = 2, column = 1, padx = (10,5), pady = 5)

# Entry field where user can manually key in their input to be stored as a variable
setname_entry = tk.Entry(output_folder_frame, textvariable = setname_string_var, width = 50, justify = tk.LEFT)
CreateToolTip(setname_entry, text = 'Type in your folder name and prefix\nfor all the resulting output filenames')
setname_entry.grid(row = 2, column = 2, padx = (5,10), pady = 5, ipady = 3)

# Create navigation buttons
output_folder_back_button = tk.Button(output_folder_frame, text = "<< Back", command = lambda : change_frame_function(output_folder_frame, known_motif_frame), width = 20)
output_folder_continue_button = tk.Button(output_folder_frame, text = "Continue >>", command = lambda : change_frame_function(output_folder_frame, homer_motif_frame), width = 20)
output_folder_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
output_folder_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################


# Function that manages the widgets in the GUI's motif enrichment analysis options by HOMER frame
# Toggled by the user choice between performing or not performing motif enrichment analysis by HOMER
def homer_motif_options_popup_function(homer_motif_options_popup_arg):
    if homer_motif_options_popup_arg == 'yes': # If the user choose to perform motif enrichment analysis by HOMER
        # Display the label and drop-down list for choosing the peakset to perform motif enrichment analysis on
        homer_motif_label.grid(sticky = "e", row = 3, column = 1, pady = (0,10), padx = (20,5))
        homer_motif_drop_down.grid(sticky = "w", row = 3, column = 2, pady = (0,10), padx = (5,20))
        # Functions to keep the frame window (and its contents) at the center of screen
        root.eval('tk::PlaceWindow . center')
        root.grab_set()
        root.focus_force()

    if homer_motif_options_popup_arg == 'no': # If the user choose not to perform motif enrichment analysis by HOMER
        # Remove the label and drop-down list for choosing the peakset to perform motif enrichment analysis on
        homer_motif_label.grid_remove()
        homer_motif_drop_down.grid_remove()
        homer_motif_string_var.set(homer_motif_options[0]) # Reset to the the topmost option of the drop-down list
        # Functions to keep the frame window (and its contents) at the center of screen
        root.eval('tk::PlaceWindow . center')
        root.grab_set()
        root.focus_force()

    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path


homer_motif_frame = tk.Frame(root) # Make a new frame

# Simple label with static text
homer_motif_label = tk.Label(homer_motif_frame, text = "Do you want HOMER to perform motif enrichment analysis?", width = 70)
homer_motif_label.grid(row = 1, column = 1, padx = 10, pady = 5, columnspan = 2)

# Radio button of the first choice that triggers a function when activated
homer_motif_yes_radio = tk.Radiobutton(homer_motif_frame, text = "Yes, I would like HOMER to perform motif enrichment analysis", padx = 5, variable = homer_motif_question_string_var, value = 'yes', command = lambda : homer_motif_options_popup_function('yes'), width = 70)
CreateToolTip(homer_motif_yes_radio, text = 'Select this if you wish to perform\nmotif enrichment analysis\nwith HOMER findMotifsGenome.pl')
homer_motif_yes_radio.grid(row = 2, column = 1, padx = 10, pady = 5, columnspan = 2)

# Simple label with static text
homer_motif_label = tk.Label(homer_motif_frame, text = "... on peak set:", width = 15)

# Drop-down options widget that triggers a function when an option is chosen
homer_motif_drop_down = tk.OptionMenu(homer_motif_frame, homer_motif_string_var, *homer_motif_options)
CreateToolTip(homer_motif_drop_down, text = 'Select the peak set(s) you want\nHOMER findMotifsGenome to perform\nmotif enrichment analysis on')
homer_motif_drop_down.config(takefocus = 1, width = 30)

# Radio button of the second choice that triggers a function when activated
homer_motif_no_radio = tk.Radiobutton(homer_motif_frame, text = "No, that would not be necessary", padx = 5, variable = homer_motif_question_string_var, value = 'no', command = lambda : homer_motif_options_popup_function('no'), width = 70)
CreateToolTip(homer_motif_no_radio, text = 'Select this to do not wish to \nperform motif enrichment analysis\nwith HOMER findMotifsGenome.pl')
homer_motif_no_radio.grid(row = 4, column = 1, padx = 10, pady = 5, columnspan = 2)

# Create navigation buttons
homer_motif_back_button = tk.Button(homer_motif_frame, text = "<< Back", command = lambda : change_frame_function(homer_motif_frame, output_folder_frame), width = 20)
homer_motif_continue_button = tk.Button(homer_motif_frame, text = "Continue >>", command = lambda : change_frame_function(homer_motif_frame, meme_motif_frame), width = 20)
homer_motif_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
homer_motif_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################


# Function that manages the widgets in the GUI's motif enrichment analysis options by MEME frame
# Toggled by the user choice between performing or not performing motif enrichment analysis by MEME
def meme_motif_options_popup_function(meme_motif_options_popup_arg):
    if meme_motif_options_popup_arg == 'yes': # If the user choose to perform motif enrichment analysis by MEME
        # Display the label and drop-down list for choosing the peakset to perform motif enrichment analysis on
        meme_motif_label.grid(sticky = "e", row = 3, column = 1, pady = (0,10), padx = (20,5))
        meme_motif_drop_down.grid(sticky = "w", row = 3, column = 2, pady = (0,10), padx = (5,20))
        # Functions to keep the frame window (and its contents) at the center of screen
        root.eval('tk::PlaceWindow . center')
        root.grab_set()
        root.focus_force()

    if meme_motif_options_popup_arg == 'no': # If the user choose not to perform motif enrichment analysis by MEME
        # Remove the label and drop-down list for choosing the peakset to perform motif enrichment analysis on
        meme_motif_label.grid_remove()
        meme_motif_drop_down.grid_remove()
        meme_motif_string_var.set(meme_motif_options[0]) # Reset to the the topmost option of the drop-down list
        # Functions to keep the frame window (and its contents) at the center of screen
        root.eval('tk::PlaceWindow . center')
        root.grab_set()
        root.focus_force()

    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path


meme_motif_frame = tk.Frame(root) # Make a new frame

# Simple label with static text
meme_motif_label = tk.Label(meme_motif_frame, text = "Do you want MEME to perform motif enrichment analysis?", width = 70)
meme_motif_label.grid(row = 1, column = 1, padx = 10, pady = 5, columnspan = 2)

# Radio button of the first choice that triggers a function when activated
meme_motif_yes_radio = tk.Radiobutton(meme_motif_frame, text = "Yes, I would like MEME to perform motif enrichment analysis", padx = 5, variable = meme_motif_question_string_var, value = 'yes', command = lambda : meme_motif_options_popup_function('yes'), width = 70)
CreateToolTip(meme_motif_yes_radio, text = 'Select this if you wish to\nperform motif enrichment analysis\nwith meme-chip')
meme_motif_yes_radio.grid(row = 2, column = 1, padx = 10, pady = 5, columnspan = 2)

meme_motif_label = tk.Label(meme_motif_frame, text = "... on peak set:", width = 15)

# Drop-down options widget that triggers a function when an option is chosen
meme_motif_drop_down = tk.OptionMenu(meme_motif_frame, meme_motif_string_var, *meme_motif_options)
CreateToolTip(meme_motif_drop_down, text = 'Select the peak set(s) you\nwant meme-chip to perform\nmotif enrichment analysis on')
meme_motif_drop_down.config(takefocus = 1, width = 30)

# Radio button of the second choice that triggers a function when activated
meme_motif_no_radio = tk.Radiobutton(meme_motif_frame, text = "No, that would not be necessary", padx = 5, variable = meme_motif_question_string_var, value = 'no', command = lambda : meme_motif_options_popup_function('no'), width = 70)
CreateToolTip(meme_motif_no_radio, text = 'Select this to do not wish to\nperform motif enrichment analysis\nwith meme-chip')
meme_motif_no_radio.grid(row = 4, column = 1, padx = 10, pady = 5, columnspan = 2)

# Create navigation buttons
meme_motif_back_button = tk.Button(meme_motif_frame, text = "<< Back", command = lambda : change_frame_function(meme_motif_frame, homer_motif_frame), width = 20)
meme_motif_continue_button = tk.Button(meme_motif_frame, text = "Continue >>", command = lambda : change_frame_function(meme_motif_frame, checkbox_frame), width = 20)
meme_motif_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
meme_motif_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


#######################################################################################################################


checkbox_frame = tk.Frame(root) # Make a new frame

# Checkbox to toggle between on and off state what triggers a function when it is either toggled on or off.
fcmerge_checkbox = tk.Checkbutton(checkbox_frame, text = ' Merged fold enrichment analysis', variable = fcmerge_var, onvalue = 1, offvalue = 0)
CreateToolTip(fcmerge_checkbox, text = 'Check this box if you want\nthe fold enrichment analysis\nfrom all replicates combined as one.\nThis option will be ignored\nwhen there are unequal\nnumber of replicates between\nChIP and control samples')
fcmerge_checkbox.grid(sticky = "w", row = 1, column = 1, padx = 10, columnspan = 2)

# Checkbox to toggle between on and off state what triggers a function when it is either toggled on or off.
goann_checkbox = tk.Checkbutton(checkbox_frame, text = ' Annotate peaks with known gene ontology terms', variable = goann_var, onvalue = 1, offvalue = 0)
CreateToolTip(goann_checkbox, text = 'Check this box if you want\neach peak in the final peaks list\nto have gene ontology annotations')
goann_checkbox.grid(sticky = "w", row = 2, column = 1, padx = 10, columnspan = 2)

# Checkbox to toggle between on and off state what triggers a function when it is either toggled on or off.
pathann_checkbox = tk.Checkbutton(checkbox_frame, text = ' Annotate peaks with known pathway terms', variable = pathann_var, onvalue = 1, offvalue = 0)
CreateToolTip(pathann_checkbox, text = 'Check this box if you want\neach peak in the final peaks list\nto have pathway annotations')
pathann_checkbox.grid(sticky = "w", row = 3, column = 1, padx = 10, columnspan = 2)

# Checkbox to toggle between on and off state what triggers a function when it is either toggled on or off.
deltemp_checkbox = tk.Checkbutton(checkbox_frame, text = ' Automatically delete large temporary files', variable = deltemp_var, onvalue = 1, offvalue = 0)
CreateToolTip(deltemp_checkbox, text = 'Check this box if you will not need\nthe large-sized intermediary files')
deltemp_checkbox.grid(sticky = "w", row = 4, column = 1, padx = 10, columnspan = 2)

# Checkbox to toggle between on and off state what triggers a function when it is either toggled on or off.
stdout_checkbox = tk.Checkbutton(checkbox_frame, text = ' Record standard outputs & errors', variable = stdout_var, onvalue = 1, offvalue = 0)
CreateToolTip(stdout_checkbox, text = 'Check this box if you want to save\npipeline standard outputs (channel 1>) and\nstandard errors (channel 2>) as text files')
stdout_checkbox.grid(sticky = "w", row = 5, column = 1, padx = 10, columnspan = 2)

# Create navigation buttons
checkbox_back_button = tk.Button(checkbox_frame, text = "<< Back", command = lambda : change_frame_function(checkbox_frame, meme_motif_frame), width = 20)
checkbox_continue_button = tk.Button(checkbox_frame, text = "Continue >>", command = lambda : change_frame_function(checkbox_frame, execute_frame), width = 20)
checkbox_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
checkbox_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################


# Function to create a new sample table based on all the samples registered in the GUI variables, and save it in the output save directory,
#   serving as a reproducible documentation for the analysis results, and as a valid input usable for a new ChIP-AP pipeline run
def print_sample_table_function():
    sample_table_output_df = pd.DataFrame.from_dict(sample_table_output_dict, orient = 'index') 
    sample_table_output_df = sample_table_output_df.transpose()
    sample_table_output_df.to_csv('{}/{}_sample_table.tsv'.format(output_dir_arg.get(), setname_string_var.get()), sep = '\t', index = False)


# Function to create a new setting table based on all program flags and arguments registered in the GUI variables, and save it in the output save directory,
#   serving as a reproducible documentation for the analysis results, and as a valid input usable for a new ChIP-AP pipeline run
def print_setting_table_function():
    suite_program_arg = [
        fastqc1_arg.get(),
        clumpify_arg.get(),
        bbduk_arg.get(),
        trimmomatic_arg.get(),
        fastqc2_arg.get(),
        bwa_mem_arg.get(),
        samtools_view_arg.get(),
        plotfingerprint_arg.get(),
        fastqc3_arg.get(),
        macs2_callpeak_arg.get(),
        gem_arg.get(),
        sicer2_arg.get(),
        homer_findPeaks_arg.get(),
        genrich_arg.get(),
        homer_mergePeaks_arg.get(),
        homer_annotatePeaks_arg.get(),
        fold_change_calculator_arg.get(),
        homer_findMotifsGenome_arg.get(),
        meme_chip_arg.get()]

    setting_table_output_dict = {'program' : suite_program_list, 'argument' : suite_program_arg}
    setting_table_output_df = pd.DataFrame.from_dict(setting_table_output_dict)
    setting_table_output_df.to_csv('{}/{}_setting_table.tsv'.format(output_dir_arg.get(), setname_string_var.get()), sep = '\t', index = False)


# Function to create a text file containing the full command line which was used to call for ChIP-AP pipeline run, and save it in the output save directory,
#   serving as a reproducible documentation for the analysis results, and as a valid input usable for a new ChIP-AP pipeline run
def print_GS_command_line_function():
    with open('{}/{}_command_line.txt'.format(output_dir_arg.get(), setname_string_var.get()), 'w') as command_line_file:
        command_line_file.write(command_line_output_string_var.get() + '\n') # The command line version where user wants to only generate the pipeline scripts


# Function to create a text file containing the full command line which was used to call for ChIP-AP pipeline run, and save it in the output save directory,
#   serving as a reproducible documentation for the analysis results, and as a valid input usable for a new ChIP-AP pipeline run
def print_GSaR_command_line_function():
    with open('{}/{}_command_line.txt'.format(output_dir_arg.get(), setname_string_var.get()), 'w') as command_line_file:
        command_line_file.write(command_line_output_string_var.get() + ' --run\n') # The command line version where user wants to start the pipeline immediately


# Function to execute the command line version where user wants to only generate the pipeline scripts
def execute_GS_command_line_function():
    subprocess.run(command_line_output_string_var.get(), shell = True)


# Function to execute the command line version where user wants to start the pipeline immediately
def execute_GSaR_command_line_function():
    subprocess.run(command_line_output_string_var.get() + ' --run', shell = True)


# Function to execute a series of commands when user activates the "Generate scripts" button
def generate_scripts_function():
    if not os.path.exists(output_dir_arg.get()):
        os.makedirs(output_dir_arg.get())
    
    print_sample_table_function() # Create a new sample table based on all the samples registered in the GUI variables, and save it in the output save directory
    print_setting_table_function() # Create a new setting table based on all program flags and arguments registered in the GUI variables, and save it in the output save directory
    print_GS_command_line_function() # Create a text file containing the full command line which was used to call for ChIP-AP pipeline run, and save it in the output save directory
    execute_GS_command_line_function() # Generate the pipeline scripts
    execute_exit_button.config(state = tk.NORMAL) # Unlock the button to exit the GUI
    cpu_count_notification_string_var.set("Scripts generated! Check them out in:\n{}".format(output_dir_arg.get())) # Notify the user
    cpu_count_notification_label.config(fg = 'green')


# Function to execute a series of commands when user activates the "Generate scripts and run" button
def generate_scripts_and_run_function():
    if not os.path.exists(output_dir_arg.get()):
        os.makedirs(output_dir_arg.get())
    
    print_sample_table_function() # Create a new sample table based on all the samples registered in the GUI variables, and save it in the output save directory
    print_setting_table_function() # Create a new setting table based on all program flags and arguments registered in the GUI variables, and save it in the output save directory
    print_GSaR_command_line_function() # Create a text file containing the full command line which was used to call for ChIP-AP pipeline run, and save it in the output save directory
    execute_exit_button.config(state = tk.NORMAL) # Unlock the button to exit the GUI
    cpu_count_notification_string_var.set("Pipeline started! Check your results later in:\n{}".format(output_dir_arg.get())) # Notify the user
    cpu_count_notification_label.config(fg = 'green')
    execute_GSaR_command_line_function() # Start the pipeline immediately


execute_frame = tk.Frame(root) # Make a new frame

# Simple label with static text
cpu_count_label = tk.Label(execute_frame, text = "CPU cores to use:", width = 20)
cpu_count_label.grid(row = 1, column = 1, padx = (10,5))

cpu_count_entry = tk.Entry(execute_frame, width = 20, textvariable = cpu_count_string_var)
CreateToolTip(cpu_count_entry, text = 'Type in the number of CPU cores\nto be used by the pipeline')
cpu_count_entry.grid(row = 1, column = 2, padx = (5,10), ipady = 3)

cpu_count_notification_label = tk.Label(execute_frame, textvariable = cpu_count_notification_string_var, width = 50)
cpu_count_notification_label.grid(row = 2, column = 1, padx = 10, columnspan = 2)

# Clickable button that triggers a function when activated
generate_scripts_button = tk.Button(execute_frame, text = 'Generate scripts', command = lambda : generate_scripts_function(), width = 25)
CreateToolTip(generate_scripts_button, text = 'Select this if you wish to\nrun the pipeline later\nby executing MASTER_script.sh\nwithin the output save folder')
generate_and_run_scripts_button = tk.Button(execute_frame, text = 'Generate and run scripts', command = lambda : generate_scripts_and_run_function(), width = 25)
CreateToolTip(generate_and_run_scripts_button, text = 'Select this if you wish to\nrun the pipeline now\nNOTE: It may take up to\nseveral hours depending\non your system')

# Clickable button that triggers a function when activated
generate_scripts_button.grid(row = 20, column = 1, sticky = "w", padx = (10,5))
generate_and_run_scripts_button.grid(row = 20, column = 2, sticky = "e", padx = (5,10))

# Create navigation buttons
execute_back_button = tk.Button(execute_frame, text = "<< Back", command = lambda : change_frame_function(execute_frame, checkbox_frame), width = 20)
execute_exit_button = tk.Button(execute_frame, text = "Exit wizard", command = lambda : exit(), width = 20, state = tk.DISABLED)
execute_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
execute_exit_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################


# Initial run of all variable tracing functions at GUI start, just right after everything finished loading.
# Traces all the initial stored values at GUI start, so GUI text notifications can inform user of what the GUI needs from the very beginning.
register_sample_function()
update_command_line_function()
rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path

root.mainloop() # Initiates the GUI looping routine


########################################################################################################################
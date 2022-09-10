#!/usr/bin/env python3
#pyright: reportUnboundVariable=false


script_version = '2.0'


import tkinter as tk
import tkinter.font as tkFont
from tkinter import filedialog
from tkinter import ttk
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
default_setting_table_file_full_path = '{}/default_settings_table.tsv'.format(genome_folder_full_path) # Path to the default setting table

# Create a non-manually-resizable GUI base frame (root)
root = tk.Tk()
root.title(chipap_program_name)
root.resizable(width = False, height = False)

line_style = ttk.Style()
line_style.configure("Line.TSeparator", background = 'black')

# Make vertical separator lines
ttk.Separator(root, style = "Line.TSeparator", orient = tk.VERTICAL).grid(row = 0, column = 15, rowspan = 6, sticky = "ns")
ttk.Separator(root, style = "Line.TSeparator", orient = tk.VERTICAL).grid(row = 0, column = 10, rowspan = 11, sticky = "ns")
ttk.Separator(root, style = "Line.TSeparator", orient = tk.VERTICAL).grid(row = 30, column = 6, rowspan = 17, sticky = "ns")
ttk.Separator(root, style = "Line.TSeparator", orient = tk.VERTICAL).grid(row = 30, column = 12, rowspan = 17, sticky = "ns")


# Set the default fonts to be used in the GUI
default_font = tkFont.nametofont("TkDefaultFont")
default_font.configure(family = 'fixed', size = 16)

text_font = tkFont.nametofont("TkTextFont")
text_font.configure(family = 'fixed', size = 16)

fixed_font = tkFont.nametofont("TkFixedFont")
fixed_font.configure(family = 'fixed', size = 16)

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
genome_ref_options = ["hg38 (Homo sapiens)", 
                        "hg19 (Homo sapiens)", 
                        "mm9 (Mus musculus)", 
                        "mm10 (Mus musculus)", 
                        "dm6 (Drosophila melanogaster)", 
                        "sacCer3 (Saccharomyces cerevisiae)",
                        "Other [!!!under construction!!!]"]

# List of supported choices of peak list to analyze for motif enrichment by HOMER(options for homer_motif_drop_down menu)
homer_motif_options = ["None",
                        "Consensus peak set", 
                        "Union peak set", 
                        "Both peak sets"]

# List of supported choices of peak list to analyze for motif enrichment by MEME (options for meme_motif_drop_down menu)
meme_motif_options = ["None",
                        "Consensus peak set", 
                        "Union peak set", 
                        "Both peak sets"]


########################################################################################################################


# Define tooltip object properties for the GUI
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
        x = x + self.widget.winfo_rootx() + 100
        y = y + cy + self.widget.winfo_rooty() + 35
        self.tipwindow = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        tw.wm_geometry("+%d+%d" % (x, y))
        label = tk.Label(tw, text = self.text, justify = tk.LEFT,
                      background = "#ffffe0", relief = tk.SOLID, borderwidth = 1,
                      font = ("tahoma", "12", "normal"))
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

    error_type_A_int_var.set(0) # Reset the error type A (samples error) detector variable

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
    check_required_input_function() # Check if the required inputs (type B) other than samples (type A) are satisfied or contain any errors
    update_command_line_function() # Update the generated ChIP-AP command line based on the latest registered variable values


# First step of ChIP samples error checking: Check to see if they are at least assigned
def check_chip_sample_assigned_function():

    if read_mode_string_var.get() == 'single' or read_mode_string_var.get() == 'paired': # Perform checks if dataset sequencing mode is single-end or paired-end
        
        if read_mode_string_var.get() == 'single': # In case of single-end dataset
            if all(chip_r1 == '' for chip_r1 in chip_list_r1): # If r1 list is all empty strings
                chip_sample_notification_string_var.set('Please assign ChIP sample') # Notify the user
                chip_sample_notification_label.config(fg = 'blue')
                error_type_A_int_var.set(1) # Set error type A (sample error) to 1
                return # Exit the function

        if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
            if all(chip_r1 == '' for chip_r1 in chip_list_r1): # If r1 list is all empty strings
                chip_sample_notification_string_var.set('Please assign ChIP sample (read 1)') # Notify the user
                chip_sample_notification_label.config(fg = 'blue')
                error_type_A_int_var.set(1) # Set error type A (sample error) to 1
                return # Exit the function

            # Check for r2 only if the input files are not aligned
            if all(chip_r2 == '' for chip_r2 in chip_list_r2) and not all('.bam' in sample for sample in (chip_list_r1 + chip_list_r1) if sample != ''): # If r2 list is all empty strings
                chip_sample_notification_string_var.set('Please assign ChIP sample (read 2)') # Notify the user
                chip_sample_notification_label.config(fg = 'blue')
                error_type_A_int_var.set(1) # Set error type A (sample error) to 1
                return # Exit the function

        # If the function reaches this point, it means there are at least one sample assigned to each of the lists necessary

        if read_mode_string_var.get() == 'single': # In case of single-end dataset
            if not bool(chip_rep1_r1_string_var.get()): # If the assigned sample(s) is not the first replicate
                chip_sample_notification_string_var.set('ChIP samples are not assigned from the top (replicate 1)') # Notify the user
                chip_sample_notification_label.config(fg = 'red') # Highlight where the problem is
                chip_rep1_r1_entry.config(bg = 'IndianRed1')
                error_type_A_int_var.set(1) # Set error type A (sample error) to 1
                return # Exit the function
        
        if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
            if not bool(chip_rep1_r1_string_var.get()) and not bool(chip_rep1_r2_string_var.get()): # If the assigned sample(s) is not the first replicate of r1 nor r2
                chip_sample_notification_string_var.set('ChIP samples are not assigned from the top (replicate 1)') # Notify the user
                chip_sample_notification_label.config(fg = 'red')
                chip_rep1_r1_entry.config(bg = 'IndianRed1') # Highlight where the problem is
                chip_rep1_r2_entry.config(bg = 'IndianRed1') # Highlight where the problem is
                error_type_A_int_var.set(1) # Set error type A (sample error) to 1
                return # Exit the function

            if not bool(chip_rep1_r1_string_var.get()): # If the assigned sample(s) is not the first replicate of r1
                chip_sample_notification_string_var.set('ChIP samples are not assigned from the top (replicate 1)') # Notify the user
                chip_sample_notification_label.config(fg = 'red')
                chip_rep1_r1_entry.config(bg = 'IndianRed1') # Highlight where the problem is
                error_type_A_int_var.set(1) # Set error type A (sample error) to 1
                return # Exit the function

            if not bool(chip_rep1_r2_string_var.get()): # If the assigned sample(s) is not the first replicate of r2
                chip_sample_notification_string_var.set('ChIP samples are not assigned from the top (replicate 1)') # Notify the user
                chip_sample_notification_label.config(fg = 'red')
                chip_rep1_r2_entry.config(bg = 'IndianRed1') # Highlight where the problem is
                error_type_A_int_var.set(1) # Set error type A (sample error) to 1
                return # Exit the function

        # If no errors found and the function did not exit prematurely because of all the tests above:
        chip_sample_notification_string_var.set('ChIP samples have been assigned') # Notify the user
        chip_sample_notification_label.config(fg = 'green')
    
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
                error_type_A_int_var.set(1) # Set error type A (sample error) to 1
                return

        if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
            if all(ctrl_r1 == '' for ctrl_r1 in ctrl_list_r1): # If r1 list is all empty strings
                ctrl_sample_notification_string_var.set('Please assign control sample (read 1)') # Notify the user
                ctrl_sample_notification_label.config(fg = 'blue')
                error_type_A_int_var.set(1) # Set error type A (sample error) to 1
                return

            # Check for r2 only if the input files are not aligned
            if all(ctrl_r2 == '' for ctrl_r2 in ctrl_list_r2) and not all('.bam' in sample for sample in (ctrl_list_r1 + ctrl_list_r1) if sample != ''): # If r2 list is all empty strings
                ctrl_sample_notification_string_var.set('Please assign control sample (read 2)') # Notify the user
                ctrl_sample_notification_label.config(fg = 'blue')
                error_type_A_int_var.set(1) # Set error type A (sample error) to 1
                return

        # If the function reaches this point, it means there are at least one sample assigned to each of the lists necessary

        if read_mode_string_var.get() == 'single': # In case of single-end dataset
            if not bool(ctrl_rep1_r1_string_var.get()): # If the assigned sample(s) is not the first replicate
                ctrl_sample_notification_string_var.set('Control samples are not assigned from the top (replicate 1)') # Notify the user
                ctrl_sample_notification_label.config(fg = 'red')
                ctrl_rep1_r1_entry.config(bg = 'IndianRed1') # Highlight where the problem is
                error_type_A_int_var.set(1) # Set error type A (sample error) to 1
                return
        
        if read_mode_string_var.get() == 'paired': # In case of paired-end dataset
            if not bool(ctrl_rep1_r1_string_var.get()) and not bool(ctrl_rep1_r2_string_var.get()): # If the assigned sample(s) is not the first replicate of r1 nor r2
                ctrl_sample_notification_string_var.set('Control samples are not assigned from the top (replicate 1)') # Notify the user
                ctrl_sample_notification_label.config(fg = 'red')
                ctrl_rep1_r1_entry.config(bg = 'IndianRed1') # Highlight where the problem is
                ctrl_rep1_r2_entry.config(bg = 'IndianRed1') # Highlight where the problem is
                error_type_A_int_var.set(1) # Set error type A (sample error) to 1
                return
            
            elif not bool(ctrl_rep1_r1_string_var.get()): # If the assigned sample(s) is not the first replicate of r1
                ctrl_sample_notification_string_var.set('Control samples are not assigned from the top (replicate 1)') # Notify the user
                ctrl_sample_notification_label.config(fg = 'red')
                ctrl_rep1_r1_entry.config(bg = 'IndianRed1') # Highlight where the problem is
                error_type_A_int_var.set(1) # Set error type A (sample error) to 1
                return

            elif not bool(ctrl_rep1_r2_string_var.get()): # If the assigned sample(s) is not the first replicate of r2
                ctrl_sample_notification_string_var.set('Control samples are not assigned from the top (replicate 1)') # Notify the user
                ctrl_sample_notification_label.config(fg = 'red')
                ctrl_rep1_r2_entry.config(bg = 'IndianRed1') # Highlight where the problem is
                error_type_A_int_var.set(1) # Set error type A (sample error) to 1
                return

        # If no errors found and the function did not exit prematurely because of all the tests above:
        ctrl_sample_notification_string_var.set('Control samples have been assigned') # Notify the user
        ctrl_sample_notification_label.config(fg = 'green')

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
            error_type_A_int_var.set(1) # Set error type A (sample error) to 1
            return # Exit the function

        elif chip_pair_error_state == 0: # If no sample pairing error is detected above 
            chip_sample_notification_string_var.set('All ChIP samples are properly paired') # Notify the user
            chip_sample_notification_label.config(fg = 'green')
        
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
            error_type_A_int_var.set(1) # Set error type A (sample error) to 1
            return # Exit the function

        elif ctrl_pair_error_state == 0: # If no sample pairing error is detected above 
            ctrl_sample_notification_string_var.set('All control samples are properly paired') # Notify the user
            ctrl_sample_notification_label.config(fg = 'green')

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
        error_type_A_int_var.set(1) # Set error type A (sample error) to 1
        return # Exit the function

    elif chip_sample_format_error_state == 0: # If no sample file extension error is detected above
        chip_sample_notification_string_var.set('All ChIP samples have valid file extension') # Notify the user
        chip_sample_notification_label.config(fg = 'green')

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
        error_type_A_int_var.set(1) # Set error type A (sample error) to 1
        return # Exit the function

    elif ctrl_sample_format_error_state == 0: # If no sample file extension error is detected above
        ctrl_sample_notification_string_var.set('All control samples have valid file extension') # Notify the user
        ctrl_sample_notification_label.config(fg = 'green')

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
        error_type_A_int_var.set(1) # Set error type A (sample error) to 1

    elif chip_file_exist_error_state == 0: # If no sample file not found error is detected above
        chip_sample_notification_string_var.set('No problem found in ChIP samples') # Notify the user
        chip_sample_notification_label.config(fg = 'green')


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
        error_type_A_int_var.set(1) # Set error type A (sample error) to 1

    elif ctrl_file_exist_error_state == 0: # If no sample file not found error is detected above
        ctrl_sample_notification_string_var.set('No problem found in control samples') # Notify the user
        ctrl_sample_notification_label.config(fg = 'green')


 # Function to check if the required inputs (type B) other than samples (type A) are satisfied or contain any errors
def check_required_input_function(*args):
    
    error_type_B_int_var.set(0) # Reset the error type B (non-sample error) detector variable

    if read_mode_string_var.get() != 'single' and read_mode_string_var.get() != 'paired': # If inputted sequencing mode is none of the accepted choices
        error_type_B_int_var.set(1) # Set error type B (non-sample error) to 1

    if peak_type_string_var.get() != 'narrow' and peak_type_string_var.get() != 'broad' and peak_type_string_var.get() != 'unsure': # If inputted peak type is none of the accepted choices 
        error_type_B_int_var.set(1) # Set error type B (non-sample error) to 1

    if bool(sample_table_string_var.get()): # Perform the check only if the sample table variable value is not empty
        if not os.path.isfile(sample_table_string_var.get()): # If the entered path leads to a non-existent file
            error_type_B_int_var.set(1) # Set error type B (non-sample error) to 1

    if bool(setting_table_string_var.get()): # Perform the check only if the setting table variable value is not empty
        if '[MODIFIED]' in setting_table_string_var.get(): # If the setting table variable value contains an indicator of setting values manual modification
            pass # Do nothing
        elif not os.path.isfile(setting_table_string_var.get()): # If the entered path leads to a non-existent file
            error_type_B_int_var.set(1) # Set error type B (non-sample error) to 1

    if not bool(genome_ref_string_var.get()): # If the reference genome build variable value is empty
        error_type_B_int_var.set(1) # Set error type B (non-sample error) to 1
    
    if not bool(genome_folder_string_var.get()): # If the reference genome directory variable value is empty
        error_type_B_int_var.set(1) # Set error type B (non-sample error) to 1
    elif not os.path.isdir(genome_folder_string_var.get()): # If the entered path leads to a non-existent directory
        error_type_B_int_var.set(1) # Set error type B (non-sample error) to 1
        genome_folder_entry.config(fg = 'red') # Highlight where the problem is
    elif os.path.isdir(genome_folder_string_var.get()): # If the entered path leads to an existing directory
        genome_folder_entry.config(fg = 'black') # Remove the highlight, if previously highlighted

    if bool(known_motif_string_var.get()): # Perform the check only if the known motif variable value is not empty
        if not os.path.isfile(known_motif_string_var.get()): # If the entered path leads to a non-existent file
            error_type_B_int_var.set(1) # Set error type B (non-sample error) to 1
            known_motif_entry.config(fg = 'red') # Highlight where the problem is
        elif os.path.isfile(known_motif_string_var.get()): # If the entered path leads to an existing file
            known_motif_entry.config(fg = 'black') # Remove the highlight, if previously highlighted

    if not bool(setname_string_var.get()): # If the dataset name variable value is empty
        error_type_B_int_var.set(1) # Set error type B (non-sample error) to 1
    
    if not bool(output_folder_string_var.get()): # If the output directory variable value is empty
        error_type_B_int_var.set(1) # Set error type B (non-sample error) to 1

    if not bool(cpu_count_string_var.get()): # If the CPU count variable value is empty
        error_type_B_int_var.set(1) # Set error type B (non-sample error) to 1
        cpu_count_entry.config(fg = 'black') # Set font color for CPU count entry field prior to any input
        cpu_count_notification_string_var.set("Maximum number of CPU cores available: {}".format(max_cpu)) # Inform the user of the acceptable range of values
        cpu_count_notification_label.config(fg = 'blue')
    
    elif bool(cpu_count_string_var.get()): # If the CPU count variable value is not empty
        try: # If the entered value can be converted into python integer
            input_cpu_count = int(cpu_count_string_var.get())

            if input_cpu_count > max_cpu: # If the CPU count variable value is larger than maximum capacity
                error_type_B_int_var.set(1) # Set error type B (non-sample error) to 1
                cpu_count_entry.config(fg = 'red') # Highlight where the problem is
                cpu_count_notification_string_var.set("Entered number exceeds available CPU cores ({})".format(max_cpu)) # Notify the user
                cpu_count_notification_label.config(fg = 'red')

            elif input_cpu_count < 1: # If the CPU count variable value is lower than minimum needed
                error_type_B_int_var.set(1) # Set error type B (non-sample error) to 1
                cpu_count_entry.config(fg = 'red') # Highlight where the problem is
                cpu_count_notification_string_var.set("Need at least one CPU core to run the pipeline") # Notify the user
                cpu_count_notification_label.config(fg = 'red')

            elif input_cpu_count <= max_cpu: # If the CPU count variable value is within acceptable range of values
                cpu_count_entry.config(fg = 'black')
                cpu_count_notification_string_var.set("") # Clear the notification
                cpu_count_notification_label.config(fg = 'blue')

        except: # If the entered value cannot be converted into python integer
            error_type_B_int_var.set(1) # Set error type B (non-sample error) to 1
            cpu_count_entry.config(fg = 'red') # Highlight where the problem is
            cpu_count_notification_string_var.set("Entered value is not an integer") # Notify the user
            cpu_count_notification_label.config(fg = 'red')

    generate_button_switch_function() # Lock or unlock the "Generate scripts" and "Generate scripts and run" buttons based on all registered GUI variables


# Function to control the locking and unlocking of the "Generate scripts" and "Generate scripts and run" buttons based on all registered GUI variables
def generate_button_switch_function():
    if error_type_A_int_var.get() == 0 and error_type_B_int_var.get() == 0: # If neither of error type A nor B was set to 1 above (no error detected)
        generate_scripts_button.config(state = tk.NORMAL) # Enable the big red button
        generate_and_run_scripts_button.config(state = tk.NORMAL) # Enable the big red button
        cpu_count_notification_string_var.set('ChIP-AP ready to go!') # Notify the user
        cpu_count_notification_label.config(fg = 'green')
    if error_type_A_int_var.get() == 1 or error_type_B_int_var.get() == 1: # If either of error type A or B was set to 1 above (error detected)
        generate_scripts_button.config(state = tk.DISABLED) # Disable the big red button
        generate_and_run_scripts_button.config(state = tk.DISABLED) # Disable the big red button


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

    if homer_motif_string_var.get() != 'None': # If user has decided which peak set to perform HOMER motif enrichment analysis on
        # Assign the HOMER motif enrichment analysis argument (consensus, union, or both) for the command line behind the corresponding flag
        homer_motif_arg.set(' --homer_motif {}'.format(homer_motif_string_var.get().split(' ')[0].lower()))
    elif homer_motif_string_var.get() == 'None': # If user has not decided which peak set to perform HOMER motif enrichment analysis on
        homer_motif_arg.set('') # Empty the HOMER motif enrichment analysis argument for the command line

    if meme_motif_string_var.get() != 'None': # If user has decided which peak set to perform MEME motif enrichment analysis on
        # Assign the MEME motif enrichment analysis argument (consensus, union, or both) for the command line behind the corresponding flag
        meme_motif_arg.set(' --meme_motif {}'.format(meme_motif_string_var.get().split(' ')[0].lower()))
    elif meme_motif_string_var.get() == 'None': # If user has not decided which peak set to perform MEME motif enrichment analysis on
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
    check_required_input_function() # Check if the required inputs other than sample errors are satisifed or contain any errors
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


error_type_A_int_var = tk.IntVar(value = 0) # To store the binary status of sample type errors
error_type_B_int_var = tk.IntVar(value = 0) # To store the binary status of non-sample type errors

max_cpu = multiprocessing.cpu_count() # Read available CPU cores
current_dir_string_var = tk.StringVar(value = default_current_dir) # Path the the last browsed folder

read_mode_string_var = tk.StringVar() # Argument for the --mode flag
peak_type_string_var = tk.StringVar() # Argument for the --peak flag
sample_table_string_var = tk.StringVar() # Path to setting table file which values are to be loaded
sample_table_notification_string_var = tk.StringVar() # GUI text notification for the sample table loading status
chip_sample_notification_string_var = tk.StringVar() # GUI text notification for the ChIP samples check status
ctrl_sample_notification_string_var = tk.StringVar() # GUI text notification for the control samples check status
setting_table_string_var = tk.StringVar(value = default_setting_table_file_full_path) # Path to setting table file which values are to be loaded
setting_table_notification_string_var = tk.StringVar(value = "Currently using default settings table") # GUI text notification for the setting table loading status
genome_ref_string_var = tk.StringVar(value = genome_ref_options[0]) # Argument for the --ref flag
genome_folder_string_var = tk.StringVar(value = genome_folder_full_path) # Argument for --genome flag
known_motif_string_var = tk.StringVar() # Argument for the --motif flag
setname_string_var = tk.StringVar() # Argument for the --setname flag
output_folder_string_var = tk.StringVar() # Argument for the --output flag
fcmerge_var = tk.IntVar() # Binary value for the --fcmerge flag
goann_var = tk.IntVar() # Binary value for the --goann flag
pathann_var = tk.IntVar() # Binary value for the --pathann flag
deltemp_var = tk.IntVar(value = 1) # Binary value for the --deltemp flag.
stdout_var = tk.IntVar() # Binary value for the 1> and 2> channel flags
homer_motif_string_var = tk.StringVar(value = homer_motif_options[0]) # Argument for --homer_motif flag
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


# Insert full ChIP-AP logo
image = Image.open(chipap_logo_full_path)
image = image.resize((300, 110), Image.ANTIALIAS)
photo = ImageTk.PhotoImage(image)

# Place the logo on the top left of the dashboard
chipap_logo = tk.Label(root, image = photo, anchor = "center")
chipap_logo.image = photo
chipap_logo.grid(row = 1, column = 1, rowspan = 4, columnspan = 3, sticky = "w",padx = 3)

# Write few sentences about ChIP-AP beside the logo
chipap_about = tk.Text(root, width = 61, height = 7, relief = tk.FLAT, bg = 'gray85',font = (None, 9))
chipap_about_text = 'Integrated Analysis Pipeline for Unbiased ChIP-seq Analysis\n\nComplete guides and walkthroughs can be found on our github\n(https://github.com/JSuryatenggara/ChIP-AP)\n\nIf you use ChIP-AP please cite our bioRxiv pre-print article\n(https://www.biorxiv.org/content/10.1101/2021.04.18.440382v1)'
chipap_about.insert(tk.END, chipap_about_text)
chipap_about.config(state = tk.DISABLED)
chipap_about.tag_configure("center", justify = tk.CENTER)
chipap_about.tag_add("center", 1.0, tk.END)
chipap_about.grid(row = 1, column = 3, rowspan = 4, columnspan = 7, sticky = "e", padx = 5, pady = 2)


########################################################################################################################


# This function is triggered when user chooses the dataset sequencing mode
def sample_button_switch_function():

    if read_mode_string_var.get() == 'single': # Sample R2 widgets are all disabled when single-end mode
        r1_state = tk.NORMAL
        r2_state = tk.DISABLED
        sample_table_button_state = tk.NORMAL

    elif read_mode_string_var.get() == 'paired':
        r1_state = tk.NORMAL
        r2_state = tk.NORMAL
        sample_table_button_state = tk.NORMAL

    else: # All sample related widgets are all disabled when user has not chosen any dataset sequencing mode
        r1_state = tk.DISABLED
        r2_state = tk.DISABLED
        sample_table_button_state = tk.DISABLED

    chip_rep1_r1_button.config(state = r1_state)
    chip_rep1_r2_button.config(state = r2_state)
    chip_rep2_r1_button.config(state = r1_state)
    chip_rep2_r2_button.config(state = r2_state)
    chip_rep3_r1_button.config(state = r1_state)
    chip_rep3_r2_button.config(state = r2_state)
    chip_rep4_r1_button.config(state = r1_state)
    chip_rep4_r2_button.config(state = r2_state)
    chip_rep5_r1_button.config(state = r1_state)
    chip_rep5_r2_button.config(state = r2_state)

    ctrl_rep1_r1_button.config(state = r1_state)
    ctrl_rep1_r2_button.config(state = r2_state)
    ctrl_rep2_r1_button.config(state = r1_state)
    ctrl_rep2_r2_button.config(state = r2_state)
    ctrl_rep3_r1_button.config(state = r1_state)
    ctrl_rep3_r2_button.config(state = r2_state)
    ctrl_rep4_r1_button.config(state = r1_state)
    ctrl_rep4_r2_button.config(state = r2_state)
    ctrl_rep5_r1_button.config(state = r1_state)
    ctrl_rep5_r2_button.config(state = r2_state)

    chip_rep1_r1_entry.config(state = r1_state)
    chip_rep1_r2_entry.config(state = r2_state)
    chip_rep2_r1_entry.config(state = r1_state)
    chip_rep2_r2_entry.config(state = r2_state)
    chip_rep3_r1_entry.config(state = r1_state)
    chip_rep3_r2_entry.config(state = r2_state)
    chip_rep4_r1_entry.config(state = r1_state)
    chip_rep4_r2_entry.config(state = r2_state)
    chip_rep5_r1_entry.config(state = r1_state)
    chip_rep5_r2_entry.config(state = r2_state)

    ctrl_rep1_r1_entry.config(state = r1_state)
    ctrl_rep1_r2_entry.config(state = r2_state)
    ctrl_rep2_r1_entry.config(state = r1_state)
    ctrl_rep2_r2_entry.config(state = r2_state)
    ctrl_rep3_r1_entry.config(state = r1_state)
    ctrl_rep3_r2_entry.config(state = r2_state)
    ctrl_rep4_r1_entry.config(state = r1_state)
    ctrl_rep4_r2_entry.config(state = r2_state)
    ctrl_rep5_r1_entry.config(state = r1_state)
    ctrl_rep5_r2_entry.config(state = r2_state)

    sample_table_button.config(state = sample_table_button_state)
    clear_sample_button.config(state = sample_table_button_state)
    sample_table_entry.config(state = sample_table_button_state)

    check_required_input_function() # Check if the required inputs (type B) other than samples (type A) are satisfied or contain any errors
    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path

    if read_mode_string_var.get(): # Prompt the user after they have chosen a dataset sequencing mode
        sample_table_notification_string_var.set('Please load sample table or assign samples manually below')
        sample_table_notification_label.config(fg = 'blue')

# Simple label with static text
read_mode_label = tk.Label(root, text = "Dataset sequencing mode:", justify = tk.LEFT, width = 30)
read_mode_label.grid(row = 1, column = 11, columnspan = 4, padx = 10, pady = (5,2))

# Radio button of the first choice that triggers a function when activated
single_end_radio = tk.Radiobutton(root, text = "Single end", padx = 10, variable = read_mode_string_var, value = 'single', width = 30, command = lambda : sample_button_switch_function())
CreateToolTip(single_end_radio, text = 'Select this if there is one\noutput file per sample')
single_end_radio.grid(row = 2, column = 11, columnspan = 4, padx = 10)

# Radio button of the second choice that triggers a function when activated
paired_end_radio = tk.Radiobutton(root, text = "Paired end", padx = 10, variable = read_mode_string_var, value = 'paired', width = 30, command = lambda : sample_button_switch_function())
CreateToolTip(paired_end_radio, text = 'Select this if there are two\noutput files per sample.\nThey are typically in pairs:\nR1 and R2 for every sample')
paired_end_radio.grid(row = 3, column = 11, columnspan = 4, padx = 10, pady = (0,5))


########################################################################################################################


# This function is triggered when user chooses the dataset peak type
def mea_switch_function(): # Disable the motif enrichment analysis (mea)-related widgets if dataset peak type is broad 
    if peak_type_string_var.get() == 'narrow':
        homer_motif_label.config(fg = 'black')
        meme_motif_label.config(fg = 'black')
        homer_motif_drop_down.config(state = tk.NORMAL)
        meme_motif_drop_down.config(state = tk.NORMAL)
    elif peak_type_string_var.get() == 'broad':
        homer_motif_string_var.set(homer_motif_options[0])
        meme_motif_string_var.set(meme_motif_options[0])
        homer_motif_label.config(fg = 'grey60')
        meme_motif_label.config(fg = 'grey60')
        homer_motif_drop_down.config(state = tk.DISABLED)
        meme_motif_drop_down.config(state = tk.DISABLED)
    elif peak_type_string_var.get() == 'unsure':
        homer_motif_label.config(fg = 'black')
        meme_motif_label.config(fg = 'black')
        homer_motif_drop_down.config(state = tk.NORMAL)
        meme_motif_drop_down.config(state = tk.NORMAL)

    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path

# Simple label with static text
peak_type_label = tk.Label(root, text = "Dataset peak type:", justify = tk.LEFT, width = 30)
peak_type_label.grid(row = 1, column = 16, columnspan = 4, padx = 10, pady = (5,2))

# Radio button of the first choice that triggers a function when activated
narrow_peak_radio = tk.Radiobutton(root, text = "Narrow peaks", padx = 10, variable = peak_type_string_var, value = 'narrow', width = 30, command = lambda : mea_switch_function())
CreateToolTip(narrow_peak_radio, text = 'Select this for ChIP-seq experiment\nusing transcription factor protein')
narrow_peak_radio.grid(row = 2, column = 16, columnspan = 4, padx = 10)

# Radio button of the second choice that triggers a function when activated
broad_peak_radio = tk.Radiobutton(root, text = "Broad peaks", padx = 10, variable = peak_type_string_var, value = 'broad', width = 30, command = lambda : mea_switch_function())
CreateToolTip(broad_peak_radio, text = 'Select this for ChIP-seq experiment\nusing chromatin modifier protein')
broad_peak_radio.grid(row = 3, column = 16, columnspan = 4, padx = 10)

# Radio button of the third choice that triggers a function when activated
unsure_peak_radio = tk.Radiobutton(root, text = "Unsure", padx = 10, variable = peak_type_string_var, value = 'unsure', width = 30, command = lambda : mea_switch_function())
CreateToolTip(unsure_peak_radio, text = 'Select this for ChIP-seq experiment\nusing protein with both possibilities')
unsure_peak_radio.grid(row = 4, column = 16, columnspan = 4, padx = 10, pady = (0,5))

# Horizontal separator line
ttk.Separator(root, style = "Line.TSeparator", orient = tk.HORIZONTAL).grid(row = 5, column = 0, columnspan = 21, sticky = "we")


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

# Dynamic label with changeable text
sample_table_notification_label = tk.Label(root, textvariable = sample_table_notification_string_var, width = 70, padx = 10, pady = 4, fg = 'blue')
sample_table_notification_string_var.set('Please specify data sequencing mode first!')
sample_table_notification_label.grid(row = 6, column = 1, padx = 10, pady = 4, columnspan = 9)

# Clickable button that triggers a function when activated
sample_table_button = tk.Button(root, text = 'Load sample table', command = lambda:sample_table_button_function(), state = tk.DISABLED, width = 15)
CreateToolTip(sample_table_button, text = 'Click here to browse and\nselect your sample table file')
sample_table_button.grid(sticky = "we", row = 7, column = 1, padx = 10, rowspan = 1, columnspan = 1)

# Entry field where user can manually key in their input to be stored as a variable
sample_table_entry = tk.Entry(root, textvariable = sample_table_string_var, width = 50, justify = tk.RIGHT, state = tk.DISABLED)
sample_table_entry.grid(sticky = "we", row = 7, column = 2, padx = (0,10), columnspan = 8, ipady = 3)

# Clickable button that triggers a function when activated
clear_sample_button = tk.Button(root, text = 'Clear samples', command = lambda:clear_sample_function('withfile'), state = tk.DISABLED, width = 15)
CreateToolTip(clear_sample_button, text = 'Click here to clear all\nassigned samples below')
clear_sample_button.grid(sticky = "we", row = 8, column = 1, padx = 10, pady = (4,5), rowspan = 1, columnspan = 1)


########################################################################################################################


# Function triggered by the "Load setting table" button, to open a browser window for choosing a local file
def setting_table_button_function():
    setting_table_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), mode = 'rb', title = 'Choose a setting table file')
    if setting_table_button_input: # If the user choose a file from the browser
        setting_table_string_var.set(setting_table_button_input.name) # Set the setting table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(setting_table_button_input.name)) # Set the last visited directory to the user-chosen file directory

    else: # If the user closed the browser without choosing any file from the browser
        clear_setting_function('withfile') # Clear all setting-related variables
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
            default_setting_button_function() # Restore default setting values
            setting_table_notification_string_var.set('Columns error, reverted to default') # Notify the user
            setting_table_notification_label.config(fg = 'red')
            rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path
            return False

        # Check first if they are both strings to avoid TypeError, then check if the table headers are 'program' and 'argument'
        if isinstance(setting_table_header[0], str) and isinstance(setting_table_header[1], str): # If both header values are string type
            if setting_table_header[0].strip().lower() != 'program' or setting_table_header[1].strip().lower() != 'argument': # If the header values are not 'program' and 'argument'
                clear_setting_function('withfile') # Clear all setting-related variables
                default_setting_button_function() # Restore default setting values
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
            default_setting_button_function() # Restore default setting values
            setting_table_notification_string_var.set('Header error, reverted to default') # Notify the user
            setting_table_notification_label.config(fg = 'red')
            rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path
            return False

    except: # If setting table loading fails
        clear_setting_function('withfile') # Clear all setting-related variables
        default_setting_button_function() # Restore default setting values
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


 # Function to restore default setting values
def default_setting_button_function():
    setting_table_string_var.set(default_setting_table_file_full_path) # Set the setting table GUI variable with the path to the default settings table file
    read_setting_table_function(default_setting_table_file_full_path) # Read the default settings table and load the contained values to their respective GUI variables
    setting_table_notification_string_var.set('Default values restored') # Notify the user
    setting_table_notification_label.config(fg = 'blue')

    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path


# Function to display the customize settings frame on top of the main frame
# If setting table loading fails (unless setting values modification present),
#   function is exited prematurely and GUI does not proceed to the customize settings frame
def setting_table_window_open():
    if '[MODIFIED]' in setting_table_string_var.get():
        pass
    else:
        if setting_table_loading_test_function() == False:
            return
        elif setting_table_loading_test_function() == True:
            pass

    global setting_table_window # Declare this customize settings frame as global so can be modified outside this function
    
    # Create a new frame on top of the main frame and display it
    setting_table_window = tk.Toplevel()
    setting_table_window.grab_set()

    # Simple labels with static texts
    fastqc1_arg_label = tk.Label(setting_table_window, text = 'fastqc1', width = 22, padx = 10, pady = 4, anchor = "e")
    clumpify_arg_label = tk.Label(setting_table_window, text = 'clumpify', width = 22, padx = 10, pady = 4, anchor = "e")
    bbduk_arg_label = tk.Label(setting_table_window, text = 'bbduk', width = 22, padx = 10, pady = 4, anchor = "e")
    trimmomatic_arg_label = tk.Label(setting_table_window, text = 'trimmomatic', width = 22, padx = 10, pady = 4, anchor = "e")
    fastqc2_arg_label = tk.Label(setting_table_window, text = 'fastqc2', width = 22, padx = 10, pady = 4, anchor = "e")
    bwa_mem_arg_label = tk.Label(setting_table_window, text = 'bwa_mem', width = 22, padx = 10, pady = 4, anchor = "e")
    samtools_view_arg_label = tk.Label(setting_table_window, text = 'samtools_view', width = 22, padx = 10, pady = 4, anchor = "e")
    plotfingerprint_arg_label = tk.Label(setting_table_window, text = 'plotfingerprint', width = 22, padx = 10, pady = 4, anchor = "e")
    fastqc3_arg_label = tk.Label(setting_table_window, text = 'fastqc3', width = 22, padx = 10, pady = 4, anchor = "e")
    macs2_callpeak_arg_label = tk.Label(setting_table_window, text = 'macs2_callpeak', width = 22, padx = 10, pady = 4, anchor = "e")
    gem_arg_label = tk.Label(setting_table_window, text = 'gem', width = 22, padx = 10, pady = 4, anchor = "e")
    sicer2_arg_label = tk.Label(setting_table_window, text = 'sicer2', width = 22, padx = 10, pady = 4, anchor = "e")
    homer_findPeaks_arg_label = tk.Label(setting_table_window, text = 'homer_findPeaks', width = 22, padx = 10, pady = 4, anchor = "e")
    genrich_arg_label = tk.Label(setting_table_window, text = 'genrich', width = 22, padx = 10, pady = 4, anchor = "e")
    homer_mergePeaks_arg_label = tk.Label(setting_table_window, text = 'homer_mergePeaks', width = 22, padx = 10, pady = 4, anchor = "e")
    homer_annotatePeaks_arg_label = tk.Label(setting_table_window, text = 'homer_annotatePeaks', width = 22, padx = 10, pady = 4, anchor = "e")
    fold_change_calculator_arg_label = tk.Label(setting_table_window, text = 'fold_change_calculator', width = 22, padx = 10, pady = 4, anchor = "e")
    homer_findMotifsGenome_label = tk.Label(setting_table_window, text = 'homer_findMotifsGenome', width = 22, padx = 10, pady = 4, anchor = "e")
    meme_chip_label = tk.Label(setting_table_window, text = 'meme_chip', width = 22, padx = 10, pady = 4, anchor = "e")

    # Display all the labels created above
    fastqc1_arg_label.grid(sticky = "we", row = 35, column = 11, padx = 10, columnspan = 1, pady = (10,0))
    clumpify_arg_label.grid(sticky = "we", row = 36, column = 11, padx = 10, columnspan = 1)
    bbduk_arg_label.grid(sticky = "we", row = 37, column = 11, padx = 10, columnspan = 1)
    trimmomatic_arg_label.grid(sticky = "we", row = 38, column = 11, padx = 10, columnspan = 1)
    fastqc2_arg_label.grid(sticky = "we", row = 39, column = 11, padx = 10, columnspan = 1)
    bwa_mem_arg_label.grid(sticky = "we", row = 40, column = 11, padx = 10, columnspan = 1)
    samtools_view_arg_label.grid(sticky = "we", row = 41, column = 11, padx = 10, columnspan = 1)
    plotfingerprint_arg_label.grid(sticky = "we", row = 42, column = 11, padx = 10, columnspan = 1)
    fastqc3_arg_label.grid(sticky = "we", row = 43, column = 11, padx = 10, columnspan = 1)
    macs2_callpeak_arg_label.grid(sticky = "we", row = 44, column = 11, padx = 10, columnspan = 1)
    gem_arg_label.grid(sticky = "we", row = 45, column = 11, padx = 10, columnspan = 1)
    sicer2_arg_label.grid(sticky = "we", row = 46, column = 11, padx = 10, columnspan = 1)
    homer_findPeaks_arg_label.grid(sticky = "we", row = 47, column = 11, padx = 10, columnspan = 1)
    genrich_arg_label.grid(sticky = "we", row = 48, column = 11, padx = 10, columnspan = 1)
    homer_mergePeaks_arg_label.grid(sticky = "we", row = 49, column = 11, padx = 10, columnspan = 1)
    homer_annotatePeaks_arg_label.grid(sticky = "we", row = 50, column = 11, padx = 10, columnspan = 1)
    fold_change_calculator_arg_label.grid(sticky = "we", row = 51, column = 11, padx = 10, columnspan = 1)
    homer_findMotifsGenome_label.grid(sticky = "we", row = 52, column = 11, padx = 10, columnspan = 1)
    meme_chip_label.grid(sticky = "we", row = 53, column = 11, padx = 10, columnspan = 1, pady = (0,10))

    # Entry fields where user can manually key in their input to be stored as a variable
    fastqc1_arg_entry = tk.Entry(setting_table_window, textvariable = fastqc1_arg, width = 50, justify = tk.LEFT)
    clumpify_arg_entry = tk.Entry(setting_table_window, textvariable = clumpify_arg, width = 50, justify = tk.LEFT)
    bbduk_arg_entry = tk.Entry(setting_table_window, textvariable = bbduk_arg, width = 50, justify = tk.LEFT)
    trimmomatic_arg_entry = tk.Entry(setting_table_window, textvariable = trimmomatic_arg, width = 50, justify = tk.LEFT)
    fastqc2_arg_entry = tk.Entry(setting_table_window, textvariable = fastqc2_arg, width = 50, justify = tk.LEFT)
    bwa_mem_arg_entry = tk.Entry(setting_table_window, textvariable = bwa_mem_arg, width = 50, justify = tk.LEFT)
    samtools_view_arg_entry = tk.Entry(setting_table_window, textvariable = samtools_view_arg, width = 50, justify = tk.LEFT)
    plotfingerprint_arg_entry = tk.Entry(setting_table_window, textvariable = plotfingerprint_arg, width = 50, justify = tk.LEFT)
    fastqc3_arg_entry = tk.Entry(setting_table_window, textvariable = fastqc3_arg, width = 50, justify = tk.LEFT)
    macs2_callpeak_arg_entry = tk.Entry(setting_table_window, textvariable = macs2_callpeak_arg, width = 50, justify = tk.LEFT)
    gem_arg_entry = tk.Entry(setting_table_window, textvariable = gem_arg, width = 50, justify = tk.LEFT)
    sicer2_arg_entry = tk.Entry(setting_table_window, textvariable = sicer2_arg, width = 50, justify = tk.LEFT)
    homer_findPeaks_arg_entry = tk.Entry(setting_table_window, textvariable = homer_findPeaks_arg, width = 50, justify = tk.LEFT)
    genrich_arg_entry = tk.Entry(setting_table_window, textvariable = genrich_arg, width = 50, justify = tk.LEFT)
    homer_mergePeaks_arg_entry = tk.Entry(setting_table_window, textvariable = homer_mergePeaks_arg, width = 50, justify = tk.LEFT)
    homer_annotatePeaks_arg_entry = tk.Entry(setting_table_window, textvariable = homer_annotatePeaks_arg, width = 50, justify = tk.LEFT)
    fold_change_calculator_arg_entry = tk.Entry(setting_table_window, textvariable = fold_change_calculator_arg, width = 50, justify = tk.LEFT)
    homer_findMotifsGenome_arg_entry = tk.Entry(setting_table_window, textvariable = homer_findMotifsGenome_arg, width = 50, justify = tk.LEFT)
    meme_chip_arg_entry = tk.Entry(setting_table_window, textvariable = meme_chip_arg, width = 50, justify = tk.LEFT)

    # Display all the entry fields created above
    fastqc1_arg_entry.grid(sticky = "we", row = 35, column = 12, padx = (0,10), columnspan = 8, ipady = 3, pady = (10,0))
    clumpify_arg_entry.grid(sticky = "we", row = 36, column = 12, padx = (0,10), columnspan = 8, ipady = 3)
    bbduk_arg_entry.grid(sticky = "we", row = 37, column = 12, padx = (0,10), columnspan = 8, ipady = 3)
    trimmomatic_arg_entry.grid(sticky = "we", row = 38, column = 12, padx = (0,10), columnspan = 8, ipady = 3)
    fastqc2_arg_entry.grid(sticky = "we", row = 39, column = 12, padx = (0,10), columnspan = 8, ipady = 3)
    bwa_mem_arg_entry.grid(sticky = "we", row = 40, column = 12, padx = (0,10), columnspan = 8, ipady = 3)
    samtools_view_arg_entry.grid(sticky = "we", row = 41, column = 12, padx = (0,10), columnspan = 8, ipady = 3)
    plotfingerprint_arg_entry.grid(sticky = "we", row = 42, column = 12, padx = (0,10), columnspan = 8, ipady = 3)
    fastqc3_arg_entry.grid(sticky = "we", row = 43, column = 12, padx = (0,10), columnspan = 8, ipady = 3)
    macs2_callpeak_arg_entry.grid(sticky = "we", row = 44, column = 12, padx = (0,10), columnspan = 8, ipady = 3)
    gem_arg_entry.grid(sticky = "we", row = 45, column = 12, padx = (0,10), columnspan = 8, ipady = 3)
    sicer2_arg_entry.grid(sticky = "we", row = 46, column = 12, padx = (0,10), columnspan = 8, ipady = 3)
    homer_findPeaks_arg_entry.grid(sticky = "we", row = 47, column = 12, padx = (0,10), columnspan = 8, ipady = 3)
    genrich_arg_entry.grid(sticky = "we", row = 48, column = 12, padx = (0,10), columnspan = 8, ipady = 3)
    homer_mergePeaks_arg_entry.grid(sticky = "we", row = 49, column = 12, padx = (0,10), columnspan = 8, ipady = 3)
    homer_annotatePeaks_arg_entry.grid(sticky = "we", row = 50, column = 12, padx = (0,10), columnspan = 8, ipady = 3)
    fold_change_calculator_arg_entry.grid(sticky = "we", row = 51, column = 12, padx = (0,10), columnspan = 8, ipady = 3)
    homer_findMotifsGenome_arg_entry.grid(sticky = "we", row = 52, column = 12, padx = (0,10), columnspan = 8, ipady = 3)
    meme_chip_arg_entry.grid(sticky = "we", row = 53, column = 12, padx = (0,10), columnspan = 8, ipady = 3, pady = (0,10))


    # Clickable button that triggers a function when activated
    window_default_setting_button = tk.Button(setting_table_window, text = 'Restore defaults', command = lambda:default_setting_button_function(), state = tk.NORMAL, width = 22)
    CreateToolTip(window_default_setting_button, text = 'Click here if you want to restore\nto ChIP-AP pipeline default settings')
    window_default_setting_button.grid(sticky = "w", row = 60, column = 11, padx = 10, pady = (4,10), rowspan = 1, columnspan = 1)

    # Clickable button that triggers a function when activated
    setting_table_window_close_button = tk.Button(setting_table_window, text = 'Accept and close', command = lambda:setting_table_window_close(), state = tk.NORMAL, width = 22)
    CreateToolTip(setting_table_window_close_button, text = 'Click here to accept the settings\nabove and close this window')
    setting_table_window_close_button.grid(sticky = "e", row = 60, column = 12, padx = 10, pady = (4,10), rowspan = 1, columnspan = 8)

    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path


# Function to close the customize settings frame on top of the main frame
def setting_table_window_close():
    current_setting_table_full_path = setting_table_string_var.get()

    # Detect setting values modification by comparing values in the entry fields to the values in the loaded settings table
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

    setting_table_window.destroy() # Close the customize settings frame

# Dynamic label with changeable text
setting_table_notification_label = tk.Label(root, textvariable = setting_table_notification_string_var, width = 70, padx = 10, pady = 4, fg = 'blue')
# setting_table_notification_string_var.set("Currently using default settings table")
setting_table_notification_label.grid(row = 6, column = 11, padx = 10, pady = 4, columnspan = 9)

# Clickable button that triggers a function when activated
setting_table_button = tk.Button(root, text = 'Load setting table', command = lambda:setting_table_button_function(), state = tk.NORMAL, width = 15)
CreateToolTip(setting_table_button, text = 'Click here to browse and\nselect your setting table file')
setting_table_button.grid(sticky = "we", row = 7, column = 11, padx = 10, rowspan = 1, columnspan = 1)

# read_setting_table_function(setting_table_string_var.get())
# Entry field where user can manually key in their input to be stored as a variable
setting_table_entry = tk.Entry(root, textvariable = setting_table_string_var, width = 50, justify = tk.RIGHT)
setting_table_entry.grid(sticky = "we", row = 7, column = 12, padx = (0,10), columnspan = 8, ipady = 3)

# Clickable button that triggers a function when activated
default_setting_button = tk.Button(root, text = 'Default settings', command = lambda:default_setting_button_function(), state = tk.NORMAL, width = 15)
CreateToolTip(default_setting_button, text = 'Click here if you want to restore\nto ChIP-AP pipeline default settings')
default_setting_button.grid(sticky = "we", row = 8, column = 11, padx = 10, pady = (4,5), rowspan = 1, columnspan = 1)

# Clickable button that triggers a function when activated
setting_table_window_open_button = tk.Button(root, text = 'Manual customization >>', command = lambda:setting_table_window_open(), state = tk.NORMAL, width = 25)
CreateToolTip(setting_table_window_open_button, text = "Warning: proceed only when\nyou know what you are doing.\nOtherwise, leave all settings\nat their default values.\nInvalid values will\nlikely cause the whole\npipeline run to break.\nCheck GitHub documentation\nto learn more about\ncustom settings.")
setting_table_window_open_button.grid(sticky = "e", row = 8, column = 12, padx = 10, pady = (4,5), rowspan = 1, columnspan = 8)

# Horizontal separator line
ttk.Separator(root, style = "Line.TSeparator", orient = tk.HORIZONTAL).grid(row = 10, column = 10, columnspan = 11, sticky = "we")


########################################################################################################################


# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def chip_rep1_r1_button_function():
    chip_rep1_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), mode = 'rb', title = 'Choose a file')
    if chip_rep1_r1_button_input: # If the user choose a file from the browser
        chip_rep1_r1_string_var.set(chip_rep1_r1_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(chip_rep1_r1_button_input.name)) # Set the last visited directory to the user-chosen file directory

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def chip_rep1_r2_button_function():
    chip_rep1_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), mode = 'rb', title = 'Choose a file')
    if chip_rep1_r2_button_input: # If the user choose a file from the browser
        chip_rep1_r2_string_var.set(chip_rep1_r2_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(chip_rep1_r2_button_input.name)) # Set the last visited directory to the user-chosen file directory

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def chip_rep2_r1_button_function():
    chip_rep2_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), mode = 'rb', title = 'Choose a file')
    if chip_rep2_r1_button_input: # If the user choose a file from the browser
        chip_rep2_r1_string_var.set(chip_rep2_r1_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(chip_rep2_r1_button_input.name)) # Set the last visited directory to the user-chosen file directory

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def chip_rep2_r2_button_function():
    chip_rep2_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), mode = 'rb', title = 'Choose a file')
    if chip_rep2_r2_button_input: # If the user choose a file from the browser
        chip_rep2_r2_string_var.set(chip_rep2_r2_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(chip_rep2_r2_button_input.name)) # Set the last visited directory to the user-chosen file directory

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def chip_rep3_r1_button_function():
    chip_rep3_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), mode = 'rb', title = 'Choose a file')
    if chip_rep3_r1_button_input: # If the user choose a file from the browser
        chip_rep3_r1_string_var.set(chip_rep3_r1_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(chip_rep3_r1_button_input.name)) # Set the last visited directory to the user-chosen file directory

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def chip_rep3_r2_button_function():
    chip_rep3_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), mode = 'rb', title = 'Choose a file')
    if chip_rep3_r2_button_input: # If the user choose a file from the browser
        chip_rep3_r2_string_var.set(chip_rep3_r2_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(chip_rep3_r2_button_input.name)) # Set the last visited directory to the user-chosen file directory

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def chip_rep4_r1_button_function():
    chip_rep4_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), mode = 'rb', title = 'Choose a file')
    if chip_rep4_r1_button_input: # If the user choose a file from the browser
        chip_rep4_r1_string_var.set(chip_rep4_r1_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(chip_rep4_r1_button_input.name)) # Set the last visited directory to the user-chosen file directory

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def chip_rep4_r2_button_function():
    chip_rep4_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), mode = 'rb', title = 'Choose a file')
    if chip_rep4_r2_button_input: # If the user choose a file from the browser
        chip_rep4_r2_string_var.set(chip_rep4_r2_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(chip_rep4_r2_button_input.name)) # Set the last visited directory to the user-chosen file directory

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def chip_rep5_r1_button_function():
    chip_rep5_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), mode = 'rb', title = 'Choose a file')
    if chip_rep5_r1_button_input: # If the user choose a file from the browser
        chip_rep5_r1_string_var.set(chip_rep5_r1_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(chip_rep5_r1_button_input.name)) # Set the last visited directory to the user-chosen file directory

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def chip_rep5_r2_button_function():
    chip_rep5_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), mode = 'rb', title = 'Choose a file')
    if chip_rep5_r2_button_input: # If the user choose a file from the browser
        chip_rep5_r2_string_var.set(chip_rep5_r2_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(chip_rep5_r2_button_input.name)) # Set the last visited directory to the user-chosen file directory


# Clickable button that triggers a function when activated and an entry field where user can manually key in their input to be stored as a variable
# For each of the five (maximum number that can be handled by ChIP-AP) first reads of the ChIP sample replicates, and each of the five second reads 
chip_rep1_r1_button = tk.Button(root, text ='ChIP rep 1 read 1', command = lambda:chip_rep1_r1_button_function(), state = tk.DISABLED, width = 15)
chip_rep1_r1_entry = tk.Entry(root, textvariable = chip_rep1_r1_string_var, width = 50, justify = tk.RIGHT, state = tk.DISABLED)
chip_rep1_r2_button = tk.Button(root, text ='ChIP rep 1 read 2', command = lambda:chip_rep1_r2_button_function(), state = tk.DISABLED, width = 15)
chip_rep1_r2_entry = tk.Entry(root, textvariable = chip_rep1_r2_string_var, width = 50, justify = tk.RIGHT, state = tk.DISABLED)
chip_rep2_r1_button = tk.Button(root, text ='ChIP rep 2 read 1', command = lambda:chip_rep2_r1_button_function(), state = tk.DISABLED, width = 15)
chip_rep2_r1_entry = tk.Entry(root, textvariable = chip_rep2_r1_string_var, width = 50, justify = tk.RIGHT, state = tk.DISABLED)
chip_rep2_r2_button = tk.Button(root, text ='ChIP rep 2 read 2', command = lambda:chip_rep2_r2_button_function(), state = tk.DISABLED, width = 15)
chip_rep2_r2_entry = tk.Entry(root, textvariable = chip_rep2_r2_string_var, width = 50, justify = tk.RIGHT, state = tk.DISABLED)
chip_rep3_r1_button = tk.Button(root, text ='ChIP rep 3 read 1', command = lambda:chip_rep3_r1_button_function(), state = tk.DISABLED, width = 15)
chip_rep3_r1_entry = tk.Entry(root, textvariable = chip_rep3_r1_string_var, width = 50, justify = tk.RIGHT, state = tk.DISABLED)
chip_rep3_r2_button = tk.Button(root, text ='ChIP rep 3 read 2', command = lambda:chip_rep3_r2_button_function(), state = tk.DISABLED, width = 15)
chip_rep3_r2_entry = tk.Entry(root, textvariable = chip_rep3_r2_string_var, width = 50, justify = tk.RIGHT, state = tk.DISABLED)
chip_rep4_r1_button = tk.Button(root, text ='ChIP rep 4 read 1', command = lambda:chip_rep4_r1_button_function(), state = tk.DISABLED, width = 15)
chip_rep4_r1_entry = tk.Entry(root, textvariable = chip_rep4_r1_string_var, width = 50, justify = tk.RIGHT, state = tk.DISABLED)
chip_rep4_r2_button = tk.Button(root, text ='ChIP rep 4 read 2', command = lambda:chip_rep4_r2_button_function(), state = tk.DISABLED, width = 15)
chip_rep4_r2_entry = tk.Entry(root, textvariable = chip_rep4_r2_string_var, width = 50, justify = tk.RIGHT, state = tk.DISABLED)
chip_rep5_r1_button = tk.Button(root, text ='ChIP rep 5 read 1', command = lambda:chip_rep5_r1_button_function(), state = tk.DISABLED, width = 15)
chip_rep5_r1_entry = tk.Entry(root, textvariable = chip_rep5_r1_string_var, width = 50, justify = tk.RIGHT, state = tk.DISABLED)
chip_rep5_r2_button = tk.Button(root, text ='ChIP rep 5 read 2', command = lambda:chip_rep5_r2_button_function(), state = tk.DISABLED, width = 15)
chip_rep5_r2_entry = tk.Entry(root, textvariable = chip_rep5_r2_string_var, width = 50, justify = tk.RIGHT, state = tk.DISABLED)

# Display all ten buttons and ten entry fields created above
chip_rep1_r1_button.grid(sticky = "we", column = 1, row = 12, padx = 10, columnspan = 1, pady = (5,0))
chip_rep1_r2_button.grid(sticky = "we", column = 1, row = 13, padx = 10, columnspan = 1, pady = (0,5))
chip_rep1_r1_entry.grid(sticky = "we", column = 2, row = 12, padx = (0,10), columnspan = 8, ipady = 3, pady = (5,0))
chip_rep1_r2_entry.grid(sticky = "we", column = 2, row = 13, padx = (0,10), columnspan = 8, ipady = 3, pady = (0,5))
chip_rep2_r1_button.grid(sticky = "we", column = 1, row = 15, padx = 10, columnspan = 1)
chip_rep2_r2_button.grid(sticky = "we", column = 1, row = 16, padx = 10, columnspan = 1, pady = (0,5))
chip_rep2_r1_entry.grid(sticky = "we", column = 2, row = 15, padx = (0,10), columnspan = 8, ipady = 3)
chip_rep2_r2_entry.grid(sticky = "we", column = 2, row = 16, padx = (0,10), columnspan = 8, ipady = 3, pady = (0,5))
chip_rep3_r1_button.grid(sticky = "we", column = 1, row = 18, padx = 10, columnspan = 1)
chip_rep3_r2_button.grid(sticky = "we", column = 1, row = 19, padx = 10, columnspan = 1, pady = (0,5))
chip_rep3_r1_entry.grid(sticky = "we", column = 2, row = 18, padx = (0,10), columnspan = 8, ipady = 3)
chip_rep3_r2_entry.grid(sticky = "we", column = 2, row = 19, padx = (0,10), columnspan = 8, ipady = 3, pady = (0,5))
chip_rep4_r1_button.grid(sticky = "we", column = 1, row = 21, padx = 10, columnspan = 1)
chip_rep4_r2_button.grid(sticky = "we", column = 1, row = 22, padx = 10, columnspan = 1, pady = (0,5))
chip_rep4_r1_entry.grid(sticky = "we", column = 2, row = 21, padx = (0,10), columnspan = 8, ipady = 3)
chip_rep4_r2_entry.grid(sticky = "we", column = 2, row = 22, padx = (0,10), columnspan = 8, ipady = 3, pady = (0,5))
chip_rep5_r1_button.grid(sticky = "we", column = 1, row = 24, padx = 10, columnspan = 1)
chip_rep5_r2_button.grid(sticky = "we", column = 1, row = 25, padx = 10, columnspan = 1)
chip_rep5_r1_entry.grid(sticky = "we", column = 2, row = 24, padx = (0,10), columnspan = 8, ipady = 3)
chip_rep5_r2_entry.grid(sticky = "we", column = 2, row = 25, padx = (0,10), columnspan = 8, ipady = 3)

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

# Dynamic label with changeable text
chip_sample_notification_label = tk.Label(root, textvariable = chip_sample_notification_string_var, width = 70, padx = 10, pady = 4)
chip_sample_notification_label.grid(sticky = "we", column = 1, row = 26, padx = 10, pady = 4, columnspan = 9)

# Horizontal separator line
ttk.Separator(root, style = "Line.TSeparator", orient = tk.HORIZONTAL).grid(row = 30, column = 1, columnspan = 11, sticky = "we")


########################################################################################################################


# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def ctrl_rep1_r1_button_function():
    ctrl_rep1_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), mode = 'rb', title = 'Choose a file')
    if ctrl_rep1_r1_button_input: # If the user choose a file from the browser
        ctrl_rep1_r1_string_var.set(ctrl_rep1_r1_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(ctrl_rep1_r1_button_input.name)) # Set the last visited directory to the user-chosen file directory

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def ctrl_rep1_r2_button_function():
    ctrl_rep1_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), mode = 'rb', title = 'Choose a file')
    if ctrl_rep1_r2_button_input: # If the user choose a file from the browser
        ctrl_rep1_r2_string_var.set(ctrl_rep1_r2_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(ctrl_rep1_r2_button_input.name)) # Set the last visited directory to the user-chosen file directory

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def ctrl_rep2_r1_button_function():
    ctrl_rep2_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), mode = 'rb', title = 'Choose a file')
    if ctrl_rep2_r1_button_input: # If the user choose a file from the browser
        ctrl_rep2_r1_string_var.set(ctrl_rep2_r1_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(ctrl_rep2_r1_button_input.name)) # Set the last visited directory to the user-chosen file directory

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def ctrl_rep2_r2_button_function():
    ctrl_rep2_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), mode = 'rb', title = 'Choose a file')
    if ctrl_rep2_r2_button_input: # If the user choose a file from the browser
        ctrl_rep2_r2_string_var.set(ctrl_rep2_r2_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(ctrl_rep2_r2_button_input.name)) # Set the last visited directory to the user-chosen file directory

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def ctrl_rep3_r1_button_function():
    ctrl_rep3_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), mode = 'rb', title = 'Choose a file')
    if ctrl_rep3_r1_button_input: # If the user choose a file from the browser
        ctrl_rep3_r1_string_var.set(ctrl_rep3_r1_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(ctrl_rep3_r1_button_input.name)) # Set the last visited directory to the user-chosen file directory

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def ctrl_rep3_r2_button_function():
    ctrl_rep3_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), mode = 'rb', title = 'Choose a file')
    if ctrl_rep3_r2_button_input: # If the user choose a file from the browser
        ctrl_rep3_r2_string_var.set(ctrl_rep3_r2_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(ctrl_rep3_r2_button_input.name)) # Set the last visited directory to the user-chosen file directory

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def ctrl_rep4_r1_button_function():
    ctrl_rep4_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), mode = 'rb', title = 'Choose a file')
    if ctrl_rep4_r1_button_input: # If the user choose a file from the browser
        ctrl_rep4_r1_string_var.set(ctrl_rep4_r1_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(ctrl_rep4_r1_button_input.name)) # Set the last visited directory to the user-chosen file directory

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def ctrl_rep4_r2_button_function():
    ctrl_rep4_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), mode = 'rb', title = 'Choose a file')
    if ctrl_rep4_r2_button_input: # If the user choose a file from the browser
        ctrl_rep4_r2_string_var.set(ctrl_rep4_r2_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(ctrl_rep4_r2_button_input.name)) # Set the last visited directory to the user-chosen file directory

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def ctrl_rep5_r1_button_function():
    ctrl_rep5_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), mode = 'rb', title = 'Choose a file')
    if ctrl_rep5_r1_button_input: # If the user choose a file from the browser
        ctrl_rep5_r1_string_var.set(ctrl_rep5_r1_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(ctrl_rep5_r1_button_input.name)) # Set the last visited directory to the user-chosen file directory

# Function triggered by the button for manually loading individual sample, to open a browser window for choosing a local file
def ctrl_rep5_r2_button_function():
    ctrl_rep5_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), mode = 'rb', title = 'Choose a file')
    if ctrl_rep5_r2_button_input: # If the user choose a file from the browser
        ctrl_rep5_r2_string_var.set(ctrl_rep5_r2_button_input.name) # Set the sample table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(ctrl_rep5_r2_button_input.name)) # Set the last visited directory to the user-chosen file directory

# Clickable button that triggers a function when activated and an entry field where user can manually key in their input to be stored as a variable
# For each of the five (maximum number that can be handled by ChIP-AP) first reads of the control sample replicates, and each of the five second reads 
ctrl_rep1_r1_button = tk.Button(root, text ='Ctrl rep 1 read 1', command = lambda:ctrl_rep1_r1_button_function(), state = tk.DISABLED, width = 15)
ctrl_rep1_r1_entry = tk.Entry(root, textvariable = ctrl_rep1_r1_string_var, width = 50, justify = tk.RIGHT, state = tk.DISABLED)
ctrl_rep1_r2_button = tk.Button(root, text ='Ctrl rep 1 read 2', command = lambda:ctrl_rep1_r2_button_function(), state = tk.DISABLED, width = 15)
ctrl_rep1_r2_entry = tk.Entry(root, textvariable = ctrl_rep1_r2_string_var, width = 50, justify = tk.RIGHT, state = tk.DISABLED)
ctrl_rep2_r1_button = tk.Button(root, text ='Ctrl rep 2 read 1', command = lambda:ctrl_rep2_r1_button_function(), state = tk.DISABLED, width = 15)
ctrl_rep2_r1_entry = tk.Entry(root, textvariable = ctrl_rep2_r1_string_var, width = 50, justify = tk.RIGHT, state = tk.DISABLED)
ctrl_rep2_r2_button = tk.Button(root, text ='Ctrl rep 2 read 2', command = lambda:ctrl_rep2_r2_button_function(), state = tk.DISABLED, width = 15)
ctrl_rep2_r2_entry = tk.Entry(root, textvariable = ctrl_rep2_r2_string_var, width = 50, justify = tk.RIGHT, state = tk.DISABLED)
ctrl_rep3_r1_button = tk.Button(root, text ='Ctrl rep 3 read 1', command = lambda:ctrl_rep3_r1_button_function(), state = tk.DISABLED, width = 15)
ctrl_rep3_r1_entry = tk.Entry(root, textvariable = ctrl_rep3_r1_string_var, width = 50, justify = tk.RIGHT, state = tk.DISABLED)
ctrl_rep3_r2_button = tk.Button(root, text ='Ctrl rep 3 read 2', command = lambda:ctrl_rep3_r2_button_function(), state = tk.DISABLED, width = 15)
ctrl_rep3_r2_entry = tk.Entry(root, textvariable = ctrl_rep3_r2_string_var, width = 50, justify = tk.RIGHT, state = tk.DISABLED)
ctrl_rep4_r1_button = tk.Button(root, text ='Ctrl rep 4 read 1', command = lambda:ctrl_rep4_r1_button_function(), state = tk.DISABLED, width = 15)
ctrl_rep4_r1_entry = tk.Entry(root, textvariable = ctrl_rep4_r1_string_var, width = 50, justify = tk.RIGHT, state = tk.DISABLED)
ctrl_rep4_r2_button = tk.Button(root, text ='Ctrl rep 4 read 2', command = lambda:ctrl_rep4_r2_button_function(), state = tk.DISABLED, width = 15)
ctrl_rep4_r2_entry = tk.Entry(root, textvariable = ctrl_rep4_r2_string_var, width = 50, justify = tk.RIGHT, state = tk.DISABLED)
ctrl_rep5_r1_button = tk.Button(root, text ='Ctrl rep 5 read 1', command = lambda:ctrl_rep5_r1_button_function(), state = tk.DISABLED, width = 15)
ctrl_rep5_r1_entry = tk.Entry(root, textvariable = ctrl_rep5_r1_string_var, width = 50, justify = tk.RIGHT, state = tk.DISABLED)
ctrl_rep5_r2_button = tk.Button(root, text ='Ctrl rep 5 read 2', command = lambda:ctrl_rep5_r2_button_function(), state = tk.DISABLED, width = 15)
ctrl_rep5_r2_entry = tk.Entry(root, textvariable = ctrl_rep5_r2_string_var, width = 50, justify = tk.RIGHT, state = tk.DISABLED)

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

# Display all ten buttons and ten entry fields created above
ctrl_rep1_r1_button.grid(sticky = "we", column = 11, row = 12, padx = 10, columnspan = 1, pady = (5,0))
ctrl_rep1_r2_button.grid(sticky = "we", column = 11, row = 13, padx = 10, columnspan = 1, pady = (0,5))
ctrl_rep1_r1_entry.grid(sticky = "we", column = 12, row = 12, padx = (0,10), columnspan = 8, ipady = 3, pady = (5,0))
ctrl_rep1_r2_entry.grid(sticky = "we", column = 12, row = 13, padx = (0,10), columnspan = 8, ipady = 3, pady = (0,5))
ctrl_rep2_r1_button.grid(sticky = "we", column = 11, row = 15, padx = 10, columnspan = 1)
ctrl_rep2_r2_button.grid(sticky = "we", column = 11, row = 16, padx = 10, columnspan = 1, pady = (0,5))
ctrl_rep2_r1_entry.grid(sticky = "we", column = 12, row = 15, padx = (0,10), columnspan = 8, ipady = 3)
ctrl_rep2_r2_entry.grid(sticky = "we", column = 12, row = 16, padx = (0,10), columnspan = 8, ipady = 3, pady = (0,5))
ctrl_rep3_r1_button.grid(sticky = "we", column = 11, row = 18, padx = 10, columnspan = 1)
ctrl_rep3_r2_button.grid(sticky = "we", column = 11, row = 19, padx = 10, columnspan = 1, pady = (0,5))
ctrl_rep3_r1_entry.grid(sticky = "we", column = 12, row = 18, padx = (0,10), columnspan = 8, ipady = 3)
ctrl_rep3_r2_entry.grid(sticky = "we", column = 12, row = 19, padx = (0,10), columnspan = 8, ipady = 3, pady = (0,5))
ctrl_rep4_r1_button.grid(sticky = "we", column = 11, row = 21, padx = 10, columnspan = 1)
ctrl_rep4_r2_button.grid(sticky = "we", column = 11, row = 22, padx = 10, columnspan = 1, pady = (0,5))
ctrl_rep4_r1_entry.grid(sticky = "we", column = 12, row = 21, padx = (0,10), columnspan = 8, ipady = 3)
ctrl_rep4_r2_entry.grid(sticky = "we", column = 12, row = 22, padx = (0,10), columnspan = 8, ipady = 3, pady = (0,5))
ctrl_rep5_r1_button.grid(sticky = "we", column = 11, row = 24, padx = 10, columnspan = 1)
ctrl_rep5_r2_button.grid(sticky = "we", column = 11, row = 25, padx = 10, columnspan = 1)
ctrl_rep5_r1_entry.grid(sticky = "we", column = 12, row = 24, padx = (0,10), columnspan = 8, ipady = 3)
ctrl_rep5_r2_entry.grid(sticky = "we", column = 12, row = 25, padx = (0,10), columnspan = 8, ipady = 3)

# Dynamic label with changeable text
ctrl_sample_notification_label = tk.Label(root, textvariable = ctrl_sample_notification_string_var, width = 70, padx = 10, pady = 4)
ctrl_sample_notification_label.grid(sticky = "we", column = 11, row = 26, padx = 10, pady = 4, columnspan = 9)

# Horizontal separator line
ttk.Separator(root, style = "Line.TSeparator", orient = tk.HORIZONTAL).grid(row = 30, column = 11, columnspan = 10, sticky = "we")


########################################################################################################################


# Function to auto-assign the genome directory to the ChIP-AP's default when the user-chosen reference genome is one of the supported builds
# Opens up the button and entry field widgets for inputting by user only when user chooses "others" reference genome, which genome files needs to be self-provided by user  
def auto_assign_genome_folder_function():
    if genome_ref_string_var.get() == "hg38 (Homo sapiens)":
        genome_folder_string_var.set(genome_folder_full_path)
        genome_folder_button.config(state = tk.DISABLED)
        genome_folder_entry.config(state = tk.DISABLED)
    elif genome_ref_string_var.get() == "hg19 (Homo sapiens)":
        genome_folder_string_var.set(genome_folder_full_path)
        genome_folder_button.config(state = tk.DISABLED)
        genome_folder_entry.config(state = tk.DISABLED)
    elif genome_ref_string_var.get() == "mm9 (Mus musculus)":
        genome_folder_string_var.set(genome_folder_full_path)
        genome_folder_button.config(state = tk.DISABLED)
        genome_folder_entry.config(state = tk.DISABLED)
    elif genome_ref_string_var.get() == "mm10 (Mus musculus)":
        genome_folder_string_var.set(genome_folder_full_path)
        genome_folder_button.config(state = tk.DISABLED)
        genome_folder_entry.config(state = tk.DISABLED)
    elif genome_ref_string_var.get() == "dm6 (Drosophila melanogaster)":
        genome_folder_string_var.set(genome_folder_full_path)
        genome_folder_button.config(state = tk.DISABLED)
        genome_folder_entry.config(state = tk.DISABLED)
    elif genome_ref_string_var.get() == "sacCer3 (Saccharomyces cerevisiae)":
        genome_folder_string_var.set(genome_folder_full_path)
        genome_folder_button.config(state = tk.DISABLED)
        genome_folder_entry.config(state = tk.DISABLED)
    else:
        genome_folder_string_var.set('')
        genome_folder_button.config(state = tk.NORMAL)
        genome_folder_entry.config(state = tk.NORMAL)

    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path


# Simple label with static text
genome_ref_label = tk.Label(root, text = "Reference genome:")
genome_ref_label.grid(row = 35, column = 1, padx = 10, columnspan = 1, pady = (5,0))

# Drop-down options widget that triggers a function when an option is chosen
genome_ref_drop_down = tk.OptionMenu(root, genome_ref_string_var, *genome_ref_options, command = lambda x = None : auto_assign_genome_folder_function())
CreateToolTip(genome_ref_drop_down, text = 'Select the genome assembly you want\nyour ChIP-seq reads to be aligned to.\nCurrently, ChIP-AP supports\nsix genome assemblies.\nOutside those supported by\nChIP-AP, you will need to\ngenerate the files yourself')
genome_ref_drop_down.config(takefocus = 1)
genome_ref_drop_down.grid(sticky = "we", row = 35, column = 2, columnspan = 4, pady = (5,0), padx = (0,10))


########################################################################################################################


# Function triggered by the "Genome folder" button, to open a browser window for choosing a local directory
def genome_folder_button_function():
    genome_folder_button_input = filedialog.askdirectory(initialdir = current_dir_string_var.get(), parent = root, title = 'Choose the Genome Directory')
    if genome_folder_button_input: # If the user choose a directory from the browser
        genome_folder_string_var.set(genome_folder_button_input) # Set the setting table GUI variable with the path to the user-chosen directory
        current_dir_string_var.set(genome_folder_button_input) # Set the last visited directory to the user-chosen file directory

    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path


# Clickable button that triggers a function when activated
genome_folder_button = tk.Button(root, text = 'Genome folder', command = lambda:genome_folder_button_function(), state = tk.DISABLED, width = 15)
CreateToolTip(genome_folder_button, text = 'Click here to browse and select\nthe directory containing your\ncustom genome reference')
genome_folder_button.grid(sticky = "we", row = 37, column = 1, padx = 10, columnspan = 1)

# Entry field where user can manually key in their input to be stored as a variable
genome_folder_entry = tk.Entry(root, textvariable = genome_folder_string_var, width = 40, justify = tk.RIGHT, state = tk.DISABLED)
genome_folder_entry.grid(sticky = "we", row = 37, column = 2, columnspan = 4, ipady = 3, padx = (0,10))


########################################################################################################################


# Function triggered by the "Known motif file" button, to open a browser window for choosing a local file
def known_motif_button_function():
    known_motif_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a .motif file')
    if known_motif_button_input: # If the user choose a file from the browser
        known_motif_string_var.set(known_motif_button_input.name) # Set the setting table GUI variable with the path to the user-chosen file
        current_dir_string_var.set(os.path.dirname(known_motif_button_input.name)) # Set the last visited directory to the user-chosen file directory

    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path


# Clickable button that triggers a function when activated
known_motif_button = tk.Button(root, text = 'Known motif file', command = lambda:known_motif_button_function(), state = tk.NORMAL, width = 15)
CreateToolTip(known_motif_button, text = 'Click here to browse and select your\n.motif file (in HOMER matrix format)')
known_motif_button.grid(sticky = "we", row = 39, column = 1, padx = 10, columnspan = 1)

# Entry field where user can manually key in their input to be stored as a variable
known_motif_entry = tk.Entry(root, textvariable = known_motif_string_var, width = 40, justify = tk.RIGHT)
known_motif_entry.grid(sticky = "we", row = 39, column = 2, columnspan = 4, ipady = 3, padx = (0,10))


########################################################################################################################


# Function triggered by the "Output save folder" button, to open a browser window for choosing a local directory
def output_folder_button_function():
    output_folder_button_input = filedialog.askdirectory(initialdir = current_dir_string_var.get(), parent = root, title = 'Choose the Output Directory')
    if output_folder_button_input: # If the user choose a directory from the browser
        output_folder_string_var.set(output_folder_button_input) # Set the setting table GUI variable with the path to the user-chosen directory
        current_dir_string_var.set(output_folder_button_input) # Set the last visited directory to the user-chosen file directory

    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path


# Clickable button that triggers a function when activated
output_folder_button = tk.Button(root, text = 'Output save folder', command = lambda:output_folder_button_function(), state = tk.NORMAL, width = 15)
CreateToolTip(output_folder_button, text = 'Click here to browse and select\nthe directory to save the results\nof your ChIP-AP pipeline run.\nTo save in a new folder,\ntype in the entry field\nyour new folder name after\nthe output save directory:\nfull_path_to_output_save_directory/\nnew_folder_name')
output_folder_button.grid(sticky = "we", column = 1, row = 41, padx = 10, columnspan = 1)

# Entry field where user can manually key in their input to be stored as a variable
output_folder_entry = tk.Entry(root, textvariable = output_folder_string_var, width = 40, justify = tk.RIGHT)
output_folder_entry.grid(sticky = "we", column = 2, row = 41, columnspan = 4, ipady = 3, padx = (0,10))


########################################################################################################################


# Simple label with static text
setname_label = tk.Label(root, text = "Dataset name:", width = 15)
setname_label.grid(row = 43, column = 1, padx = 10, columnspan = 1, pady = (0,5))

# Entry field where user can manually key in their input to be stored as a variable
setname_entry = tk.Entry(root, textvariable = setname_string_var, width = 40, justify = tk.LEFT)
CreateToolTip(setname_entry, text = 'Type in your folder name and prefix\nfor all the resulting output filenames')
setname_entry.grid(sticky = "we", row = 43, column = 2, columnspan = 4, ipady = 3, pady = (0,5), padx = (0,10))


########################################################################################################################


# Simple label with static text
homer_motif_label = tk.Label(root, text = "HOMER motif enrichment:", anchor = "w")
homer_motif_label.grid(row = 35, column = 7, padx = 5, columnspan = 5, sticky = "sw", pady = (5,0))

# Drop-down options widget that triggers a function when an option is chosen
homer_motif_drop_down = tk.OptionMenu(root, homer_motif_string_var, *homer_motif_options)
CreateToolTip(homer_motif_drop_down, text = 'Select the peak set(s) you want\nHOMER findMotifsGenome to perform\nmotif enrichment analysis on')
homer_motif_drop_down.config(takefocus = 1)
homer_motif_drop_down.grid(sticky = "we", row = 37, column = 7, padx = 5, columnspan = 5)


########################################################################################################################


# Simple label with static text
meme_motif_label = tk.Label(root, text = "MEME motif enrichment:", anchor = "w")
meme_motif_label.grid(row = 39, column = 7, padx = 5, columnspan = 5, sticky = "sw")

# Drop-down options widget that triggers a function when an option is chosen
meme_motif_drop_down = tk.OptionMenu(root, meme_motif_string_var, *meme_motif_options)
CreateToolTip(meme_motif_drop_down, text = 'Select the peak set(s) you\nwant meme-chip to perform\nmotif enrichment analysis on')
meme_motif_drop_down.config(takefocus = 1)
meme_motif_drop_down.grid(sticky = "we", row = 41, column = 7, padx = 5, columnspan = 5)


########################################################################################################################


# Checkbox to toggle between on and off state what triggers a function when it is either toggled on or off.
fcmerge_checkbox = tk.Checkbutton(root, text = ' Merged fold enrichment analysis', variable = fcmerge_var, onvalue = 1, offvalue = 0)
CreateToolTip(fcmerge_checkbox, text = 'Check this box if you want\nthe fold enrichment analysis\nfrom all replicates combined as one.\nThis option will be ignored\nwhen there are unequal\nnumber of replicates between\nChIP and control samples')
fcmerge_checkbox.grid(sticky = "w", row = 35, column = 13, columnspan = 8, padx = (0,10), pady = (5,0))


########################################################################################################################


# Checkbox to toggle between on and off state what triggers a function when it is either toggled on or off.
goann_checkbox = tk.Checkbutton(root, text = ' Annotate peaks with known gene ontology terms', variable = goann_var, onvalue = 1, offvalue = 0)
CreateToolTip(goann_checkbox, text = 'Check this box if you want\neach peak in the final peaks list\nto have gene ontology annotations')
goann_checkbox.grid(sticky = "w", row = 37, column = 13, columnspan = 8, padx = (0,10))


########################################################################################################################


# Checkbox to toggle between on and off state what triggers a function when it is either toggled on or off.
pathann_checkbox = tk.Checkbutton(root, text = ' Annotate peaks with known pathway terms', variable = pathann_var, onvalue = 1, offvalue = 0)
CreateToolTip(pathann_checkbox, text = 'Check this box if you want\neach peak in the final peaks list\nto have pathway annotations')
pathann_checkbox.grid(sticky = "w", row = 39, column = 13, columnspan = 8, padx = (0,10))


########################################################################################################################


# Checkbox to toggle between on and off state what triggers a function when it is either toggled on or off.
deltemp_checkbox = tk.Checkbutton(root, text = ' Delete large temporary files', variable = deltemp_var, onvalue = 1, offvalue = 0)
CreateToolTip(deltemp_checkbox, text = 'Check this box if you will not need\nthe large-sized intermediary files')
deltemp_checkbox.grid(sticky = "w", row = 41, column = 13, columnspan = 8, padx = (0,10))


########################################################################################################################


# Checkbox to toggle between on and off state what triggers a function when it is either toggled on or off.
stdout_checkbox = tk.Checkbutton(root, text = ' Record standard outputs & errors', variable = stdout_var, onvalue = 1, offvalue = 0)
CreateToolTip(stdout_checkbox, text = 'Check this box if you want to save\npipeline standard outputs (channel 1>) and\nstandard errors (channel 2>) as text files')
stdout_checkbox.grid(sticky = "w", row = 43, column = 13, columnspan = 8, padx = (0,10), pady = (0,5))

# Horizontal separator line
ttk.Separator(root, style = "Line.TSeparator", orient = tk.HORIZONTAL).grid(row = 46, column = 0, columnspan = 21, sticky = "we")


########################################################################################################################

# Simple label with static text
cpu_count_label = tk.Label(root, text = "CPU cores to use:", width = 20, justify = tk.RIGHT)
cpu_count_label.grid(row = 62, column = 1, padx = 10, columnspan = 1, pady = 5)

# Entry field where user can manually key in their input to be stored as a variable
cpu_count_entry = tk.Entry(root, width = 5, textvariable = cpu_count_string_var)
CreateToolTip(cpu_count_entry, text = 'Type in the number of CPU cores\nto be used by the pipeline')
cpu_count_entry.grid(sticky = "w", row = 62, column = 2, padx = (0,10), ipady = 3, pady = 5)

# Dynamic label with changeable text
cpu_count_notification_label = tk.Label(root, textvariable = cpu_count_notification_string_var, width = 65, anchor = "w")
cpu_count_notification_label.grid(row = 62, column = 3, columnspan = 9, padx = 10)


########################################################################################################################

# Simple label with static text
command_line_label = tk.Label(root, text = 'ChIP-AP command line:', width = 50)
command_line_label.grid(sticky = "we", column = 1, row = 47, padx = 10, columnspan = 19, pady = 5)

command_line_font = tkFont.Font(size = 12) # Define the font size to be used for the displayed full command line

# Dynamic label with changeable text
command_line_output_label = tk.Label(root, textvariable = command_line_output_string_var, relief = tk.SOLID, borderwidth = 1, bg = "white",
                                    width = 50, anchor = "w", justify = tk.LEFT, padx = 10, pady = 4, wraplength = 1450)
command_line_output_label['font'] = command_line_font
command_line_output_label.grid(sticky = "we", column = 1, row = 48, padx = 10, columnspan = 19, ipady = 4)


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

    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path
    print_sample_table_function() # Create a new sample table based on all the samples registered in the GUI variables, and save it in the output save directory
    print_setting_table_function() # Create a new setting table based on all program flags and arguments registered in the GUI variables, and save it in the output save directory
    print_GS_command_line_function() # Create a text file containing the full command line which was used to call for ChIP-AP pipeline run, and save it in the output save directory
    execute_GS_command_line_function() # Generate the pipeline scripts
    cpu_count_notification_string_var.set("Scripts generated! You can close this window now") # Notify the user
    cpu_count_notification_label.config(fg = 'green')


# Function to execute a series of commands when user activates the "Generate scripts and run" button
def generate_scripts_and_run_function():
    if not os.path.exists(output_dir_arg.get()):
        os.makedirs(output_dir_arg.get())
    
    rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path
    print_sample_table_function() # Create a new sample table based on all the samples registered in the GUI variables, and save it in the output save directory
    print_setting_table_function() # Create a new setting table based on all program flags and arguments registered in the GUI variables, and save it in the output save directory
    print_GSaR_command_line_function() # Create a text file containing the full command line which was used to call for ChIP-AP pipeline run, and save it in the output save directory
    cpu_count_notification_string_var.set("Pipeline started! You can close this window now") # Notify the user
    cpu_count_notification_label.config(fg = 'green')
    execute_GSaR_command_line_function() # Start the pipeline immediately


# Clickable button that triggers a function when activated
generate_scripts_button = tk.Button(root, text = 'Generate scripts', command = lambda:generate_scripts_function(), state = tk.DISABLED, width = 18)
CreateToolTip(generate_scripts_button, text = 'Select this if you wish to\nrun the pipeline later\nby executing MASTER_script.sh\nwithin the output save folder')
generate_scripts_button.grid(row = 62, column = 12, columnspan = 9, sticky = "w", padx = 10, pady = 5)

# Clickable button that triggers a function when activated
generate_and_run_scripts_button = tk.Button(root, text = 'Generate and run scripts', command = lambda:generate_scripts_and_run_function(), state = tk.DISABLED, width = 25)
CreateToolTip(generate_and_run_scripts_button, text = 'Select this if you wish to\nrun the pipeline now\nNOTE: It may take up to\nseveral hours depending\non your system')
generate_and_run_scripts_button.grid(row = 62, column = 12, columnspan = 9, sticky = "e", padx = 10, pady = 5)


########################################################################################################################


# Traced variables. The assigned function will be triggered every time the variable's value changes.
# For variables which value only affects the resulting ChIP-AP command line and nothing else, update_command_line_function is used.
# For variables with more complex rules and greater influence over the values of other variables, other functions are used.
# See the comments of all the tracing functions (above) to learn what things are done when any of these traced variables change in value.
read_mode_string_var.trace('w', register_sample_function)
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
sample_table_string_var.trace('w', check_traced_input_function)
sample_table_string_var.trace('w', sample_table_entry_trace_load_function)
setting_table_string_var.trace('w', check_traced_input_function)
setting_table_string_var.trace('w', setting_table_entry_trace_load_function)
genome_ref_string_var.trace('w', check_traced_input_function)
genome_folder_string_var.trace('w', check_traced_input_function)
known_motif_string_var.trace('w', check_traced_input_function)
setname_string_var.trace('w', check_traced_input_function)
output_folder_string_var.trace('w', check_traced_input_function)
fcmerge_var.trace('w', update_command_line_function)
goann_var.trace('w', update_command_line_function)
pathann_var.trace('w', update_command_line_function)
deltemp_var.trace('w', update_command_line_function)
stdout_var.trace('w', update_command_line_function)
homer_motif_string_var.trace('w', update_command_line_function)
meme_motif_string_var.trace('w', update_command_line_function)
cpu_count_string_var.trace('w', check_traced_input_function)

# Initial run of all variable tracing functions at GUI start, just right after everything finished loading.
# Traces all the initial stored values at GUI start, so GUI text notifications can inform user of what the GUI needs from the very beginning.
register_sample_function()
check_required_input_function() # Check if the required inputs (type B) other than samples (type A) are satisfied or contain any errors
update_command_line_function() # Update the generated ChIP-AP command line based on the latest registered variable values
rearrange_entry_field_function() # Rearrange the entry fields such that they always display the last n characters of entered path

root.mainloop() # Initiates the GUI looping routine


########################################################################################################################
#!/usr/bin/env python3
#pyright: reportUnboundVariable=false

from sys import stdout
import tkinter as tk
import tkinter.font as tkFont
from tkinter import filedialog
from tkinter import ttk
import os
import time
import multiprocessing
import subprocess
# try:
#     import pandas as pd
# except ImportError:
#     print("Dependency missing, downloading and installing now")
#     import pip
#     pip.main(['install', '--user', 'pandas'])
#     time.sleep(5) # Sleep for 3 seconds

import pandas as pd

current_pipeline = 'chipap_v4.1.py'
default_current_dir = os.path.expanduser('~')
genome_folder_full_path = os.path.expanduser('~/genomes')
setting_table_file_full_path = '{}/default_settings_table'.format(genome_folder_full_path)

root = tk.Tk()
root.title(current_pipeline)
root.resizable(width = False, height = False)

default_font = tkFont.nametofont("TkDefaultFont")
default_font.configure(family = 'fixed', size = 20)

text_font = tkFont.nametofont("TkTextFont")
text_font.configure(family = 'fixed', size = 20)

fixed_font = tkFont.nametofont("TkFixedFont")
fixed_font.configure(family = 'fixed', size = 20)

valid_extension_list = ['.fastq', '.fq', '.fastq.gz', '.fq.gz', '.bam']

suite_program_list = [
    'fastqc1',
    'clumpify',
    'bbduk',
    'trimmomatic',
    'fastqc2',
    'bwa_mem',
    'samtools_view',
    'plotfingerprint',
    'fastqc3',
    'macs2_callpeak',
    'gem',
    'sicer2',
    'homer_findPeaks',
    'genrich',
    'homer_mergePeaks',
    'homer_annotatePeaks',
    'fold_change_calculator']

genome_ref_options = ["hg38 (Homo sapiens build 38)", 
                        "hg19 (Homo sapiens build 19)", 
                        "mm9 (Mus musculus build 9)", 
                        "mm10 (Mus musculus build 10)", 
                        "dm6 (Drosophila melanogaster build 6)", 
                        "sacCer3 (Saccharomyces cerevisiae build 3)",
                        "other [!!!under construction!!!]"]


########################################################################################################################


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


def change_frame_function(from_frame, to_frame):
    if from_frame == sample_table_frame and to_frame == chip_list_frame and sample_table_question_string_var.get() == 'yes':
        if sample_table_loading_test_function()  == False:
            return
        elif sample_table_loading_test_function() == True:
            pass

    if from_frame == setting_table_frame and to_frame == setting_value_frame and sample_table_question_string_var.get() == 'yes':
        if setting_table_loading_test_function() == False:
            return
        elif setting_table_loading_test_function() == True:
            pass
    
    from_frame.grid_remove()
    to_frame.grid(row = 0, column = 0, sticky = "news")
    root.eval('tk::PlaceWindow . center')
    root.grab_set()
    root.focus_force()

    if to_frame == chip_list_frame:
        display_chip_widget_function()

    if to_frame == ctrl_list_frame:
        display_ctrl_widget_function()

    root.eval('tk::PlaceWindow . center')
    root.grab_set()
    root.focus_force()
    
    rearrange_entry_field_function()


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


def register_sample_function(*args):
        
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

    global chip_list_r1
    global chip_list_r2
    global ctrl_list_r1
    global ctrl_list_r2
    
    chip_list_r1 = []
    chip_list_r2 = []
    ctrl_list_r1 = []
    ctrl_list_r2 = []

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

    global sample_table_output_dict
    
    sample_table_output_dict = {'chip_read_1' : chip_list_r1,
                                'chip_read_2' : chip_list_r2,
                                'ctrl_read_1' : ctrl_list_r1,
                                'ctrl_read_2' : ctrl_list_r2}

    check_chip_sample_assigned_function()
    check_ctrl_sample_assigned_function()
    check_required_input_function()
    update_command_line_function()
    

def check_chip_sample_assigned_function():

    if read_mode_string_var.get() == 'single' or read_mode_string_var.get() == 'paired':
        
        if read_mode_string_var.get() == 'single':
            if all(chip_r1 == '' for chip_r1 in chip_list_r1):
                chip_sample_notification_string_var.set('Please assign ChIP sample')
                chip_sample_notification_label.config(fg = 'blue')
                chip_list_frame_continue_button.config(state = tk.DISABLED)
                return

        if read_mode_string_var.get() == 'paired':
            if all(chip_r1 == '' for chip_r1 in chip_list_r1):
                chip_sample_notification_string_var.set('Please assign ChIP sample (read 1)')
                chip_sample_notification_label.config(fg = 'blue')
                chip_list_frame_continue_button.config(state = tk.DISABLED)
                return

            # Check for r2 only if the input files are not aligned
            if all(chip_r2 == '' for chip_r2 in chip_list_r2) and not all('.bam' in sample for sample in (chip_list_r1 + chip_list_r1) if sample != ''):
                chip_sample_notification_string_var.set('Please assign ChIP sample (read 2)')
                chip_sample_notification_label.config(fg = 'blue')
                chip_list_frame_continue_button.config(state = tk.DISABLED)
                return

        if read_mode_string_var.get() == 'single':
            if not bool(chip_rep1_r1_string_var.get()):
                chip_sample_notification_string_var.set('ChIP samples are not assigned from the top (replicate 1)')
                chip_sample_notification_label.config(fg = 'red')
                chip_rep1_r1_entry.config(bg = 'IndianRed1')
                chip_list_frame_continue_button.config(state = tk.DISABLED)
                return
        
        if read_mode_string_var.get() == 'paired':
            if not bool(chip_rep1_r1_string_var.get()) and not bool(chip_rep1_r2_string_var.get()):
                chip_sample_notification_string_var.set('ChIP samples are not assigned from the top (replicate 1)')
                chip_sample_notification_label.config(fg = 'red')
                chip_rep1_r1_entry.config(bg = 'IndianRed1')
                chip_rep1_r2_entry.config(bg = 'IndianRed1')
                chip_list_frame_continue_button.config(state = tk.DISABLED)
                return

            if not bool(chip_rep1_r1_string_var.get()):
                chip_sample_notification_string_var.set('ChIP samples are not assigned from the top (replicate 1)')
                chip_sample_notification_label.config(fg = 'red')
                chip_rep1_r1_entry.config(bg = 'IndianRed1')
                chip_list_frame_continue_button.config(state = tk.DISABLED)
                return

            if not bool(chip_rep1_r2_string_var.get()):
                chip_sample_notification_string_var.set('ChIP samples are not assigned from the top (replicate 1)')
                chip_sample_notification_label.config(fg = 'red')
                chip_rep1_r2_entry.config(bg = 'IndianRed1')
                chip_list_frame_continue_button.config(state = tk.DISABLED)
                return

        chip_sample_notification_string_var.set('ChIP samples have been assigned')
        chip_sample_notification_label.config(fg = 'green')
        chip_list_frame_continue_button.config(state = tk.NORMAL)
    
        check_chip_sample_file_pair_function()

    else:
        return


def check_ctrl_sample_assigned_function():

    if read_mode_string_var.get() == 'single' or read_mode_string_var.get() == 'paired':
        
        if read_mode_string_var.get() == 'single':
            if all(ctrl_r1 == '' for ctrl_r1 in ctrl_list_r1):
                ctrl_sample_notification_string_var.set('Please assign control sample')
                ctrl_sample_notification_label.config(fg = 'blue')
                ctrl_list_frame_continue_button.config(state = tk.DISABLED)
                return

        if read_mode_string_var.get() == 'paired':
            if all(ctrl_r1 == '' for ctrl_r1 in ctrl_list_r1):
                ctrl_sample_notification_string_var.set('Please assign control sample (read 1)')
                ctrl_sample_notification_label.config(fg = 'blue')
                ctrl_list_frame_continue_button.config(state = tk.DISABLED)
                return

            # Check for r2 only if the input files are not aligned
            if all(ctrl_r2 == '' for ctrl_r2 in ctrl_list_r2) and not all('.bam' in sample for sample in (ctrl_list_r1 + ctrl_list_r1) if sample != ''):
                ctrl_sample_notification_string_var.set('Please assign control sample (read 2)')
                ctrl_sample_notification_label.config(fg = 'blue')
                ctrl_list_frame_continue_button.config(state = tk.DISABLED)
                return

        if read_mode_string_var.get() == 'single':
            if not bool(ctrl_rep1_r1_string_var.get()):
                ctrl_sample_notification_string_var.set('Control samples are not assigned from the top (replicate 1)')
                ctrl_sample_notification_label.config(fg = 'red')
                ctrl_rep1_r1_entry.config(bg = 'IndianRed1')
                ctrl_list_frame_continue_button.config(state = tk.DISABLED)
                return
        
        if read_mode_string_var.get() == 'paired':
            if not bool(ctrl_rep1_r1_string_var.get()) and not bool(ctrl_rep1_r2_string_var.get()):
                ctrl_sample_notification_string_var.set('Control samples are not assigned from the top (replicate 1)')
                ctrl_sample_notification_label.config(fg = 'red')
                ctrl_rep1_r1_entry.config(bg = 'IndianRed1')
                ctrl_rep1_r2_entry.config(bg = 'IndianRed1')
                ctrl_list_frame_continue_button.config(state = tk.DISABLED)
                return

            if not bool(ctrl_rep1_r1_string_var.get()):
                ctrl_sample_notification_string_var.set('Control samples are not assigned from the top (replicate 1)')
                ctrl_sample_notification_label.config(fg = 'red')
                ctrl_rep1_r1_entry.config(bg = 'IndianRed1')
                ctrl_list_frame_continue_button.config(state = tk.DISABLED)
                return

            if not bool(ctrl_rep1_r2_string_var.get()):
                ctrl_sample_notification_string_var.set('Control samples are not assigned from the top (replicate 1)')
                ctrl_sample_notification_label.config(fg = 'red')
                ctrl_rep1_r2_entry.config(bg = 'IndianRed1')
                ctrl_list_frame_continue_button.config(state = tk.DISABLED)
                return

        ctrl_sample_notification_string_var.set('Control samples have been assigned')
        ctrl_sample_notification_label.config(fg = 'green')
        ctrl_list_frame_continue_button.config(state = tk.NORMAL)
    
        check_ctrl_sample_file_pair_function()

    else:
        return


def check_chip_sample_file_pair_function():

    if read_mode_string_var.get() == 'paired':
    
        chip_pair_error_state = 0

        if bool(chip_rep1_r1_string_var.get()) != bool(chip_rep1_r2_string_var.get()):
            chip_pair_error_state = 1
            
            if not bool(chip_rep1_r1_string_var.get()):
                chip_rep1_r1_entry.config(bg = 'IndianRed1')
            
            if not bool(chip_rep1_r2_string_var.get()):
                chip_rep1_r2_entry.config(bg = 'IndianRed1')

        if bool(chip_rep2_r1_string_var.get()) != bool(chip_rep2_r2_string_var.get()):
            chip_pair_error_state = 1

            if not bool(chip_rep2_r1_string_var.get()):
                chip_rep2_r1_entry.config(bg = 'IndianRed1')
            
            if not bool(chip_rep2_r2_string_var.get()):
                chip_rep2_r2_entry.config(bg = 'IndianRed1')

        if bool(chip_rep3_r1_string_var.get()) != bool(chip_rep3_r2_string_var.get()):
            chip_pair_error_state = 1

            if not bool(chip_rep3_r1_string_var.get()):
                chip_rep3_r1_entry.config(bg = 'IndianRed1')
            
            if not bool(chip_rep3_r2_string_var.get()):
                chip_rep3_r2_entry.config(bg = 'IndianRed1')

        if bool(chip_rep4_r1_string_var.get()) != bool(chip_rep4_r2_string_var.get()):
            chip_pair_error_state = 1

            if not bool(chip_rep4_r1_string_var.get()):
                chip_rep4_r1_entry.config(bg = 'IndianRed1')
            
            if not bool(chip_rep4_r2_string_var.get()):
                chip_rep4_r2_entry.config(bg = 'IndianRed1')

        if bool(chip_rep5_r1_string_var.get()) != bool(chip_rep5_r2_string_var.get()):
            chip_pair_error_state = 1

            if not bool(chip_rep5_r1_string_var.get()):
                chip_rep5_r1_entry.config(bg = 'IndianRed1')
            
            if not bool(chip_rep5_r2_string_var.get()):
                chip_rep5_r2_entry.config(bg = 'IndianRed1')
    
        if chip_pair_error_state == 1:
            chip_sample_notification_string_var.set('One or more of ChIP samples have missing pair')
            chip_sample_notification_label.config(fg = 'red')
            chip_list_frame_continue_button.config(state = tk.DISABLED)
            return

        elif chip_pair_error_state == 0:
            chip_sample_notification_string_var.set('All ChIP samples are properly paired')
            chip_sample_notification_label.config(fg = 'green')
        
        check_chip_sample_validity_function()

    else:
        check_chip_sample_validity_function()


def check_ctrl_sample_file_pair_function():

    if read_mode_string_var.get() == 'paired':
    
        ctrl_pair_error_state = 0

        if bool(ctrl_rep1_r1_string_var.get()) != bool(ctrl_rep1_r2_string_var.get()):
            ctrl_pair_error_state = 1
            
            if not bool(ctrl_rep1_r1_string_var.get()):
                ctrl_rep1_r1_entry.config(bg = 'IndianRed1')
            
            if not bool(ctrl_rep1_r2_string_var.get()):
                ctrl_rep1_r2_entry.config(bg = 'IndianRed1')

        if bool(ctrl_rep2_r1_string_var.get()) != bool(ctrl_rep2_r2_string_var.get()):
            ctrl_pair_error_state = 1

            if not bool(ctrl_rep2_r1_string_var.get()):
                ctrl_rep2_r1_entry.config(bg = 'IndianRed1')
            
            if not bool(ctrl_rep2_r2_string_var.get()):
                ctrl_rep2_r2_entry.config(bg = 'IndianRed1')

        if bool(ctrl_rep3_r1_string_var.get()) != bool(ctrl_rep3_r2_string_var.get()):
            ctrl_pair_error_state = 1

            if not bool(ctrl_rep3_r1_string_var.get()):
                ctrl_rep3_r1_entry.config(bg = 'IndianRed1')
            
            if not bool(ctrl_rep3_r2_string_var.get()):
                ctrl_rep3_r2_entry.config(bg = 'IndianRed1')

        if bool(ctrl_rep4_r1_string_var.get()) != bool(ctrl_rep4_r2_string_var.get()):
            ctrl_pair_error_state = 1

            if not bool(ctrl_rep4_r1_string_var.get()):
                ctrl_rep4_r1_entry.config(bg = 'IndianRed1')
            
            if not bool(ctrl_rep4_r2_string_var.get()):
                ctrl_rep4_r2_entry.config(bg = 'IndianRed1')

        if bool(ctrl_rep5_r1_string_var.get()) != bool(ctrl_rep5_r2_string_var.get()):
            ctrl_pair_error_state = 1

            if not bool(ctrl_rep5_r1_string_var.get()):
                ctrl_rep5_r1_entry.config(bg = 'IndianRed1')
            
            if not bool(ctrl_rep5_r2_string_var.get()):
                ctrl_rep5_r2_entry.config(bg = 'IndianRed1')
    
        if ctrl_pair_error_state == 1:
            ctrl_sample_notification_string_var.set('One or more of control samples have missing pair')
            ctrl_sample_notification_label.config(fg = 'red')
            ctrl_list_frame_continue_button.config(state = tk.DISABLED)
            return

        elif ctrl_pair_error_state == 0:
            ctrl_sample_notification_string_var.set('All control samples are properly paired')
            ctrl_sample_notification_label.config(fg = 'green')
        
        check_ctrl_sample_validity_function()

    else:
        check_ctrl_sample_validity_function()

    
def check_chip_sample_validity_function():

    chip_sample_format_error_state = 0

    if read_mode_string_var.get() == 'single':

        if bool(chip_rep1_r1_string_var.get()) and not any(chip_rep1_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep1_r1_entry.config(fg = 'red')
            chip_sample_format_error_state = 1

        if bool(chip_rep2_r1_string_var.get()) and not any(chip_rep2_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep2_r1_entry.config(fg = 'red')
            chip_sample_format_error_state = 1

        if bool(chip_rep3_r1_string_var.get()) and not any(chip_rep3_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep3_r1_entry.config(fg = 'red')
            chip_sample_format_error_state = 1

        if bool(chip_rep4_r1_string_var.get()) and not any(chip_rep4_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep4_r1_entry.config(fg = 'red')
            chip_sample_format_error_state = 1

        if bool(chip_rep5_r1_string_var.get()) and not any(chip_rep5_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep5_r1_entry.config(fg = 'red')
            chip_sample_format_error_state = 1


    if read_mode_string_var.get() == 'paired':

        if bool(chip_rep1_r1_string_var.get()) and not any(chip_rep1_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep1_r1_entry.config(fg = 'red')
            chip_sample_format_error_state = 1

        if bool(chip_rep1_r2_string_var.get()) and not any(chip_rep1_r2_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep1_r2_entry.config(fg = 'red')
            chip_sample_format_error_state = 1

        if bool(chip_rep2_r1_string_var.get()) and not any(chip_rep2_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep2_r1_entry.config(fg = 'red')
            chip_sample_format_error_state = 1

        if bool(chip_rep2_r2_string_var.get()) and not any(chip_rep2_r2_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep2_r2_entry.config(fg = 'red')
            chip_sample_format_error_state = 1

        if bool(chip_rep3_r1_string_var.get()) and not any(chip_rep3_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep3_r1_entry.config(fg = 'red')
            chip_sample_format_error_state = 1

        if bool(chip_rep3_r2_string_var.get()) and not any(chip_rep3_r2_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep3_r2_entry.config(fg = 'red')
            chip_sample_format_error_state = 1

        if bool(chip_rep4_r1_string_var.get()) and not any(chip_rep4_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep4_r1_entry.config(fg = 'red')
            chip_sample_format_error_state = 1

        if bool(chip_rep4_r2_string_var.get()) and not any(chip_rep4_r2_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep4_r2_entry.config(fg = 'red')
            chip_sample_format_error_state = 1

        if bool(chip_rep5_r1_string_var.get()) and not any(chip_rep5_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep5_r1_entry.config(fg = 'red')
            chip_sample_format_error_state = 1

        if bool(chip_rep5_r2_string_var.get()) and not any(chip_rep5_r2_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            chip_rep5_r2_entry.config(fg = 'red')
            chip_sample_format_error_state = 1


    if chip_sample_format_error_state == 1:
        chip_sample_notification_string_var.set('One or more ChIP samples do not have a valid file extension')
        chip_sample_notification_label.config(fg = 'red')
        chip_list_frame_continue_button.config(state = tk.DISABLED)
        return

    elif chip_sample_format_error_state == 0:
        chip_sample_notification_string_var.set('All ChIP samples have valid file extension')
        chip_sample_notification_label.config(fg = 'green')
    
    check_chip_sample_file_exist_function()


def check_ctrl_sample_validity_function():

    ctrl_sample_format_error_state = 0

    if read_mode_string_var.get() == 'single':

        if bool(ctrl_rep1_r1_string_var.get()) and not any(ctrl_rep1_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep1_r1_entry.config(fg = 'red')
            ctrl_sample_format_error_state = 1

        if bool(ctrl_rep2_r1_string_var.get()) and not any(ctrl_rep2_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep2_r1_entry.config(fg = 'red')
            ctrl_sample_format_error_state = 1

        if bool(ctrl_rep3_r1_string_var.get()) and not any(ctrl_rep3_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep3_r1_entry.config(fg = 'red')
            ctrl_sample_format_error_state = 1

        if bool(ctrl_rep4_r1_string_var.get()) and not any(ctrl_rep4_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep4_r1_entry.config(fg = 'red')
            ctrl_sample_format_error_state = 1

        if bool(ctrl_rep5_r1_string_var.get()) and not any(ctrl_rep5_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep5_r1_entry.config(fg = 'red')
            ctrl_sample_format_error_state = 1


    if read_mode_string_var.get() == 'paired':

        if bool(ctrl_rep1_r1_string_var.get()) and not any(ctrl_rep1_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep1_r1_entry.config(fg = 'red')
            ctrl_sample_format_error_state = 1

        if bool(ctrl_rep1_r2_string_var.get()) and not any(ctrl_rep1_r2_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep1_r2_entry.config(fg = 'red')
            ctrl_sample_format_error_state = 1

        if bool(ctrl_rep2_r1_string_var.get()) and not any(ctrl_rep2_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep2_r1_entry.config(fg = 'red')
            ctrl_sample_format_error_state = 1

        if bool(ctrl_rep2_r2_string_var.get()) and not any(ctrl_rep2_r2_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep2_r2_entry.config(fg = 'red')
            ctrl_sample_format_error_state = 1

        if bool(ctrl_rep3_r1_string_var.get()) and not any(ctrl_rep3_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep3_r1_entry.config(fg = 'red')
            ctrl_sample_format_error_state = 1

        if bool(ctrl_rep3_r2_string_var.get()) and not any(ctrl_rep3_r2_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep3_r2_entry.config(fg = 'red')
            ctrl_sample_format_error_state = 1

        if bool(ctrl_rep4_r1_string_var.get()) and not any(ctrl_rep4_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep4_r1_entry.config(fg = 'red')
            ctrl_sample_format_error_state = 1

        if bool(ctrl_rep4_r2_string_var.get()) and not any(ctrl_rep4_r2_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep4_r2_entry.config(fg = 'red')
            ctrl_sample_format_error_state = 1

        if bool(ctrl_rep5_r1_string_var.get()) and not any(ctrl_rep5_r1_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep5_r1_entry.config(fg = 'red')
            ctrl_sample_format_error_state = 1

        if bool(ctrl_rep5_r2_string_var.get()) and not any(ctrl_rep5_r2_string_var.get().endswith(valid_extension) for valid_extension in valid_extension_list):
            ctrl_rep5_r2_entry.config(fg = 'red')
            ctrl_sample_format_error_state = 1


    if ctrl_sample_format_error_state == 1:
        ctrl_sample_notification_string_var.set('One or more control samples do not have a valid file extension')
        ctrl_sample_notification_label.config(fg = 'red')
        ctrl_list_frame_continue_button.config(state = tk.DISABLED)
        return

    elif ctrl_sample_format_error_state == 0:
        ctrl_sample_notification_string_var.set('All control samples have valid file extension')
        ctrl_sample_notification_label.config(fg = 'green')
    
    check_ctrl_sample_file_exist_function()


def check_chip_sample_file_exist_function():

    chip_file_exist_error_state = 0

    if chip_rep1_r1_string_var.get():
        if not os.path.isfile(chip_rep1_r1_string_var.get()):
            chip_file_exist_error_state = 1
            chip_rep1_r1_entry.config(fg = 'red')

    if chip_rep1_r2_string_var.get():
        if not os.path.isfile(chip_rep1_r2_string_var.get()):
            chip_file_exist_error_state = 1
            chip_rep1_r2_entry.config(fg = 'red')

    if chip_rep2_r1_string_var.get():
        if not os.path.isfile(chip_rep2_r1_string_var.get()):
            chip_file_exist_error_state = 1
            chip_rep2_r1_entry.config(fg = 'red')

    if chip_rep2_r2_string_var.get():
        if not os.path.isfile(chip_rep2_r2_string_var.get()):
            chip_file_exist_error_state = 1
            chip_rep2_r2_entry.config(fg = 'red')

    if chip_rep3_r1_string_var.get():
        if not os.path.isfile(chip_rep3_r1_string_var.get()):
            chip_file_exist_error_state = 1
            chip_rep3_r1_entry.config(fg = 'red')

    if chip_rep3_r2_string_var.get():
        if not os.path.isfile(chip_rep3_r2_string_var.get()):
            chip_file_exist_error_state = 1
            chip_rep3_r2_entry.config(fg = 'red')

    if chip_rep4_r1_string_var.get():
        if not os.path.isfile(chip_rep4_r1_string_var.get()):
            chip_file_exist_error_state = 1
            chip_rep4_r1_entry.config(fg = 'red')

    if chip_rep4_r2_string_var.get():
        if not os.path.isfile(chip_rep4_r2_string_var.get()):
            chip_file_exist_error_state = 1
            chip_rep4_r2_entry.config(fg = 'red')

    if chip_rep5_r1_string_var.get():
        if not os.path.isfile(chip_rep5_r1_string_var.get()):
            chip_file_exist_error_state = 1
            chip_rep5_r1_entry.config(fg = 'red')

    if chip_rep5_r2_string_var.get():
        if not os.path.isfile(chip_rep5_r2_string_var.get()):
            chip_file_exist_error_state = 1
            chip_rep5_r2_entry.config(fg = 'red')

    if chip_file_exist_error_state == 1:
        chip_sample_notification_string_var.set('One or more ChIP sample files do not exist')
        chip_sample_notification_label.config(fg = 'red')
        chip_list_frame_continue_button.config(state = tk.DISABLED)

    elif chip_file_exist_error_state == 0:
        chip_sample_notification_string_var.set('No problem found in ChIP samples')
        chip_sample_notification_label.config(fg = 'green')
        chip_list_frame_continue_button.config(state = tk.NORMAL)


def check_ctrl_sample_file_exist_function():
    
    ctrl_file_exist_error_state = 0

    if ctrl_rep1_r1_string_var.get():
        if not os.path.isfile(ctrl_rep1_r1_string_var.get()):
            ctrl_file_exist_error_state = 1
            ctrl_rep1_r1_entry.config(fg = 'red')

    if ctrl_rep1_r2_string_var.get():
        if not os.path.isfile(ctrl_rep1_r2_string_var.get()):
            ctrl_file_exist_error_state = 1
            ctrl_rep1_r2_entry.config(fg = 'red')

    if ctrl_rep2_r1_string_var.get():
        if not os.path.isfile(ctrl_rep2_r1_string_var.get()):
            ctrl_file_exist_error_state = 1
            ctrl_rep2_r1_entry.config(fg = 'red')

    if ctrl_rep2_r2_string_var.get():
        if not os.path.isfile(ctrl_rep2_r2_string_var.get()):
            ctrl_file_exist_error_state = 1
            ctrl_rep2_r2_entry.config(fg = 'red')

    if ctrl_rep3_r1_string_var.get():
        if not os.path.isfile(ctrl_rep3_r1_string_var.get()):
            ctrl_file_exist_error_state = 1
            ctrl_rep3_r1_entry.config(fg = 'red')

    if ctrl_rep3_r2_string_var.get():
        if not os.path.isfile(ctrl_rep3_r2_string_var.get()):
            ctrl_file_exist_error_state = 1
            ctrl_rep3_r2_entry.config(fg = 'red')

    if ctrl_rep4_r1_string_var.get():
        if not os.path.isfile(ctrl_rep4_r1_string_var.get()):
            ctrl_file_exist_error_state = 1
            ctrl_rep4_r1_entry.config(fg = 'red')

    if ctrl_rep4_r2_string_var.get():
        if not os.path.isfile(ctrl_rep4_r2_string_var.get()):
            ctrl_file_exist_error_state = 1
            ctrl_rep4_r2_entry.config(fg = 'red')

    if ctrl_rep5_r1_string_var.get():
        if not os.path.isfile(ctrl_rep5_r1_string_var.get()):
            ctrl_file_exist_error_state = 1
            ctrl_rep5_r1_entry.config(fg = 'red')

    if ctrl_rep5_r2_string_var.get():
        if not os.path.isfile(ctrl_rep5_r2_string_var.get()):
            ctrl_file_exist_error_state = 1
            ctrl_rep5_r2_entry.config(fg = 'red')

    if ctrl_file_exist_error_state == 1:
        ctrl_sample_notification_string_var.set('One or more control sample files do not exist')
        ctrl_sample_notification_label.config(fg = 'red')
        ctrl_list_frame_continue_button.config(state = tk.DISABLED)

    elif ctrl_file_exist_error_state == 0:
        ctrl_sample_notification_string_var.set('No problem found in control samples')
        ctrl_sample_notification_label.config(fg = 'green')
        ctrl_list_frame_continue_button.config(state = tk.NORMAL)


def check_required_input_function(*args):
    
    if read_mode_string_var.get() != 'single' and read_mode_string_var.get() != 'paired':
        read_mode_continue_button.config(state = tk.DISABLED)
    else:
        read_mode_continue_button.config(state = tk.NORMAL)

    if peak_type_string_var.get() != 'narrow' and peak_type_string_var.get() != 'broad':
        peak_type_continue_button.config(state = tk.DISABLED)
    else:
        peak_type_continue_button.config(state = tk.NORMAL)


    if not bool(sample_table_question_string_var.get()):
        sample_table_continue_button.config(state = tk.DISABLED)
    else:
        if sample_table_question_string_var.get() == 'yes':
            if not bool(sample_table_string_var.get()):
                sample_table_continue_button.config(state = tk.DISABLED)
            else:
                sample_table_continue_button.config(state = tk.NORMAL)
        if sample_table_question_string_var.get() == 'no':    
            sample_table_continue_button.config(state = tk.NORMAL)


    if not bool(setting_table_question_string_var.get()):
        setting_table_continue_button.config(state = tk.DISABLED)
    else:
        if setting_table_question_string_var.get() == 'yes':
            if not bool(setting_table_string_var.get()):
                setting_table_continue_button.config(state = tk.DISABLED)
            else:
                setting_table_continue_button.config(state = tk.NORMAL)
        if setting_table_question_string_var.get() == 'no':    
            setting_table_continue_button.config(state = tk.NORMAL)


    if not bool(genome_ref_string_var.get()):
        genome_ref_continue_button.config(state = tk.DISABLED)
    elif not bool(genome_folder_string_var.get()):
        genome_ref_continue_button.config(state = tk.DISABLED)
    else:
        genome_ref_continue_button.config(state = tk.NORMAL)


    if bool(known_motif_string_var.get()):
        if not os.path.isfile(known_motif_string_var.get()):
            known_motif_continue_button.config(state = tk.DISABLED)
        else:
            known_motif_continue_button.config(state = tk.NORMAL)
    else:
        known_motif_continue_button.config(state = tk.NORMAL)


    if not bool(setname_string_var.get()):
        output_folder_continue_button.config(state = tk.DISABLED)
    elif not bool(output_folder_string_var.get()):
        output_folder_continue_button.config(state = tk.DISABLED)       
    else:
        output_folder_continue_button.config(state = tk.NORMAL)


    if not bool(cpu_count_string_var.get()):
        generate_scripts_button.config(state = tk.DISABLED)
        generate_and_run_scripts_button.config(state = tk.DISABLED)
        cpu_count_notification_string_var.set("Maximum number of CPU cores available: {}".format(max_cpu))
        cpu_count_notification_label.config(fg = 'blue')
    
    else:
        if int(cpu_count_string_var.get()) > max_cpu:
            generate_scripts_button.config(state = tk.DISABLED)
            generate_and_run_scripts_button.config(state = tk.DISABLED)
            cpu_count_notification_string_var.set("Entered number exceeds available CPU cores ({})".format(max_cpu))
            cpu_count_notification_label.config(fg = 'red')

        elif int(cpu_count_string_var.get()) < 1:
            generate_scripts_button.config(state = tk.DISABLED)
            generate_and_run_scripts_button.config(state = tk.DISABLED)
            cpu_count_notification_string_var.set("Need at least one CPU core to run the pipeline")
            cpu_count_notification_label.config(fg = 'red')

        elif int(cpu_count_string_var.get()) <= max_cpu:
            generate_scripts_button.config(state = tk.NORMAL)
            generate_and_run_scripts_button.config(state = tk.NORMAL)
            cpu_count_notification_string_var.set("ChIP-AP ready for action!")
            cpu_count_notification_label.config(fg = 'green')
    

def update_command_line_function(*args):

    if read_mode_string_var.get() == 'single' or read_mode_string_var.get() == 'paired':
        read_mode_arg.set(' --mode {}'.format(read_mode_string_var.get()))
    else:
        read_mode_arg.set('')

    if peak_type_string_var.get() == 'narrow' or peak_type_string_var.get() == 'broad':
        peak_type_arg.set(' --peak {}'.format(peak_type_string_var.get()))
    else:
        peak_type_arg.set('')

    if bool(output_folder_string_var.get()):
        output_folder_arg.set(' --output {}'.format(output_folder_string_var.get()))
    else:
        output_folder_arg.set('')

    if bool(setname_string_var.get()):
        setname_arg.set(' --setname {}'.format(setname_string_var.get()))
    else:
        setname_arg.set('')

    output_dir.set('{}/{}'.format(os.path.abspath(output_folder_string_var.get()), setname_string_var.get()))

    if bool(genome_ref_string_var.get()):
        genome_ref_arg.set(' --ref {}'.format(genome_ref_string_var.get().split(' ')[0]))
    else:
        genome_ref_arg.set('')

    if bool(genome_folder_string_var.get()):
        genome_folder_arg.set(' --genome {}'.format(genome_folder_string_var.get()))
    else:
        genome_folder_arg.set('')

    if bool(output_folder_string_var.get()) and bool(setname_string_var.get()):
        sample_table_arg.set(' --sample_table {}/{}_sample_table.tsv'.format(output_dir.get(), setname_string_var.get()))
    else:
        sample_table_arg.set('')

    if bool(output_folder_string_var.get()) and bool(setname_string_var.get()):
        setting_table_arg.set(' --custom_setting_table {}/{}_setting_table.tsv'.format(output_dir.get(), setname_string_var.get()))
    else:
        setting_table_arg.set('')

    if bool(known_motif_string_var.get()):
        known_motif_arg.set(' --motif {}'.format(known_motif_string_var.get()))
    else:
        known_motif_arg.set('')

    if fcmerge_var.get() == 1:
        fcmerge_arg.set(' --fcmerge')
    else:
        fcmerge_arg.set('')

    if goann_var.get() == 1:
        goann_arg.set(' --goann')
    else:
        goann_arg.set('')

    if pathann_var.get() == 1:
        pathann_arg.set(' --pathann')
    else:
        pathann_arg.set('')

    if deltemp_var.get() == 1:
        deltemp_arg.set(' --deltemp')
    else:
        deltemp_arg.set('')

    if bool(cpu_count_string_var.get()):
        cpu_count_arg.set(' --thread {}'.format(cpu_count_string_var.get()))
    else:
        cpu_count_arg.set('')

    if stdout_var.get() == 1 and bool(output_folder_string_var.get()) and bool(setname_string_var.get()):
        stdout_arg.set(' 1> {}/{}.out'.format(output_dir.get(), setname_string_var.get()))
    else:
        stdout_arg.set('')

    if stderr_var.get() == 1 and bool(output_folder_string_var.get()) and bool(setname_string_var.get()):
        stderr_arg.set(' 2> {}/{}.err'.format(output_dir.get(), setname_string_var.get()))
    else:
        stderr_arg.set('')

    command_line_output_string_var.set('{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}'.format(current_pipeline,
                                                                                read_mode_arg.get(),
                                                                                peak_type_arg.get(),
                                                                                output_folder_arg.get(),
                                                                                setname_arg.get(),
                                                                                genome_ref_arg.get(),
                                                                                genome_folder_arg.get(),
                                                                                sample_table_arg.get(),
                                                                                setting_table_arg.get(),
                                                                                known_motif_arg.get(),
                                                                                fcmerge_arg.get(),
                                                                                goann_arg.get(),
                                                                                pathann_arg.get(),
                                                                                deltemp_arg.get(),
                                                                                cpu_count_arg.get(),
                                                                                stdout_arg.get(),
                                                                                stderr_arg.get()))


def check_traced_input_function(*args):
    check_required_input_function()
    update_command_line_function()


########################################################################################################################


max_cpu = multiprocessing.cpu_count()

current_dir_string_var = tk.StringVar(value = default_current_dir)

read_mode_string_var = tk.StringVar()
peak_type_string_var = tk.StringVar()

chip_rep1_r1_string_var = tk.StringVar()
chip_rep1_r2_string_var = tk.StringVar()
chip_rep2_r1_string_var = tk.StringVar()
chip_rep2_r2_string_var = tk.StringVar()
chip_rep3_r1_string_var = tk.StringVar()
chip_rep3_r2_string_var = tk.StringVar()
chip_rep4_r1_string_var = tk.StringVar()
chip_rep4_r2_string_var = tk.StringVar()
chip_rep5_r1_string_var = tk.StringVar()
chip_rep5_r2_string_var = tk.StringVar()

ctrl_rep1_r1_string_var = tk.StringVar()
ctrl_rep1_r2_string_var = tk.StringVar()
ctrl_rep2_r1_string_var = tk.StringVar()
ctrl_rep2_r2_string_var = tk.StringVar()
ctrl_rep3_r1_string_var = tk.StringVar()
ctrl_rep3_r2_string_var = tk.StringVar()
ctrl_rep4_r1_string_var = tk.StringVar()
ctrl_rep4_r2_string_var = tk.StringVar()
ctrl_rep5_r1_string_var = tk.StringVar()
ctrl_rep5_r2_string_var = tk.StringVar()

sample_table_question_string_var = tk.StringVar()
sample_table_string_var = tk.StringVar()
sample_table_notification_string_var = tk.StringVar()
chip_sample_notification_string_var = tk.StringVar()
ctrl_sample_notification_string_var = tk.StringVar()
chip_rep_number_int_var = tk.IntVar(value = 1)
ctrl_rep_number_int_var = tk.IntVar(value = 1)

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

setting_table_question_string_var = tk.StringVar()
setting_table_string_var = tk.StringVar()
setting_table_notification_string_var = tk.StringVar()

genome_ref_string_var = tk.StringVar(value = genome_ref_options[0])
genome_folder_string_var = tk.StringVar(value = genome_folder_full_path)
known_motif_string_var = tk.StringVar()
setname_string_var = tk.StringVar()
output_folder_string_var = tk.StringVar()

fcmerge_var = tk.IntVar()
goann_var = tk.IntVar()
pathann_var = tk.IntVar()
deltemp_var = tk.IntVar(value = 1)
stdout_var = tk.IntVar()
stderr_var = tk.IntVar()
cpu_count_string_var = tk.StringVar()
cpu_count_notification_string_var = tk.StringVar()

command_line_output_string_var = tk.StringVar()

read_mode_arg = tk.StringVar()
peak_type_arg = tk.StringVar()
output_folder_arg = tk.StringVar()
setname_arg = tk.StringVar()
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

output_dir = tk.StringVar()

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
setting_table_question_string_var.trace('w', check_required_input_function)
setting_table_string_var.trace('w', check_traced_input_function)
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
stderr_var.trace('w', update_command_line_function)
cpu_count_string_var.trace('w', check_traced_input_function)


########################################################################################################################


opening_frame = tk.Frame(root)

opening_label = tk.Label(opening_frame, text = "Welcome to ChIP-AP")
opening_label.grid(row = 1, column = 1, padx = 10, pady = (5,0), columnspan = 2)

opening_frame_exit_button = tk.Button(opening_frame, text = "Exit wizard", command = lambda : exit(), width = 20)
opening_frame_continue_button = tk.Button(opening_frame, text = "Continue >>", command = lambda : change_frame_function(opening_frame, read_mode_frame), width = 20)

opening_frame_exit_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
opening_frame_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))

opening_frame.grid(row = 0, column = 0, sticky = "news")


########################################################################################################################


read_mode_frame = tk.Frame(root)

read_mode_label = tk.Label(read_mode_frame, text = "Dataset sequencing mode:", justify = tk.LEFT, width = 30)
read_mode_label.grid(row = 1, column = 1, padx = 10, pady = 5, columnspan = 2)

single_end_radio = tk.Radiobutton(read_mode_frame, text = "Single end", padx = 5, variable = read_mode_string_var, value = 'single', width = 30)
CreateToolTip(single_end_radio, text = 'Select this if there is one sequencer output file per sample.')
single_end_radio.grid(row = 2, column = 1, padx = 10, pady = 2, columnspan = 2)

paired_end_radio = tk.Radiobutton(read_mode_frame, text = "Paired end", padx = 5, variable = read_mode_string_var, value = 'paired', width = 30)
CreateToolTip(paired_end_radio, text = 'Select this if there are two sequencer output files per sample.\nThey are typically in pairs R1 and R2 for every sample.')
paired_end_radio.grid(row = 3, column = 1, padx = 10, pady = 2, columnspan = 2)

read_mode_back_button = tk.Button(read_mode_frame, text = "<< Back", command = lambda : change_frame_function(read_mode_frame, opening_frame), width = 20)
read_mode_continue_button = tk.Button(read_mode_frame, text = "Continue >>", command = lambda : change_frame_function(read_mode_frame, peak_type_frame), width = 20)

read_mode_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
read_mode_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################


peak_type_frame = tk.Frame(root)

peak_type_label = tk.Label(peak_type_frame, text = "Dataset peak type:", justify = tk.LEFT, width = 30)
peak_type_label.grid(row = 1, column = 1, padx = 10, pady = 5, columnspan = 2)

narrow_peak_radio = tk.Radiobutton(peak_type_frame, text = "Narrow peaks", padx = 5, variable = peak_type_string_var, value = 'narrow', width = 30)
CreateToolTip(narrow_peak_radio, text = 'Select this for ChIP-seq experiment with transcription factor protein.')
narrow_peak_radio.grid(row = 2, column = 1, padx = 10, pady = 2, columnspan = 2)

broad_peak_radio = tk.Radiobutton(peak_type_frame, text = "Broad peaks", padx = 5, variable = peak_type_string_var, value = 'broad', width = 30)
CreateToolTip(broad_peak_radio, text = 'Select this for ChIP-seq experiment with chromatin modifier protein.')
broad_peak_radio.grid(row = 3, column = 1, padx = 10, pady = 2, columnspan = 2)

peak_type_back_button = tk.Button(peak_type_frame, text = "<< Back", command = lambda : change_frame_function(peak_type_frame, read_mode_frame), width = 20)
peak_type_continue_button = tk.Button(peak_type_frame, text = "Continue >>", command = lambda : change_frame_function(peak_type_frame, sample_table_frame), width = 20)

peak_type_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
peak_type_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################


def sample_table_button_function():
    sample_table_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = sample_table_frame, mode = 'rb', title = 'Choose a sample table file')
    if sample_table_button_input:
        sample_table_string_var.set(sample_table_button_input.name)
        current_dir_string_var.set(os.path.dirname(sample_table_button_input.name))
        sample_table_loading_test_function()

    else:
        clear_sample_function()
        sample_table_notification_string_var.set('No sample table was selected. No samples were loaded.')
        sample_table_notification_label.config(fg = 'red')
        
        rearrange_entry_field_function()
        return

        
def sample_table_loading_test_function():
    try:
        sample_table_absolute_path = sample_table_string_var.get()
        sample_table_df = pd.read_csv(sample_table_absolute_path, delimiter='\t')
        sample_table_df.fillna('', inplace = True)

    except:
        chip_r1_sample_list = []
        ctrl_r1_sample_list = []
        chip_r2_sample_list = []
        ctrl_r2_sample_list = []
        clear_sample_function()
        sample_table_notification_string_var.set('Sample table loading error! No samples were loaded.')
        sample_table_notification_label.config(fg = 'red')
        return False

    try:
        chip_r1_sample_list = [chip_r1_sample if str(chip_r1_sample) != 'nan' else '' for chip_r1_sample in sample_table_df['chip_read_1']]
        ctrl_r1_sample_list = [ctrl_r1_sample if str(ctrl_r1_sample) != 'nan' else '' for ctrl_r1_sample in sample_table_df['ctrl_read_1']]
        chip_r2_sample_list = [chip_r2_sample if str(chip_r2_sample) != 'nan' else '' for chip_r2_sample in sample_table_df['chip_read_2']]
        ctrl_r2_sample_list = [ctrl_r2_sample if str(ctrl_r2_sample) != 'nan' else '' for ctrl_r2_sample in sample_table_df['ctrl_read_2']]

    except:
        chip_r1_sample_list = []
        ctrl_r1_sample_list = []
        chip_r2_sample_list = []
        ctrl_r2_sample_list = []
        clear_sample_function()
        sample_table_notification_string_var.set('Sample table reading error! No samples were loaded.')
        sample_table_notification_label.config(fg = 'red')
        return False


    if read_mode_string_var.get() == 'single':

        if not all(chip_r2 == '' for chip_r2 in chip_r2_sample_list) or not all(ctrl_r2 == '' for ctrl_r2 in ctrl_r2_sample_list):
            clear_sample_function()
            sample_table_notification_string_var.set('Error! Paired samples detected. No samples were loaded.')
            sample_table_notification_label.config(fg = 'red')
            return False
        
        else:
            try:
                chip_rep1_r1_string_var.set(chip_r1_sample_list[0])
                chip_rep2_r1_string_var.set(chip_r1_sample_list[1])
                chip_rep3_r1_string_var.set(chip_r1_sample_list[2])
                chip_rep4_r1_string_var.set(chip_r1_sample_list[3])
                chip_rep5_r1_string_var.set(chip_r1_sample_list[4])

            except:
                pass

            try:
                ctrl_rep1_r1_string_var.set(ctrl_r1_sample_list[0])
                ctrl_rep2_r1_string_var.set(ctrl_r1_sample_list[1])
                ctrl_rep3_r1_string_var.set(ctrl_r1_sample_list[2])
                ctrl_rep4_r1_string_var.set(ctrl_r1_sample_list[3])
                ctrl_rep5_r1_string_var.set(ctrl_r1_sample_list[4])

            except:
                pass
            
            sample_table_notification_string_var.set('Sample table loading successful.')
            sample_table_notification_label.config(fg = 'green')
            rearrange_entry_field_function()
            return True

    if read_mode_string_var.get() == 'paired':
        try:
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

        except:
            pass

        try:
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

        except:
            pass

        sample_table_notification_string_var.set('Sample table loading successful.')
        sample_table_notification_label.config(fg = 'green')
        
        rearrange_entry_field_function()
        return True


def sample_table_entry_popup_function(sample_table_entry_popup_arg):
    if sample_table_entry_popup_arg == 'yes':
        sample_table_entry.grid(row = 3, column = 1, padx = (10,5), ipady = 3)
        sample_table_button.grid(row = 3, column = 2, padx = (5,10), pady = 2)
        chip_rep_number_label.grid_remove()
        chip_rep_number_drop_down.grid_remove()
        ctrl_rep_number_label.grid_remove()
        ctrl_rep_number_drop_down.grid_remove()
        clear_sample_function()
        sample_table_notification_string_var.set('Please load your sample table file')
        sample_table_notification_label.config(fg = 'blue')
        root.eval('tk::PlaceWindow . center')
        root.grab_set()
        root.focus_force()

    if sample_table_entry_popup_arg == 'no':
        sample_table_entry.grid_remove()
        sample_table_button.grid_remove()
        chip_rep_number_label.grid(row = 5, column = 1, padx = (10,5), sticky = "e")
        chip_rep_number_drop_down.grid(row = 5, column = 2, padx = (5,10), sticky = "w")
        ctrl_rep_number_label.grid(row = 6, column = 1, padx = (10,5), sticky = "e")
        ctrl_rep_number_drop_down.grid(row = 6, column = 2, padx = (5,10), sticky = "w")
        clear_sample_function()
        sample_table_notification_string_var.set('Please choose the number of sample replicates')
        sample_table_notification_label.config(fg = 'blue')
        root.eval('tk::PlaceWindow . center')
        root.grab_set()
        root.focus_force()

    rearrange_entry_field_function()


def clear_sample_function():

    sample_table_string_var.set('')

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

    rearrange_entry_field_function()


sample_table_frame = tk.Frame(root)

sample_table_label = tk.Label(sample_table_frame, text = "Do you want to use your sample table?", width = 70)
sample_table_label.grid(row = 1, column = 1, padx = 10, pady = 5, columnspan = 2)

sample_table_yes_radio = tk.Radiobutton(sample_table_frame, text = "Yes, I would like use my sample table", padx = 5, variable = sample_table_question_string_var, value = 'yes', command = lambda : sample_table_entry_popup_function('yes'), width = 50)
CreateToolTip(sample_table_yes_radio, text = 'Select this to load a pre-populated sample table in ChIP-AP format (see GitHub documentation).\nThe ChIP and control samples list will be filled automatically based on the loaded sample table.')
sample_table_yes_radio.grid(row = 2, column = 1, padx = 10, pady = 2, columnspan = 2)

sample_table_entry = tk.Entry(sample_table_frame, textvariable = sample_table_string_var, width = 50, justify = tk.RIGHT)
sample_table_entry.xview_moveto(1)
sample_table_button = tk.Button(sample_table_frame, text = 'Browse', width = 20, command = lambda : sample_table_button_function())
CreateToolTip(sample_table_button, text = 'Click here to browse and select your sample table file.')

sample_table_no_radio = tk.Radiobutton(sample_table_frame, text = "No, I would like to assign my samples manually", padx = 5, variable = sample_table_question_string_var, value = 'no', command = lambda : sample_table_entry_popup_function('no'), width = 50)
CreateToolTip(sample_table_no_radio, text = 'Select this if you want to type in, or browse and select all your sample files yourself.')
sample_table_no_radio.grid(row = 4, column = 1, padx = 10, pady = 2, columnspan = 2)

replicate_number_choice = [1, 2, 3, 4, 5]

chip_rep_number_label = tk.Label(sample_table_frame, text = "Number of ChIP replicate(s):")

chip_rep_number_drop_down = tk.OptionMenu(sample_table_frame, chip_rep_number_int_var, *replicate_number_choice)
chip_rep_number_drop_down.config(width = 3, takefocus = 1)

chip_rep_number_drop_down_menu = root.nametowidget(chip_rep_number_drop_down.menuname)
chip_rep_number_drop_down_menu.config(font = default_font)

ctrl_rep_number_label = tk.Label(sample_table_frame, text = "Number of control replicate(s):")

ctrl_rep_number_drop_down = tk.OptionMenu(sample_table_frame, ctrl_rep_number_int_var, *replicate_number_choice)
ctrl_rep_number_drop_down.config(width = 3, takefocus = 1)

ctrl_rep_number_drop_down_menu = root.nametowidget(ctrl_rep_number_drop_down.menuname)
ctrl_rep_number_drop_down_menu.config(font = default_font)

sample_table_notification_label = tk.Label(sample_table_frame, textvariable = sample_table_notification_string_var, width = 70, padx = 5, pady = 5)
sample_table_notification_label.grid(row = 7, column = 1, padx = 10, pady = 5, columnspan = 2)

sample_table_back_button = tk.Button(sample_table_frame, text = "<< Back", command = lambda : change_frame_function(sample_table_frame, peak_type_frame), width = 20)
sample_table_continue_button = tk.Button(sample_table_frame, text = "Continue >>", command = lambda : change_frame_function(sample_table_frame, chip_list_frame), width = 20)

sample_table_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
sample_table_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################


def chip_rep1_r1_button_function():
    chip_rep1_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = chip_list_frame, mode = 'rb', title = 'Choose a file')
    if chip_rep1_r1_button_input:
        chip_rep1_r1_string_var.set(chip_rep1_r1_button_input.name)
        current_dir_string_var.set(os.path.dirname(chip_rep1_r1_button_input.name))
    rearrange_entry_field_function()

def chip_rep1_r2_button_function():
    chip_rep1_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = chip_list_frame, mode = 'rb', title = 'Choose a file')
    if chip_rep1_r2_button_input:
        chip_rep1_r2_string_var.set(chip_rep1_r2_button_input.name)
        current_dir_string_var.set(os.path.dirname(chip_rep1_r2_button_input.name))
    rearrange_entry_field_function()

def chip_rep2_r1_button_function():
    chip_rep2_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = chip_list_frame, mode = 'rb', title = 'Choose a file')
    if chip_rep2_r1_button_input:
        chip_rep2_r1_string_var.set(chip_rep2_r1_button_input.name)
        current_dir_string_var.set(os.path.dirname(chip_rep2_r1_button_input.name))
    rearrange_entry_field_function()

def chip_rep2_r2_button_function():
    chip_rep2_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = chip_list_frame, mode = 'rb', title = 'Choose a file')
    if chip_rep2_r2_button_input:
        chip_rep2_r2_string_var.set(chip_rep2_r2_button_input.name)
        current_dir_string_var.set(os.path.dirname(chip_rep2_r2_button_input.name))
    rearrange_entry_field_function()

def chip_rep3_r1_button_function():
    chip_rep3_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = chip_list_frame, mode = 'rb', title = 'Choose a file')
    if chip_rep3_r1_button_input:
        chip_rep3_r1_string_var.set(chip_rep3_r1_button_input.name)
        current_dir_string_var.set(os.path.dirname(chip_rep3_r1_button_input.name))
    rearrange_entry_field_function()

def chip_rep3_r2_button_function():
    chip_rep3_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = chip_list_frame, mode = 'rb', title = 'Choose a file')
    if chip_rep3_r2_button_input:
        chip_rep3_r2_string_var.set(chip_rep3_r2_button_input.name)
        current_dir_string_var.set(os.path.dirname(chip_rep3_r2_button_input.name))
    rearrange_entry_field_function()

def chip_rep4_r1_button_function():
    chip_rep4_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = chip_list_frame, mode = 'rb', title = 'Choose a file')
    if chip_rep4_r1_button_input:
        chip_rep4_r1_string_var.set(chip_rep4_r1_button_input.name)
        current_dir_string_var.set(os.path.dirname(chip_rep4_r1_button_input.name))
    rearrange_entry_field_function()

def chip_rep4_r2_button_function():
    chip_rep4_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = chip_list_frame, mode = 'rb', title = 'Choose a file')
    if chip_rep4_r2_button_input:
        chip_rep4_r2_string_var.set(chip_rep4_r2_button_input.name)
        current_dir_string_var.set(os.path.dirname(chip_rep4_r2_button_input.name))
    rearrange_entry_field_function()

def chip_rep5_r1_button_function():
    chip_rep5_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = chip_list_frame, mode = 'rb', title = 'Choose a file')
    if chip_rep5_r1_button_input:
        chip_rep5_r1_string_var.set(chip_rep5_r1_button_input.name)
        current_dir_string_var.set(os.path.dirname(chip_rep5_r1_button_input.name))
    rearrange_entry_field_function()

def chip_rep5_r2_button_function():
    chip_rep5_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = chip_list_frame, mode = 'rb', title = 'Choose a file')
    if chip_rep5_r2_button_input:
        chip_rep5_r2_string_var.set(chip_rep5_r2_button_input.name)
        current_dir_string_var.set(os.path.dirname(chip_rep5_r2_button_input.name))
    rearrange_entry_field_function()


def display_chip_widget_function():
    for chip_widget in chip_widget_list:
        chip_widget.grid_remove()

    if sample_table_question_string_var.get() == 'yes':
        if bool(chip_rep1_r1_string_var.get()):
            chip_rep1_r1_button.grid(column = 1, row = 1, padx = (10,5))
            chip_rep1_r1_entry.grid(column = 2, row = 1, padx = (5,10), ipady = 3)
            
            if read_mode_string_var.get() == 'paired':
                chip_rep1_r2_button.grid(column = 1, row = 2, padx = (10,5))
                chip_rep1_r2_entry.grid(column = 2, row = 2, padx = (5,10), ipady = 3)

        if bool(chip_rep2_r1_string_var.get()):
            chip_rep2_r1_button.grid(column = 1, row = 3, padx = (10,5), pady = (10,0))
            chip_rep2_r1_entry.grid(column = 2, row = 3, padx = (5,10), ipady = 3, pady = (10,0))
            
            if read_mode_string_var.get() == 'paired':
                chip_rep2_r2_button.grid(column = 1, row = 4, padx = (10,5))
                chip_rep2_r2_entry.grid(column = 2, row = 4, padx = (5,10), ipady = 3)
    
        if bool(chip_rep3_r1_string_var.get()):
            chip_rep3_r1_button.grid(column = 1, row = 5, padx = (10,5), pady = (10,0))
            chip_rep3_r1_entry.grid(column = 2, row = 5, padx = (5,10), ipady = 3, pady = (10,0))     
            
            if read_mode_string_var.get() == 'paired':
                chip_rep3_r2_button.grid(column = 1, row = 6, padx = (10,5))
                chip_rep3_r2_entry.grid(column = 2, row = 6, padx = (5,10), ipady = 3)

        if bool(chip_rep4_r1_string_var.get()):
            chip_rep4_r1_button.grid(column = 1, row = 7, padx = (10,5), pady = (10,0))
            chip_rep4_r1_entry.grid(column = 2, row = 7, padx = (5,10), ipady = 3, pady = (10,0))
            
            if read_mode_string_var.get() == 'paired':
                chip_rep4_r2_button.grid(column = 1, row = 8, padx = (10,5))
                chip_rep4_r2_entry.grid(column = 2, row = 8, padx = (5,10), ipady = 3)

        if bool(chip_rep5_r1_string_var.get()):
            chip_rep5_r1_button.grid(column = 1, row = 9, padx = (10,5), pady = (10,0))
            chip_rep5_r1_entry.grid(column = 2, row = 9, padx = (5,10), ipady = 3, pady = (10,0))
            
            if read_mode_string_var.get() == 'paired':
                chip_rep5_r2_button.grid(column = 1, row = 10, padx = (10,5))
                chip_rep5_r2_entry.grid(column = 2, row = 10, padx = (5,10), ipady = 3)

    if sample_table_question_string_var.get() == 'no':
        if chip_rep_number_int_var.get() >= 1:
            chip_rep1_r1_button.grid(column = 1, row = 1, padx = (10,5))
            chip_rep1_r1_entry.grid(column = 2, row = 1, padx = (5,10), ipady = 3)
            
            if read_mode_string_var.get() == 'paired':
                chip_rep1_r2_button.grid(column = 1, row = 2, padx = (10,5))
                chip_rep1_r2_entry.grid(column = 2, row = 2, padx = (5,10), ipady = 3)

        if chip_rep_number_int_var.get() >= 2:
            chip_rep2_r1_button.grid(column = 1, row = 3, padx = (10,5), pady = (10,0))
            chip_rep2_r1_entry.grid(column = 2, row = 3, padx = (5,10), ipady = 3, pady = (10,0))
            
            if read_mode_string_var.get() == 'paired':
                chip_rep2_r2_button.grid(column = 1, row = 4, padx = (10,5))
                chip_rep2_r2_entry.grid(column = 2, row = 4, padx = (5,10), ipady = 3)
    
        if chip_rep_number_int_var.get() >= 3:
            chip_rep3_r1_button.grid(column = 1, row = 5, padx = (10,5), pady = (10,0))
            chip_rep3_r1_entry.grid(column = 2, row = 5, padx = (5,10), ipady = 3, pady = (10,0))
            
            if read_mode_string_var.get() == 'paired':
                chip_rep3_r2_button.grid(column = 1, row = 6, padx = (10,5))
                chip_rep3_r2_entry.grid(column = 2, row = 6, padx = (5,10), ipady = 3)

        if chip_rep_number_int_var.get() >= 4:
            chip_rep4_r1_button.grid(column = 1, row = 7, padx = (10,5), pady = (10,0))
            chip_rep4_r1_entry.grid(column = 2, row = 7, padx = (5,10), ipady = 3, pady = (10,0))

            if read_mode_string_var.get() == 'paired':
                chip_rep4_r2_button.grid(column = 1, row = 8, padx = (10,5))
                chip_rep4_r2_entry.grid(column = 2, row = 8, padx = (5,10), ipady = 3)

        if chip_rep_number_int_var.get() >= 5:
            chip_rep5_r1_button.grid(column = 1, row = 9, padx = (10,5), pady = (10,0))
            chip_rep5_r1_entry.grid(column = 2, row = 9, padx = (5,10), ipady = 3, pady = (10,0))

            if read_mode_string_var.get() == 'paired':
                chip_rep5_r2_button.grid(column = 1, row = 10, padx = (10,5))
                chip_rep5_r2_entry.grid(column = 2, row = 10, padx = (5,10), ipady = 3)

    rearrange_entry_field_function()


chip_list_frame = tk.Frame(root)

chip_list_label = tk.Label(chip_list_frame, text = "Please type in, or browse and select your ChIP sample files below", justify = tk.LEFT, width = 70)
chip_list_label.grid(row = 0, column = 1, padx = 10, pady = (5,10), columnspan = 2)

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

CreateToolTip(chip_rep1_r1_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam file extension accepted)')
CreateToolTip(chip_rep1_r2_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam file extension accepted)')
CreateToolTip(chip_rep2_r1_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam file extension accepted)')
CreateToolTip(chip_rep2_r2_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam file extension accepted)')
CreateToolTip(chip_rep3_r1_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam file extension accepted)')
CreateToolTip(chip_rep3_r2_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam file extension accepted)')
CreateToolTip(chip_rep4_r1_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam file extension accepted)')
CreateToolTip(chip_rep4_r2_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam file extension accepted)')
CreateToolTip(chip_rep5_r1_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam file extension accepted)')
CreateToolTip(chip_rep5_r2_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam file extension accepted)')

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

chip_sample_notification_label = tk.Label(chip_list_frame, textvariable = chip_sample_notification_string_var, width = 70, padx = 5, pady = 5)
chip_sample_notification_label.grid(row = 15, column = 1, padx = 10, pady = 5, columnspan = 2)

chip_list_frame_back_button = tk.Button(chip_list_frame, text = "<< Back", command = lambda : change_frame_function(chip_list_frame, sample_table_frame), width = 20)
chip_list_frame_continue_button = tk.Button(chip_list_frame, text = "Continue >>", command = lambda : change_frame_function(chip_list_frame, ctrl_list_frame), width = 20)

chip_list_frame_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
chip_list_frame_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################


def ctrl_rep1_r1_button_function():
    ctrl_rep1_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = ctrl_list_frame, mode = 'rb', title = 'Choose a file')
    if ctrl_rep1_r1_button_input:
        ctrl_rep1_r1_string_var.set(ctrl_rep1_r1_button_input.name)
        current_dir_string_var.set(os.path.dirname(ctrl_rep1_r1_button_input.name))
    rearrange_entry_field_function()

def ctrl_rep1_r2_button_function():
    ctrl_rep1_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = ctrl_list_frame, mode = 'rb', title = 'Choose a file')
    if ctrl_rep1_r2_button_input:
        ctrl_rep1_r2_string_var.set(ctrl_rep1_r2_button_input.name)
        current_dir_string_var.set(os.path.dirname(ctrl_rep1_r2_button_input.name))
    rearrange_entry_field_function()

def ctrl_rep2_r1_button_function():
    ctrl_rep2_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = ctrl_list_frame, mode = 'rb', title = 'Choose a file')
    if ctrl_rep2_r1_button_input:
        ctrl_rep2_r1_string_var.set(ctrl_rep2_r1_button_input.name)
        current_dir_string_var.set(os.path.dirname(ctrl_rep2_r1_button_input.name))
    rearrange_entry_field_function()

def ctrl_rep2_r2_button_function():
    ctrl_rep2_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = ctrl_list_frame, mode = 'rb', title = 'Choose a file')
    if ctrl_rep2_r2_button_input:
        ctrl_rep2_r2_string_var.set(ctrl_rep2_r2_button_input.name)
        current_dir_string_var.set(os.path.dirname(ctrl_rep2_r2_button_input.name))
    rearrange_entry_field_function()

def ctrl_rep3_r1_button_function():
    ctrl_rep3_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = ctrl_list_frame, mode = 'rb', title = 'Choose a file')
    if ctrl_rep3_r1_button_input:
        ctrl_rep3_r1_string_var.set(ctrl_rep3_r1_button_input.name)
        current_dir_string_var.set(os.path.dirname(ctrl_rep3_r1_button_input.name))
    rearrange_entry_field_function()

def ctrl_rep3_r2_button_function():
    ctrl_rep3_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = ctrl_list_frame, mode = 'rb', title = 'Choose a file')
    if ctrl_rep3_r2_button_input:
        ctrl_rep3_r2_string_var.set(ctrl_rep3_r2_button_input.name)
        current_dir_string_var.set(os.path.dirname(ctrl_rep3_r2_button_input.name))
    rearrange_entry_field_function()

def ctrl_rep4_r1_button_function():
    ctrl_rep4_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = ctrl_list_frame, mode = 'rb', title = 'Choose a file')
    if ctrl_rep4_r1_button_input:
        ctrl_rep4_r1_string_var.set(ctrl_rep4_r1_button_input.name)
        current_dir_string_var.set(os.path.dirname(ctrl_rep4_r1_button_input.name))
    rearrange_entry_field_function()

def ctrl_rep4_r2_button_function():
    ctrl_rep4_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = ctrl_list_frame, mode = 'rb', title = 'Choose a file')
    if ctrl_rep4_r2_button_input:
        ctrl_rep4_r2_string_var.set(ctrl_rep4_r2_button_input.name)
        current_dir_string_var.set(os.path.dirname(ctrl_rep4_r2_button_input.name))
    rearrange_entry_field_function()

def ctrl_rep5_r1_button_function():
    ctrl_rep5_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = ctrl_list_frame, mode = 'rb', title = 'Choose a file')
    if ctrl_rep5_r1_button_input:
        ctrl_rep5_r1_string_var.set(ctrl_rep5_r1_button_input.name)
        current_dir_string_var.set(os.path.dirname(ctrl_rep5_r1_button_input.name))
    rearrange_entry_field_function()

def ctrl_rep5_r2_button_function():
    ctrl_rep5_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = ctrl_list_frame, mode = 'rb', title = 'Choose a file')
    if ctrl_rep5_r2_button_input:
        ctrl_rep5_r2_string_var.set(ctrl_rep5_r2_button_input.name)
        current_dir_string_var.set(os.path.dirname(ctrl_rep5_r2_button_input.name))
    rearrange_entry_field_function()


def display_ctrl_widget_function():
    for ctrl_widget in ctrl_widget_list:
        ctrl_widget.grid_remove()

    if sample_table_question_string_var.get() == 'yes':
        if bool(ctrl_rep1_r1_string_var.get()):
            ctrl_rep1_r1_button.grid(column = 1, row = 1, padx = (10,5))
            ctrl_rep1_r1_entry.grid(column = 2, row = 1, padx = (5,10), ipady = 3)
        
            if read_mode_string_var.get() == 'paired':
                ctrl_rep1_r2_button.grid(column = 1, row = 2, padx = (10,5))
                ctrl_rep1_r2_entry.grid(column = 2, row = 2, padx = (5,10), ipady = 3)

        if bool(ctrl_rep2_r1_string_var.get()):
            ctrl_rep2_r1_button.grid(column = 1, row = 3, padx = (10,5), pady = (10,0))
            ctrl_rep2_r1_entry.grid(column = 2, row = 3, padx = (5,10), ipady = 3, pady = (10,0))
            
            if read_mode_string_var.get() == 'paired':
                ctrl_rep2_r2_button.grid(column = 1, row = 4, padx = (10,5))
                ctrl_rep2_r2_entry.grid(column = 2, row = 4, padx = (5,10), ipady = 3)
    
        if bool(ctrl_rep3_r1_string_var.get()):
            ctrl_rep3_r1_button.grid(column = 1, row = 5, padx = (10,5), pady = (10,0))
            ctrl_rep3_r1_entry.grid(column = 2, row = 5, padx = (5,10), ipady = 3, pady = (10,0))
            
            if read_mode_string_var.get() == 'paired':
                ctrl_rep3_r2_button.grid(column = 1, row = 6, padx = (10,5))
                ctrl_rep3_r2_entry.grid(column = 2, row = 6, padx = (5,10), ipady = 3)

        if bool(ctrl_rep4_r1_string_var.get()):
            ctrl_rep4_r1_button.grid(column = 1, row = 7, padx = (10,5), pady = (10,0))
            ctrl_rep4_r1_entry.grid(column = 2, row = 7, padx = (5,10), ipady = 3, pady = (10,0))

            if read_mode_string_var.get() == 'paired':
                ctrl_rep4_r2_button.grid(column = 1, row = 8, padx = (10,5))
                ctrl_rep4_r2_entry.grid(column = 2, row = 8, padx = (5,10), ipady = 3)

        if bool(ctrl_rep5_r1_string_var.get()):
            ctrl_rep5_r1_button.grid(column = 1, row = 9, padx = (10,5), pady = (10,0))
            ctrl_rep5_r1_entry.grid(column = 2, row = 9, padx = (5,10), ipady = 3, pady = (10,0))

            if read_mode_string_var.get() == 'paired':
                ctrl_rep5_r2_button.grid(column = 1, row = 10, padx = (10,5))
                ctrl_rep5_r2_entry.grid(column = 2, row = 10, padx = (5,10), ipady = 3)

    if sample_table_question_string_var.get() == 'no':
        if ctrl_rep_number_int_var.get() >= 1:
            ctrl_rep1_r1_button.grid(column = 1, row = 1, padx = (10,5))
            ctrl_rep1_r1_entry.grid(column = 2, row = 1, padx = (5,10), ipady = 3)
            
            if read_mode_string_var.get() == 'paired':
                ctrl_rep1_r2_button.grid(column = 1, row = 2, padx = (10,5))
                ctrl_rep1_r2_entry.grid(column = 2, row = 2, padx = (5,10), ipady = 3)

        if ctrl_rep_number_int_var.get() >= 2:
            ctrl_rep2_r1_button.grid(column = 1, row = 3, padx = (10,5), pady = (10,0))
            ctrl_rep2_r1_entry.grid(column = 2, row = 3, padx = (5,10), ipady = 3, pady = (10,0))
            
            if read_mode_string_var.get() == 'paired':
                ctrl_rep2_r2_button.grid(column = 1, row = 4, padx = (10,5))
                ctrl_rep2_r2_entry.grid(column = 2, row = 4, padx = (5,10), ipady = 3)
    
        if ctrl_rep_number_int_var.get() >= 3:
            ctrl_rep3_r1_button.grid(column = 1, row = 5, padx = (10,5), pady = (10,0))
            ctrl_rep3_r1_entry.grid(column = 2, row = 5, padx = (5,10), ipady = 3, pady = (10,0))
            
            if read_mode_string_var.get() == 'paired':
                ctrl_rep3_r2_button.grid(column = 1, row = 6, padx = (10,5))
                ctrl_rep3_r2_entry.grid(column = 2, row = 6, padx = (5,10), ipady = 3)

        if ctrl_rep_number_int_var.get() >= 4:
            ctrl_rep4_r1_button.grid(column = 1, row = 7, padx = (10,5), pady = (10,0))
            ctrl_rep4_r1_entry.grid(column = 2, row = 7, padx = (5,10), ipady = 3, pady = (10,0))

            if read_mode_string_var.get() == 'paired':
                ctrl_rep4_r2_button.grid(column = 1, row = 8, padx = (10,5))
                ctrl_rep4_r2_entry.grid(column = 2, row = 8, padx = (5,10), ipady = 3)

        if ctrl_rep_number_int_var.get() >= 5:
            ctrl_rep5_r1_button.grid(column = 1, row = 9, padx = (10,5), pady = (10,0))
            ctrl_rep5_r1_entry.grid(column = 2, row = 9, padx = (5,10), ipady = 3, pady = (10,0))

            if read_mode_string_var.get() == 'paired':
                ctrl_rep5_r2_button.grid(column = 1, row = 10, padx = (10,5))
                ctrl_rep5_r2_entry.grid(column = 2, row = 10, padx = (5,10), ipady = 3)

    rearrange_entry_field_function()


ctrl_list_frame = tk.Frame(root)

ctrl_list_label = tk.Label(ctrl_list_frame, text = "Please type in, or browse and select your control sample files below", justify = tk.LEFT, width = 70)
ctrl_list_label.grid(row = 0, column = 1, padx = 10, pady = (5,10), columnspan = 2)

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

CreateToolTip(ctrl_rep1_r1_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam file extension accepted)')
CreateToolTip(ctrl_rep1_r2_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam file extension accepted)')
CreateToolTip(ctrl_rep2_r1_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam file extension accepted)')
CreateToolTip(ctrl_rep2_r2_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam file extension accepted)')
CreateToolTip(ctrl_rep3_r1_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam file extension accepted)')
CreateToolTip(ctrl_rep3_r2_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam file extension accepted)')
CreateToolTip(ctrl_rep4_r1_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam file extension accepted)')
CreateToolTip(ctrl_rep4_r2_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam file extension accepted)')
CreateToolTip(ctrl_rep5_r1_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam file extension accepted)')
CreateToolTip(ctrl_rep5_r2_button, text = 'Click here to browse and select your file\n(.fastq, .fastq.gz, .fq, .fq.gz, or .bam file extension accepted)')

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

ctrl_sample_notification_label = tk.Label(ctrl_list_frame, textvariable = ctrl_sample_notification_string_var, width = 70, padx = 5, pady = 5)
ctrl_sample_notification_label.grid(row = 15, column = 1, padx = 10, pady = 5, columnspan = 2)

ctrl_list_frame_back_button = tk.Button(ctrl_list_frame, text = "<< Back", command = lambda : change_frame_function(ctrl_list_frame, chip_list_frame), width = 20)
ctrl_list_frame_continue_button = tk.Button(ctrl_list_frame, text = "Continue >>", command = lambda : change_frame_function(ctrl_list_frame, setting_table_frame), width = 20)

ctrl_list_frame_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
ctrl_list_frame_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################


def setting_table_button_function():
    setting_table_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), mode = 'rb', title = 'Choose a setting table file')
    if setting_table_button_input:
        setting_table_string_var.set(setting_table_button_input.name)
        current_dir_string_var.set(os.path.dirname(setting_table_button_input.name))
        setting_table_loading_test_function()

    else:
        clear_setting_function()
        setting_table_question_string_var.set('no')
        setting_table_entry_popup_function('no')
        setting_table_notification_string_var.set('No custom settings table was selected, reverted to default')
        setting_table_notification_label.config(fg = 'red')

        rearrange_entry_field_function()


def setting_table_loading_test_function():
    try:
        setting_table_absolute_path = setting_table_string_var.get()
        setting_table_df = pd.read_csv(setting_table_absolute_path, delimiter='\t')
        setting_table_df.fillna('', inplace = True)
        setting_table_header = setting_table_df.columns.values.tolist() 

        # Check the formatting of the custom settings table, to ensure correct program-argument readings.
        # Check if the table consists of two columns 
        if len(setting_table_header) != 2:
            clear_setting_function()
            setting_table_question_string_var.set('no')
            setting_table_entry_popup_function('no')
            setting_table_notification_string_var.set('Columns error, reverted to default')
            setting_table_notification_label.config(fg = 'red')
            
        else:
            pass

        # Check if the table headers are 'program' and 'argument'. Check first if they are both strings to avoid TypeError.
        if isinstance(setting_table_header[0], str) and isinstance(setting_table_header[1], str):
            if setting_table_header[0].strip().lower() != 'program' or setting_table_header[1].strip().lower() != 'argument':
                clear_setting_function()
                setting_table_question_string_var.set('no')
                setting_table_entry_popup_function('no')
                setting_table_notification_string_var.set('Header error, reverted to default')
                setting_table_notification_label.config(fg = 'red')

            else:
                setting_table_notification_string_var.set('Custom settings table loading successful')
                setting_table_notification_label.config(fg = 'green')

        else:
            clear_setting_function()
            setting_table_question_string_var.set('no')
            setting_table_entry_popup_function('no')
            setting_table_notification_string_var.set('Header error, reverted to default')
            setting_table_notification_label.config(fg = 'red')

    except:
        clear_setting_function()
        setting_table_question_string_var.set('no')
        setting_table_entry_popup_function('no')
        setting_table_notification_string_var.set('Custom settings table loading error, reverted to default')
        setting_table_notification_label.config(fg = 'red')

    read_setting_table_function(setting_table_string_var.get())
    rearrange_entry_field_function()

    if setting_table_notification_string_var.get() == 'Custom settings table loading successful':
        return True
    else:
        return False


def read_setting_table_function(read_setting_table_arg):
    setting_table_absolute_path = read_setting_table_arg
    setting_table_df = pd.read_csv(setting_table_absolute_path, delimiter='\t')
    setting_table_df.fillna('', inplace = True)
    setting_table_header = setting_table_df.columns.values.tolist() 

    # Parse the location of 'program' and 'argument' columns
    setting_table_program_colnum    = setting_table_header.index('program')
    setting_table_argument_colnum   = setting_table_header.index('argument')
    setting_table_array             = setting_table_df.values.tolist()

    argument_dict = {}

    for suite_program in suite_program_list:
        argument_dict[suite_program] = []

    # For each entry in the 'program' column that matched with any of the program name in suite_program_list, 
    #   assign the argument as the value and program as the key in dictionary argument_dict
    for setting_table_array_row in setting_table_array:
        if setting_table_array_row[setting_table_program_colnum] in suite_program_list:
            current_setting_table_program = setting_table_array_row[setting_table_program_colnum]
            current_setting_table_argument = setting_table_array_row[setting_table_argument_colnum]
            argument_dict[current_setting_table_program].append(current_setting_table_argument)

    # Finally, join all arguments value within each program key with a single space and 
    #   assign the joined string into their own variable for easier calling later downstream
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


def setting_table_entry_popup_function(setting_table_entry_popup_arg):
    if setting_table_entry_popup_arg == 'yes':
        setting_table_entry.grid(row = 3, column = 1, padx = (5,10), ipady = 3)
        setting_table_button.grid(row = 3, column = 2, padx = (10,5), pady = 2)
        clear_setting_function()
        setting_table_notification_string_var.set('Please load your custom settings table file')
        setting_table_notification_label.config(fg = 'blue')
        root.eval('tk::PlaceWindow . center')
        root.grab_set()
        root.focus_force()

    if setting_table_entry_popup_arg == 'no':
        setting_table_entry.grid_remove()
        setting_table_button.grid_remove()
        setting_table_string_var.set(setting_table_file_full_path)
        read_setting_table_function(setting_table_file_full_path)
        setting_table_notification_string_var.set('Currently using default settings table')
        setting_table_notification_label.config(fg = 'blue')
        root.eval('tk::PlaceWindow . center')
        root.grab_set()
        root.focus_force()

    rearrange_entry_field_function()


def clear_setting_function():

    setting_table_string_var.set('')

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

    rearrange_entry_field_function()


setting_table_frame = tk.Frame(root)

setting_table_label = tk.Label(setting_table_frame, text = "Do you want to use your custom settings table?", width = 70)
setting_table_label.grid(row = 1, column = 1, padx = 10, pady = 5, columnspan = 2)

setting_table_yes_radio = tk.Radiobutton(setting_table_frame, text = "Yes, I would like use my custom settings table", padx = 5, variable = setting_table_question_string_var, value = 'yes', command = lambda : setting_table_entry_popup_function('yes'), width = 45)
CreateToolTip(setting_table_yes_radio, text = 'Select this to load a custom setting table in ChIP-AP format (see GitHub documentation).\nThe ChIP-AP pipeline will then be run based on the custom settings from loaded setting table.')
setting_table_yes_radio.grid(row = 2, column = 1, padx = 10, pady = 2, columnspan = 2)

setting_table_entry = tk.Entry(setting_table_frame, textvariable = setting_table_string_var, width = 50, justify = tk.RIGHT)
setting_table_button = tk.Button(setting_table_frame, text = 'Browse', command = lambda : setting_table_button_function(), width = 20)
CreateToolTip(setting_table_button, text = 'Click here to browse and select your setting table file.')

setting_table_no_radio = tk.Radiobutton(setting_table_frame, text = "No, I would like to use the default settings", padx = 5, variable = setting_table_question_string_var, value = 'no', command = lambda : setting_table_entry_popup_function('no'), width = 45)
CreateToolTip(setting_table_no_radio, text = 'Select this if you want to use the ChIP-AP pipeline default settings.')
setting_table_no_radio.grid(row = 4, column = 1, padx = 10, pady = 2, columnspan = 2)

setting_table_notification_label = tk.Label(setting_table_frame, textvariable = setting_table_notification_string_var, width = 70, padx = 5, pady = 5)
setting_table_notification_label.grid(row = 5, column = 1, padx = 10, pady = 5, columnspan = 2)

setting_table_back_button = tk.Button(setting_table_frame, text = "<< Back", command = lambda : change_frame_function(setting_table_frame, ctrl_list_frame), width = 20)
setting_table_continue_button = tk.Button(setting_table_frame, text = "Continue >>", command = lambda : change_frame_function(setting_table_frame, setting_value_frame), width = 20)

setting_table_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
setting_table_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################


setting_value_frame = tk.Frame(root)

setting_value_label = tk.Label(setting_value_frame, text = "Please modify the setting values below as needed", width = 70)
CreateToolTip(setting_value_label, text = "Warning: change the values below only when you know what you are doing.\nOtherwise, please just leave them at their default values.\nInvalid values will cause the whole pipeline run to break.\nCheck GitHub documentation to learn more about custom settings.")
setting_value_label.grid(row = 0, column = 1, pady = 5, columnspan = 2)

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

setting_value_back_button = tk.Button(setting_value_frame, text = "<< Back", command = lambda : change_frame_function(setting_value_frame, setting_table_frame), width = 22)
setting_value_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))

default_setting_button = tk.Button(setting_value_frame, text = "Restore defaults", command = lambda : setting_table_entry_popup_function('no'), width = 22)
CreateToolTip(default_setting_button, text = 'Click here if you want to restore to ChIP-AP pipeline default settings')
default_setting_button.grid(row = 21, column = 1, padx = 5, pady = (5,10), columnspan = 2)

setting_value_continue_button = tk.Button(setting_value_frame, text = "Continue >>", command = lambda : change_frame_function(setting_value_frame, genome_ref_frame), width = 22)
setting_value_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################


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

    rearrange_entry_field_function()


def genome_folder_button_function():
    genome_folder_button_input = filedialog.askdirectory(initialdir = current_dir_string_var.get(), parent = root, title = 'Choose the Genome Directory')
    if genome_folder_button_input:
        genome_folder_string_var.set(genome_folder_button_input)
        current_dir_string_var.set(genome_folder_button_input)

    rearrange_entry_field_function()


genome_ref_frame = tk.Frame(root)

genome_ref_label = tk.Label(genome_ref_frame, text = "Reference genome:")
genome_ref_label.grid(row = 1, column = 1, padx = (10,5), pady = (10,5))

genome_ref_drop_down = tk.OptionMenu(genome_ref_frame, genome_ref_string_var, *genome_ref_options, command = lambda x = None : auto_assign_genome_folder_function())
CreateToolTip(genome_ref_drop_down, text = 'Select the genome assembly you want your ChIP-seq reads to be aligned to. Currently, ChIP-AP supports six genome assemblies.\nOutside those supported by ChIP-AP, you will need to generate the reference files yourself')
genome_ref_drop_down.config(width = 45, takefocus = 1)

genome_ref_drop_down_menu = root.nametowidget(genome_ref_drop_down.menuname)
genome_ref_drop_down_menu.config(font = default_font)                   

genome_ref_drop_down.grid(row = 1, column = 2, padx = (5,10), pady = (10,5))

genome_folder_button = tk.Button(genome_ref_frame, text = 'Genome folder', command = lambda : genome_folder_button_function(), state = tk.NORMAL, width = 20)
CreateToolTip(genome_folder_button, text = 'Click here to browse and select the directory containing your custom genome reference')
genome_folder_entry = tk.Entry(genome_ref_frame, textvariable = genome_folder_string_var, width = 50, justify = tk.RIGHT)

genome_ref_back_button = tk.Button(genome_ref_frame, text = "<< Back", command = lambda : change_frame_function(genome_ref_frame, setting_value_frame), width = 20)
genome_ref_continue_button = tk.Button(genome_ref_frame, text = "Continue >>", command = lambda : change_frame_function(genome_ref_frame, known_motif_frame), width = 20)

genome_ref_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
genome_ref_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################


def known_motif_button_function():
    known_motif_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a .motif file')
    if known_motif_button_input:
        known_motif_string_var.set(known_motif_button_input.name)
        current_dir_string_var.set(os.path.dirname(known_motif_button_input.name))

    rearrange_entry_field_function()


known_motif_frame = tk.Frame(root)

known_motif_button = tk.Button(known_motif_frame, text = 'Known motif file', command = lambda : known_motif_button_function(), state = tk.NORMAL, width = 20)
CreateToolTip(known_motif_button, text = 'Click here to browse and select your .motif file (in HOMER matrix format)')
known_motif_button.grid(row = 1, column = 1, padx = (10,5), pady = (10,5))

known_motif_entry = tk.Entry(known_motif_frame, textvariable = known_motif_string_var, width = 50, justify = tk.RIGHT)
known_motif_entry.grid(row = 1, column = 2, padx = (5,10), pady = (10,5), ipady = 3)

known_motif_back_button = tk.Button(known_motif_frame, text = "<< Back", command = lambda : change_frame_function(known_motif_frame, genome_ref_frame), width = 20)
known_motif_continue_button = tk.Button(known_motif_frame, text = "Continue >>", command = lambda : change_frame_function(known_motif_frame, output_folder_frame), width = 20)

known_motif_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
known_motif_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################


def output_folder_button_function():
    output_folder_button_input = filedialog.askdirectory(initialdir = current_dir_string_var.get(), parent = root, title = 'Choose the Output Directory')
    if output_folder_button_input:
        output_folder_string_var.set(output_folder_button_input)
        current_dir_string_var.set(output_folder_button_input)

    rearrange_entry_field_function()


output_folder_frame = tk.Frame(root)

output_folder_button = tk.Button(output_folder_frame, text = 'Output save folder', command = lambda : output_folder_button_function(), state = tk.NORMAL, width = 20)
CreateToolTip(output_folder_button, text = 'Click here to browse and select the directory to save the results of your ChIP-AP pipeline run.\nTo save in a new folder, type in the entry field your new folder name after the output save directory:\nfull_path_to_output_save_directory/new_folder_name')
output_folder_button.grid(row = 1, column = 1, padx = (10,5), pady = (10,5))

output_folder_entry = tk.Entry(output_folder_frame, textvariable = output_folder_string_var, width = 50, justify = tk.RIGHT)
output_folder_entry.grid(row = 1, column = 2, padx = (5,10), pady = (10,5), ipady = 3)

setname_label = tk.Label(output_folder_frame, text = "Dataset prefix:", width = 20)
setname_label.grid(row = 2, column = 1, padx = (10,5), pady = 5)

setname_entry = tk.Entry(output_folder_frame, textvariable = setname_string_var, width = 50, justify = tk.LEFT)
CreateToolTip(setname_entry, text = 'Type in your desired folder name and prefix for all the resulting output filenames.')
setname_entry.grid(row = 2, column = 2, padx = (5,10), pady = 5, ipady = 3)

output_folder_back_button = tk.Button(output_folder_frame, text = "<< Back", command = lambda : change_frame_function(output_folder_frame, known_motif_frame), width = 20)
output_folder_continue_button = tk.Button(output_folder_frame, text = "Continue >>", command = lambda : change_frame_function(output_folder_frame, checkbox_frame), width = 20)

output_folder_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
output_folder_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################


checkbox_frame = tk.Frame(root)

fcmerge_checkbox = tk.Checkbutton(checkbox_frame, text = ' Merged fold enrichment analysis', variable = fcmerge_var, onvalue = 1, offvalue = 0)
CreateToolTip(fcmerge_checkbox, text = 'Check this box if you want the fold enrichment analysis from all replicates combined as one.\nThis option is automatically chosen when there are unequal number of replicates between ChIP and control samples.')
fcmerge_checkbox.grid(sticky = "w", row = 1, column = 1, padx = 10, columnspan = 2)

goann_checkbox = tk.Checkbutton(checkbox_frame, text = ' Annotate peaks with known gene ontology terms', variable = goann_var, onvalue = 1, offvalue = 0)
CreateToolTip(goann_checkbox, text = 'Check this box if you want each peak in the final detected peaks list to have gene ontology annotations')
goann_checkbox.grid(sticky = "w", row = 2, column = 1, padx = 10, columnspan = 2)

pathann_checkbox = tk.Checkbutton(checkbox_frame, text = ' Annotate peaks with known pathway terms', variable = pathann_var, onvalue = 1, offvalue = 0)
CreateToolTip(pathann_checkbox, text = 'Check this box if you want each peak in the final detected peaks list to have pathway annotations')
pathann_checkbox.grid(sticky = "w", row = 3, column = 1, padx = 10, columnspan = 2)

deltemp_checkbox = tk.Checkbutton(checkbox_frame, text = ' Automatically delete large temporary files', variable = deltemp_var, onvalue = 1, offvalue = 0)
CreateToolTip(deltemp_checkbox, text = 'Check this box if you will not need the large-sized intermediary files.')
deltemp_checkbox.grid(sticky = "w", row = 4, column = 1, padx = 10, columnspan = 2)

stdout_checkbox = tk.Checkbutton(checkbox_frame, text = ' Record standard outputs', variable = stdout_var, onvalue = 1, offvalue = 0)
CreateToolTip(stdout_checkbox, text = 'Check this box if you want to save pipeline outputs (channel 1>) in a text file.')
stdout_checkbox.grid(sticky = "w", row = 5, column = 1, padx = 10, columnspan = 2)

stderr_checkbox = tk.Checkbutton(checkbox_frame, text = ' Record standard errors', variable = stderr_var, onvalue = 1, offvalue = 0)
CreateToolTip(stderr_checkbox, text = 'Check this box if you want to save pipeline errors (channel 2>) in a text file.')
stderr_checkbox.grid(sticky = "w", row = 6, column = 1, padx = 10, columnspan = 2)

checkbox_back_button = tk.Button(checkbox_frame, text = "<< Back", command = lambda : change_frame_function(checkbox_frame, output_folder_frame), width = 20)
checkbox_continue_button = tk.Button(checkbox_frame, text = "Continue >>", command = lambda : change_frame_function(checkbox_frame, execute_frame), width = 20)

checkbox_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
checkbox_continue_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################


def print_sample_table_function():
    
    sample_table_output_df = pd.DataFrame.from_dict(sample_table_output_dict, orient = 'index') 
    sample_table_output_df = sample_table_output_df.transpose()

    sample_table_output_df.to_csv('{}/{}_sample_table.tsv'.format(output_dir.get(), setname_string_var.get()), sep = '\t', index = False)


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
        fold_change_calculator_arg.get()]

    setting_table_output_dict = {'program' : suite_program_list, 'argument' : suite_program_arg}

    setting_table_output_df = pd.DataFrame.from_dict(setting_table_output_dict)
    
    setting_table_output_df.to_csv('{}/{}_setting_table.tsv'.format(output_dir.get(), setname_string_var.get()), sep = '\t', index = False)


def print_GS_command_line_function():
    with open('{}/{}_command_line.txt'.format(output_dir.get(), setname_string_var.get()), 'w') as command_line_file:
        command_line_file.write(command_line_output_string_var.get() + '\n')


def print_GSaR_command_line_function():
    with open('{}/{}_command_line.txt'.format(output_dir.get(), setname_string_var.get()), 'w') as command_line_file:
        command_line_file.write(command_line_output_string_var.get() + ' --run\n')


def execute_GS_command_line_function():
    subprocess.run(command_line_output_string_var.get(), shell = True)


def execute_GSaR_command_line_function():
    subprocess.run(command_line_output_string_var.get() + ' --run', shell = True)


def generate_scripts_function():
    if not os.path.exists(output_dir.get()):
        os.makedirs(output_dir.get())
    
    print_sample_table_function()
    print_setting_table_function()
    print_GS_command_line_function()
    execute_GS_command_line_function()
    execute_exit_button.config(state = tk.NORMAL)
    cpu_count_notification_string_var.set("Scripts generated! Check them out in:\n{}".format(output_dir.get()))
    cpu_count_notification_label.config(fg = 'green')

def generate_scripts_and_run_function():
    if not os.path.exists(output_dir.get()):
        os.makedirs(output_dir.get())
    
    print_sample_table_function()
    print_setting_table_function()
    print_GSaR_command_line_function()
    execute_exit_button.config(state = tk.NORMAL)
    cpu_count_notification_string_var.set("Pipeline started! Check your results later in:\n{}".format(output_dir.get()))
    cpu_count_notification_label.config(fg = 'green')
    execute_GSaR_command_line_function()


execute_frame = tk.Frame(root)

cpu_count_label = tk.Label(execute_frame, text = "CPU cores to use:", width = 20)
cpu_count_label.grid(row = 1, column = 1, padx = (10,5))

cpu_count_entry = tk.Entry(execute_frame, width = 20, textvariable = cpu_count_string_var)
CreateToolTip(cpu_count_entry, text = 'Type in your desired number of CPU cores to be used by the pipeline')
cpu_count_entry.grid(row = 1, column = 2, padx = (5,10), ipady = 3)

cpu_count_notification_label = tk.Label(execute_frame, textvariable = cpu_count_notification_string_var, width = 50)
cpu_count_notification_label.grid(row = 2, column = 1, padx = 10, columnspan = 2)

generate_scripts_button = tk.Button(execute_frame, text = 'Generate scripts', command = lambda : generate_scripts_function(), width = 25)
CreateToolTip(generate_scripts_button, text = 'Select this if you wish to run the pipeline later\nby executing MASTER_script.sh within the output save folder')
generate_and_run_scripts_button = tk.Button(execute_frame, text = 'Generate and run scripts', command = lambda : generate_scripts_and_run_function(), width = 25)
CreateToolTip(generate_and_run_scripts_button, text = 'Select this if you wish to run the pipeline now\nNOTE: It may take up to several hours depending on your system')

generate_scripts_button.grid(row = 20, column = 1, sticky = "w", padx = (10,5))
generate_and_run_scripts_button.grid(row = 20, column = 2, sticky = "e", padx = (5,10))

execute_back_button = tk.Button(execute_frame, text = "<< Back", command = lambda : change_frame_function(execute_frame, checkbox_frame), width = 20)
execute_exit_button = tk.Button(execute_frame, text = "Exit wizard", command = lambda : exit(), width = 20, state = tk.DISABLED)

execute_back_button.grid(sticky = "w", row = 21, column = 1, padx = (10,5), pady = (5,10))
execute_exit_button.grid(sticky = "e", row = 21, column = 2, padx = (5,10), pady = (5,10))


########################################################################################################################


register_sample_function()
update_command_line_function()
rearrange_entry_field_function()

if __name__ == "__main__":
    root.mainloop()


########################################################################################################################
#!/usr/bin/env python3
#pyright: reportUnboundVariable=false


# script_version = '2.0'


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

chipap_icon_full_path = os.path.expanduser('{}/ChIP-AP_icon_GUI.png'.format(sys.path[0]))
chipap_logo_full_path = os.path.expanduser('{}/ChIP-AP_logo_GUI.jpg'.format(sys.path[0]))
default_current_dir = os.path.expanduser('~')
genome_folder_full_path = os.path.expanduser('~/genomes')
default_setting_table_file_full_path = '{}/default_settings_table.tsv'.format(genome_folder_full_path)


root = tk.Tk()
root.title(chipap_program_name)
root.resizable(width = False, height = False)

line_style = ttk.Style()
line_style.configure("Line.TSeparator", background = 'black')

ttk.Separator(root, style = "Line.TSeparator", orient = tk.VERTICAL).grid(row = 0, column = 15, rowspan = 6, sticky = "ns")
ttk.Separator(root, style = "Line.TSeparator", orient = tk.VERTICAL).grid(row = 0, column = 10, rowspan = 11, sticky = "ns")
ttk.Separator(root, style = "Line.TSeparator", orient = tk.VERTICAL).grid(row = 30, column = 6, rowspan = 17, sticky = "ns")
ttk.Separator(root, style = "Line.TSeparator", orient = tk.VERTICAL).grid(row = 30, column = 12, rowspan = 17, sticky = "ns")

default_font = tkFont.nametofont("TkDefaultFont")
default_font.configure(family = 'fixed', size = 16)

text_font = tkFont.nametofont("TkTextFont")
text_font.configure(family = 'fixed', size = 16)

fixed_font = tkFont.nametofont("TkFixedFont")
fixed_font.configure(family = 'fixed', size = 16)


root.tk.call('wm', 'iconphoto', root._w, tk.PhotoImage(file = chipap_icon_full_path))


argument_dict = {}

valid_extension_list = [".fastq",
                        ".fq",
                        ".fastq.gz",
                        ".fq.gz",
                        ".bam"]

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

genome_ref_options = ["hg38 (Homo sapiens)", 
                        "hg19 (Homo sapiens)", 
                        "mm9 (Mus musculus)", 
                        "mm10 (Mus musculus)", 
                        "dm6 (Drosophila melanogaster)", 
                        "sacCer3 (Saccharomyces cerevisiae)",
                        "Other [!!!under construction!!!]"]

homer_motif_options = ["None",
                        "Consensus peak set", 
                        "Union peak set", 
                        "Both peak sets"]

meme_motif_options = ["None",
                        "Consensus peak set", 
                        "Union peak set", 
                        "Both peak sets"]


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
    error_type_A_int_var.set(0)

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
                error_type_A_int_var.set(1)
                return

        if read_mode_string_var.get() == 'paired':
            if all(chip_r1 == '' for chip_r1 in chip_list_r1):
                chip_sample_notification_string_var.set('Please assign ChIP sample (read 1)')
                chip_sample_notification_label.config(fg = 'blue')
                error_type_A_int_var.set(1)
                return

            # Check for r2 only if the input files are not aligned
            if all(chip_r2 == '' for chip_r2 in chip_list_r2) and not all('.bam' in sample for sample in (chip_list_r1 + chip_list_r1) if sample != ''):
                chip_sample_notification_string_var.set('Please assign ChIP sample (read 2)')
                chip_sample_notification_label.config(fg = 'blue')
                error_type_A_int_var.set(1)
                return

        if read_mode_string_var.get() == 'single':
            if not bool(chip_rep1_r1_string_var.get()):
                chip_sample_notification_string_var.set('ChIP samples are not assigned from the top (replicate 1)')
                chip_sample_notification_label.config(fg = 'red')
                chip_rep1_r1_entry.config(bg = 'IndianRed1')
                error_type_A_int_var.set(1)
                return
        
        if read_mode_string_var.get() == 'paired':
            if not bool(chip_rep1_r1_string_var.get()) and not bool(chip_rep1_r2_string_var.get()):
                chip_sample_notification_string_var.set('ChIP samples are not assigned from the top (replicate 1)')
                chip_sample_notification_label.config(fg = 'red')
                chip_rep1_r1_entry.config(bg = 'IndianRed1')
                chip_rep1_r2_entry.config(bg = 'IndianRed1')
                error_type_A_int_var.set(1)
                return

            if not bool(chip_rep1_r1_string_var.get()):
                chip_sample_notification_string_var.set('ChIP samples are not assigned from the top (replicate 1)')
                chip_sample_notification_label.config(fg = 'red')
                chip_rep1_r1_entry.config(bg = 'IndianRed1')
                error_type_A_int_var.set(1)
                return

            if not bool(chip_rep1_r2_string_var.get()):
                chip_sample_notification_string_var.set('ChIP samples are not assigned from the top (replicate 1)')
                chip_sample_notification_label.config(fg = 'red')
                chip_rep1_r2_entry.config(bg = 'IndianRed1')
                error_type_A_int_var.set(1)
                return

        chip_sample_notification_string_var.set('ChIP samples have been assigned')
        chip_sample_notification_label.config(fg = 'green')
    
        check_chip_sample_file_pair_function()

    else:
        return


def check_ctrl_sample_assigned_function():

    if read_mode_string_var.get() == 'single' or read_mode_string_var.get() == 'paired':
        
        if read_mode_string_var.get() == 'single':
            if all(ctrl_r1 == '' for ctrl_r1 in ctrl_list_r1):
                ctrl_sample_notification_string_var.set('Please assign control sample')
                ctrl_sample_notification_label.config(fg = 'blue')
                error_type_A_int_var.set(1)
                return

        if read_mode_string_var.get() == 'paired':
            if all(ctrl_r1 == '' for ctrl_r1 in ctrl_list_r1):
                ctrl_sample_notification_string_var.set('Please assign control sample (read 1)')
                ctrl_sample_notification_label.config(fg = 'blue')
                error_type_A_int_var.set(1)
                return

            # Check for r2 only if the input files are not aligned
            if all(ctrl_r2 == '' for ctrl_r2 in ctrl_list_r2) and not all('.bam' in sample for sample in (ctrl_list_r1 + ctrl_list_r1) if sample != ''):
                ctrl_sample_notification_string_var.set('Please assign control sample (read 2)')
                ctrl_sample_notification_label.config(fg = 'blue')
                error_type_A_int_var.set(1)
                return

        if read_mode_string_var.get() == 'single':
            if not bool(ctrl_rep1_r1_string_var.get()):
                ctrl_sample_notification_string_var.set('Control samples are not assigned from the top (replicate 1)')
                ctrl_sample_notification_label.config(fg = 'red')
                ctrl_rep1_r1_entry.config(bg = 'IndianRed1')
                error_type_A_int_var.set(1)
                return
        
        if read_mode_string_var.get() == 'paired':
            if not bool(ctrl_rep1_r1_string_var.get()) and not bool(ctrl_rep1_r2_string_var.get()):
                ctrl_sample_notification_string_var.set('Control samples are not assigned from the top (replicate 1)')
                ctrl_sample_notification_label.config(fg = 'red')
                ctrl_rep1_r1_entry.config(bg = 'IndianRed1')
                ctrl_rep1_r2_entry.config(bg = 'IndianRed1')
                error_type_A_int_var.set(1)
                return
            
            elif not bool(ctrl_rep1_r1_string_var.get()):
                ctrl_sample_notification_string_var.set('Control samples are not assigned from the top (replicate 1)')
                ctrl_sample_notification_label.config(fg = 'red')
                ctrl_rep1_r1_entry.config(bg = 'IndianRed1')
                error_type_A_int_var.set(1)
                return

            elif not bool(ctrl_rep1_r2_string_var.get()):
                ctrl_sample_notification_string_var.set('Control samples are not assigned from the top (replicate 1)')
                ctrl_sample_notification_label.config(fg = 'red')
                ctrl_rep1_r2_entry.config(bg = 'IndianRed1')
                error_type_A_int_var.set(1)
                return

        ctrl_sample_notification_string_var.set('Control samples have been assigned')
        ctrl_sample_notification_label.config(fg = 'green')

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
            error_type_A_int_var.set(1)
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
            error_type_A_int_var.set(1)
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
        error_type_A_int_var.set(1)
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
        error_type_A_int_var.set(1)
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
        error_type_A_int_var.set(1)

    elif chip_file_exist_error_state == 0:
        chip_sample_notification_string_var.set('No problem found in ChIP samples')
        chip_sample_notification_label.config(fg = 'green')


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
        error_type_A_int_var.set(1)

    elif ctrl_file_exist_error_state == 0:
        ctrl_sample_notification_string_var.set('No problem found in control samples')
        ctrl_sample_notification_label.config(fg = 'green')


def check_required_input_function(*args):
    
    error_type_B_int_var.set(0)

    cpu_count_entry.config(fg = 'black')

    if read_mode_string_var.get() != 'single' and read_mode_string_var.get() != 'paired':
        error_type_B_int_var.set(1)

    if peak_type_string_var.get() != 'narrow' and peak_type_string_var.get() != 'broad' and peak_type_string_var.get() != 'unsure':
        error_type_B_int_var.set(1)

    if bool(sample_table_string_var.get()):
        if not os.path.isfile(sample_table_string_var.get()):
            error_type_B_int_var.set(1)

    if bool(setting_table_string_var.get()):
        if '[MODIFIED]' in setting_table_string_var.get():
            pass
        elif not os.path.isfile(setting_table_string_var.get()):
            error_type_B_int_var.set(1)

    if not bool(genome_ref_string_var.get()):
        error_type_B_int_var.set(1)
    elif not bool(genome_folder_string_var.get()):
        error_type_B_int_var.set(1)

    if bool(known_motif_string_var.get()):
        if not os.path.isfile(known_motif_string_var.get()):
            error_type_B_int_var.set(1)

    if not bool(setname_string_var.get()):
        error_type_B_int_var.set(1)
    elif not bool(output_folder_string_var.get()):
        error_type_B_int_var.set(1)

    if not bool(cpu_count_string_var.get()):
        error_type_B_int_var.set(1)
        cpu_count_entry.config(fg = 'black')
        cpu_count_notification_string_var.set("Maximum number of CPU cores available: {}".format(max_cpu))
        cpu_count_notification_label.config(fg = 'blue')
    
    else:
        if int(cpu_count_string_var.get()) > max_cpu:
            error_type_B_int_var.set(1)
            cpu_count_entry.config(fg = 'red')
            cpu_count_notification_string_var.set("Entered number exceeds available CPU cores ({})".format(max_cpu))
            cpu_count_notification_label.config(fg = 'red')

        elif int(cpu_count_string_var.get()) < 1:
            error_type_B_int_var.set(1)
            cpu_count_entry.config(fg = 'red')
            cpu_count_notification_string_var.set("Need at least one CPU core to run the pipeline")
            cpu_count_notification_label.config(fg = 'red')

        elif int(cpu_count_string_var.get()) <= max_cpu:
            cpu_count_entry.config(fg = 'black')
            cpu_count_notification_string_var.set("")
            cpu_count_notification_label.config(fg = 'blue')

    generate_button_switch_function()


def generate_button_switch_function():
    if error_type_A_int_var.get() == 0 and error_type_B_int_var.get() == 0:
        generate_scripts_button.config(state = tk.NORMAL)
        generate_and_run_scripts_button.config(state = tk.NORMAL)
        cpu_count_notification_string_var.set('ChIP-AP ready to go!')
        cpu_count_notification_label.config(fg = 'green')
    else:
        generate_scripts_button.config(state = tk.DISABLED)
        generate_and_run_scripts_button.config(state = tk.DISABLED)


def update_command_line_function(*args):

    if read_mode_string_var.get() == 'single' or read_mode_string_var.get() == 'paired':
        read_mode_arg.set(' --mode {}'.format(read_mode_string_var.get()))
    else:
        read_mode_arg.set('')

    if peak_type_string_var.get() == 'narrow' or peak_type_string_var.get() == 'broad' or peak_type_string_var.get() == 'unsure':
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

    if homer_motif_string_var.get() != 'None':
        homer_motif_arg.set(' --homer_motif {}'.format(homer_motif_string_var.get().split(' ')[0].lower()))
    elif homer_motif_string_var.get() == 'None':
        homer_motif_arg.set('')

    if meme_motif_string_var.get() != 'None':
        meme_motif_arg.set(' --meme_motif {}'.format(meme_motif_string_var.get().split(' ')[0].lower()))
    elif meme_motif_string_var.get() == 'None':
        meme_motif_arg.set('')
        
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
        stderr_arg.set(' 2> {}/{}.err'.format(output_dir.get(), setname_string_var.get()))
    else:
        stdout_arg.set('')
        stderr_arg.set('')


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


def check_traced_input_function(*args):
    check_required_input_function()
    update_command_line_function()


def sample_table_entry_trace_load_function(*args):
    if bool(sample_table_string_var.get()):
        if not os.path.isfile(sample_table_string_var.get()):
            sample_table_notification_string_var.set('Sample table not found')
            sample_table_notification_label.config(fg = 'red')
            clear_sample_function('valuesonly')
        elif os.path.isfile(sample_table_string_var.get()):
            sample_table_loading_test_function()


def setting_table_entry_trace_load_function(*args):
    if bool(setting_table_string_var.get()):
        if '[MODIFIED]' in setting_table_string_var.get():
            pass
        elif not os.path.isfile(setting_table_string_var.get()):
            setting_table_notification_string_var.set('Setting table not found')
            setting_table_notification_label.config(fg = 'red')
            clear_setting_function('valuesonly')
            read_setting_table_function(default_setting_table_file_full_path)
        elif os.path.isfile(setting_table_string_var.get()):
            setting_table_loading_test_function()


########################################################################################################################

max_cpu = multiprocessing.cpu_count()

current_dir_string_var = tk.StringVar(value = default_current_dir)

error_type_A_int_var = tk.IntVar(value = 0)
error_type_B_int_var = tk.IntVar(value = 0)

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

sample_table_string_var = tk.StringVar()
sample_table_notification_string_var = tk.StringVar()
chip_sample_notification_string_var = tk.StringVar()
ctrl_sample_notification_string_var = tk.StringVar()

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

setting_table_string_var = tk.StringVar(value = default_setting_table_file_full_path)
setting_table_notification_string_var = tk.StringVar(value = "Currently using default settings table")

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
# stderr_var = tk.IntVar()

homer_motif_string_var = tk.StringVar(value = homer_motif_options[0])
meme_motif_string_var = tk.StringVar(value = meme_motif_options[0])

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
homer_motif_arg = tk.StringVar()
meme_motif_arg = tk.StringVar()

output_dir = tk.StringVar()

########################################################################################################################

image = Image.open(chipap_logo_full_path)
image = image.resize((300, 110), Image.ANTIALIAS)
photo = ImageTk.PhotoImage(image)

chipap_logo = tk.Label(root, image = photo, anchor = "center")
chipap_logo.image = photo
chipap_logo.grid(row = 1, column = 1, rowspan = 4, columnspan = 3, sticky = "w",padx = 3)

chipap_about = tk.Text(root, width = 61, height = 7, relief = tk.FLAT, bg = 'gray85',font = (None, 9))
chipap_about_text = 'Integrated Analysis Pipeline for Unbiased ChIP-seq Analysis\n\nComplete guides and walkthroughs can be found on our github\n(https://github.com/JSuryatenggara/ChIP-AP)\n\nIf you use ChIP-AP please cite our bioRxiv pre-print article\n(https://www.biorxiv.org/content/10.1101/2021.04.18.440382v1)'
chipap_about.insert(tk.END, chipap_about_text)
chipap_about.config(state = tk.DISABLED)
chipap_about.tag_configure("center", justify = tk.CENTER)
chipap_about.tag_add("center", 1.0, tk.END)
chipap_about.grid(row = 1, column = 3, rowspan = 4, columnspan = 7, sticky = "e", padx = 5, pady = 2)

########################################################################################################################

def sample_button_switch_function():

    if read_mode_string_var.get() == 'single':
        r1_state = tk.NORMAL
        r2_state = tk.DISABLED
        sample_table_button_state = tk.NORMAL

    elif read_mode_string_var.get() == 'paired':
        r1_state = tk.NORMAL
        r2_state = tk.NORMAL
        sample_table_button_state = tk.NORMAL

    else:
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

    check_required_input_function()
    rearrange_entry_field_function()

    if read_mode_string_var.get():
        sample_table_notification_string_var.set('Please load sample table or assign samples manually below')
        sample_table_notification_label.config(fg = 'blue')

read_mode_label = tk.Label(root, text = "Dataset sequencing mode:", justify = tk.LEFT, width = 30)
read_mode_label.grid(row = 1, column = 11, columnspan = 4, padx = 10, pady = (5,2))

single_end_radio = tk.Radiobutton(root, text = "Single end", padx = 10, variable = read_mode_string_var, value = 'single', width = 30, command = lambda : sample_button_switch_function())
CreateToolTip(single_end_radio, text = 'Select this if there is one\noutput file per sample')
single_end_radio.grid(row = 2, column = 11, columnspan = 4, padx = 10)

paired_end_radio = tk.Radiobutton(root, text = "Paired end", padx = 10, variable = read_mode_string_var, value = 'paired', width = 30, command = lambda : sample_button_switch_function())
CreateToolTip(paired_end_radio, text = 'Select this if there are two\noutput files per sample.\nThey are typically in pairs:\nR1 and R2 for every sample')
paired_end_radio.grid(row = 3, column = 11, columnspan = 4, padx = 10, pady = (0,5))

########################################################################################################################

def mea_switch_function():
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

    rearrange_entry_field_function()

peak_type_label = tk.Label(root, text = "Dataset peak type:", justify = tk.LEFT, width = 30)
peak_type_label.grid(row = 1, column = 16, columnspan = 4, padx = 10, pady = (5,2))

narrow_peak_radio = tk.Radiobutton(root, text = "Narrow peaks", padx = 10, variable = peak_type_string_var, value = 'narrow', width = 30, command = lambda : mea_switch_function())
CreateToolTip(narrow_peak_radio, text = 'Select this for ChIP-seq experiment\nusing transcription factor protein')
narrow_peak_radio.grid(row = 2, column = 16, columnspan = 4, padx = 10)

broad_peak_radio = tk.Radiobutton(root, text = "Broad peaks", padx = 10, variable = peak_type_string_var, value = 'broad', width = 30, command = lambda : mea_switch_function())
CreateToolTip(broad_peak_radio, text = 'Select this for ChIP-seq experiment\nusing chromatin modifier protein')
broad_peak_radio.grid(row = 3, column = 16, columnspan = 4, padx = 10)

unsure_peak_radio = tk.Radiobutton(root, text = "Unsure", padx = 10, variable = peak_type_string_var, value = 'unsure', width = 30, command = lambda : mea_switch_function())
CreateToolTip(unsure_peak_radio, text = 'Select this for ChIP-seq experiment\nusing protein with both possibilities')
unsure_peak_radio.grid(row = 4, column = 16, columnspan = 4, padx = 10, pady = (0,5))

ttk.Separator(root, style = "Line.TSeparator", orient = tk.HORIZONTAL).grid(row = 5, column = 0, columnspan = 21, sticky = "we")

########################################################################################################################

def sample_table_button_function():
    sample_table_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a sample table file')
    if sample_table_button_input:
        sample_table_string_var.set(sample_table_button_input.name)
        current_dir_string_var.set(os.path.dirname(sample_table_button_input.name))
        sample_table_loading_test_function()

    else:
        clear_sample_function('withfile')
        sample_table_notification_string_var.set('No sample table was selected. No samples were loaded')
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
        clear_sample_function('withfile')
        sample_table_notification_string_var.set('Sample table loading error! No samples were loaded')
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
        clear_sample_function('withfile')
        sample_table_notification_string_var.set('Sample table reading error! No samples were loaded')
        sample_table_notification_label.config(fg = 'red')
        return False


    if read_mode_string_var.get() == 'single':

        if not all(chip_r2 == '' for chip_r2 in chip_r2_sample_list) or not all(ctrl_r2 == '' for ctrl_r2 in ctrl_r2_sample_list):
            clear_sample_function('withfile')
            sample_table_notification_string_var.set('Error! Paired samples detected. No samples were loaded')
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
            
            sample_table_notification_string_var.set('Sample table loading successful')
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

        sample_table_notification_string_var.set('Sample table loading successful')
        sample_table_notification_label.config(fg = 'green')
        
        rearrange_entry_field_function()


def clear_sample_function(clear_sample_function_arg):
    
    if clear_sample_function_arg == 'withfile':
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


sample_table_notification_label = tk.Label(root, textvariable = sample_table_notification_string_var, width = 70, padx = 10, pady = 4, fg = 'blue')
sample_table_notification_string_var.set('Please specify data sequencing mode first!')
sample_table_notification_label.grid(row = 6, column = 1, padx = 10, pady = 4, columnspan = 9)

sample_table_button = tk.Button(root, text = 'Load sample table', command = lambda:sample_table_button_function(), state = tk.DISABLED, width = 15)
CreateToolTip(sample_table_button, text = 'Click here to browse and\nselect your sample table file')
sample_table_button.grid(sticky = "we", row = 7, column = 1, padx = 10, rowspan = 1, columnspan = 1)

sample_table_entry = tk.Entry(root, textvariable = sample_table_string_var, width = 50, justify = tk.RIGHT, state = tk.DISABLED)
sample_table_entry.grid(sticky = "we", row = 7, column = 2, padx = (0,10), columnspan = 8, ipady = 3)

clear_sample_button = tk.Button(root, text = 'Clear samples', command = lambda:clear_sample_function('withfile'), state = tk.DISABLED, width = 15)
CreateToolTip(clear_sample_button, text = 'Click here to clear all\nassigned samples below')
clear_sample_button.grid(sticky = "we", row = 8, column = 1, padx = 10, pady = (4,5), rowspan = 1, columnspan = 1)

########################################################################################################################

def setting_table_button_function():
    setting_table_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), mode = 'rb', title = 'Choose a setting table file')
    if setting_table_button_input:
        setting_table_string_var.set(setting_table_button_input.name)
        current_dir_string_var.set(os.path.dirname(setting_table_button_input.name))
        setting_table_loading_test_function()

    else:
        return


def setting_table_loading_test_function():
    try:
        setting_table_absolute_path = setting_table_string_var.get()
        setting_table_df = pd.read_csv(setting_table_absolute_path, delimiter='\t')
        setting_table_df.fillna('', inplace = True)
        setting_table_header = setting_table_df.columns.values.tolist() 

        # Check the formatting of the custom settings table, to ensure correct program-argument readings.
        # Check if the table consists of two columns 
        if len(setting_table_header) != 2:
            clear_setting_function('withfile')
            default_setting_button_function()
            setting_table_notification_string_var.set('Columns error, reverted to default')
            setting_table_notification_label.config(fg = 'red')
            
        else:
            pass

        # Check if the table headers are 'program' and 'argument'. Check first if they are both strings to avoid TypeError.
        if isinstance(setting_table_header[0], str) and isinstance(setting_table_header[1], str):
            if setting_table_header[0].strip().lower() != 'program' or setting_table_header[1].strip().lower() != 'argument':
                clear_setting_function('withfile')
                default_setting_button_function()
                setting_table_notification_string_var.set('Header error, reverted to default')
                setting_table_notification_label.config(fg = 'red')

            else:
                if setting_table_string_var.get() == default_setting_table_file_full_path:
                    setting_table_notification_string_var.set('Currently using default settings table')
                    setting_table_notification_label.config(fg = 'blue')
                else:
                    setting_table_notification_string_var.set('Custom settings table loading successful')
                    setting_table_notification_label.config(fg = 'green')

        else:
            clear_setting_function('withfile')
            default_setting_button_function()
            setting_table_notification_string_var.set('Header error, reverted to default')
            setting_table_notification_label.config(fg = 'red')

    except:
        clear_setting_function('withfile')
        default_setting_button_function()
        setting_table_notification_string_var.set('Custom settings table loading error, reverted to default')
        setting_table_notification_label.config(fg = 'red')

    read_setting_table_function(setting_table_string_var.get())
    rearrange_entry_field_function()

    if setting_table_notification_string_var.get() == 'Custom settings table loading successful':
        return True
    if setting_table_notification_string_var.get() == 'Currently using default settings table':
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

    global argument_dict
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
    homer_findMotifsGenome_arg.set(' '.join(argument_dict['homer_findMotifsGenome']))
    meme_chip_arg.set(' '.join(argument_dict['meme_chip']))

def clear_setting_function(clear_setting_function_arg):

    if clear_setting_function_arg == 'withfile':
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
    homer_findMotifsGenome_arg.set('')
    meme_chip_arg.set('')

    rearrange_entry_field_function()


def default_setting_button_function():
    setting_table_notification_string_var.set('Default values restored')
    setting_table_notification_label.config(fg = 'blue')
    setting_table_string_var.set(default_setting_table_file_full_path)
    read_setting_table_function(default_setting_table_file_full_path)

    rearrange_entry_field_function()


def setting_table_window_open():
    if '[MODIFIED]' in setting_table_string_var.get():
        pass
    else:
        if setting_table_loading_test_function() == False:
            return
        elif setting_table_loading_test_function() == True:
            pass

    global setting_table_window

    setting_table_window = tk.Toplevel()
    setting_table_window.grab_set()

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

    window_default_setting_button = tk.Button(setting_table_window, text = 'Restore defaults', command = lambda:default_setting_button_function(), state = tk.NORMAL, width = 22)
    CreateToolTip(window_default_setting_button, text = 'Click here if you want to restore\nto ChIP-AP pipeline default settings')
    window_default_setting_button.grid(sticky = "w", row = 60, column = 11, padx = 10, pady = (4,10), rowspan = 1, columnspan = 1)

    setting_table_window_close_button = tk.Button(setting_table_window, text = 'Accept and close', command = lambda:setting_table_window_close(), state = tk.NORMAL, width = 22)
    CreateToolTip(setting_table_window_close_button, text = 'Click here to accept the settings\nabove and close this window')
    setting_table_window_close_button.grid(sticky = "e", row = 60, column = 12, padx = 10, pady = (4,10), rowspan = 1, columnspan = 8)

    rearrange_entry_field_function()


def setting_table_window_close():
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

    if setting_changed == 1:
        setting_table_notification_string_var.set('Settings were manually modified by user')
        setting_table_notification_label.config(fg = 'orange2')
        if '[MODIFIED]' not in setting_table_string_var.get(): 
            setting_table_string_var.set('{}{}'.format(current_setting_table_full_path, '[MODIFIED]'))
    elif setting_changed == 0:
        pass

    setting_table_window.destroy()


setting_table_notification_label = tk.Label(root, textvariable = setting_table_notification_string_var, width = 70, padx = 10, pady = 4, fg = 'blue')
setting_table_notification_string_var.set("Currently using default settings table")
setting_table_notification_label.grid(row = 6, column = 11, padx = 10, pady = 4, columnspan = 9)

setting_table_button = tk.Button(root, text = 'Load setting table', command = lambda:setting_table_button_function(), state = tk.NORMAL, width = 15)
CreateToolTip(setting_table_button, text = 'Click here to browse and\nselect your setting table file')
setting_table_button.grid(sticky = "we", row = 7, column = 11, padx = 10, rowspan = 1, columnspan = 1)

read_setting_table_function(setting_table_string_var.get())
setting_table_entry = tk.Entry(root, textvariable = setting_table_string_var, width = 50, justify = tk.RIGHT)
setting_table_entry.grid(sticky = "we", row = 7, column = 12, padx = (0,10), columnspan = 8, ipady = 3)

default_setting_button = tk.Button(root, text = 'Default settings', command = lambda:default_setting_button_function(), state = tk.NORMAL, width = 15)
CreateToolTip(default_setting_button, text = 'Click here if you want to restore\nto ChIP-AP pipeline default settings')
default_setting_button.grid(sticky = "we", row = 8, column = 11, padx = 10, pady = (4,5), rowspan = 1, columnspan = 1)

setting_table_window_open_button = tk.Button(root, text = 'Manual customization >>', command = lambda:setting_table_window_open(), state = tk.NORMAL, width = 25)
CreateToolTip(setting_table_window_open_button, text = "Warning: proceed only when\nyou know what you are doing.\nOtherwise, leave all settings\nat their default values.\nInvalid values will\nlikely cause the whole\npipeline run to break.\nCheck GitHub documentation\nto learn more about\ncustom settings.")
setting_table_window_open_button.grid(sticky = "e", row = 8, column = 12, padx = 10, pady = (4,5), rowspan = 1, columnspan = 8)
    
ttk.Separator(root, style = "Line.TSeparator", orient = tk.HORIZONTAL).grid(row = 10, column = 10, columnspan = 11, sticky = "we")

########################################################################################################################

def chip_rep1_r1_button_function():
    chip_rep1_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a file')
    if chip_rep1_r1_button_input:
        chip_rep1_r1_string_var.set(chip_rep1_r1_button_input.name)
        current_dir_string_var.set(os.path.dirname(chip_rep1_r1_button_input.name))

def chip_rep1_r2_button_function():
    chip_rep1_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a file')
    if chip_rep1_r2_button_input:
        chip_rep1_r2_string_var.set(chip_rep1_r2_button_input.name)
        current_dir_string_var.set(os.path.dirname(chip_rep1_r2_button_input.name))

def chip_rep2_r1_button_function():
    chip_rep2_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a file')
    if chip_rep2_r1_button_input:
        chip_rep2_r1_string_var.set(chip_rep2_r1_button_input.name)
        current_dir_string_var.set(os.path.dirname(chip_rep2_r1_button_input.name))

def chip_rep2_r2_button_function():
    chip_rep2_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a file')
    if chip_rep2_r2_button_input:
        chip_rep2_r2_string_var.set(chip_rep2_r2_button_input.name)
        current_dir_string_var.set(os.path.dirname(chip_rep2_r2_button_input.name))

def chip_rep3_r1_button_function():
    chip_rep3_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a file')
    if chip_rep3_r1_button_input:
        chip_rep3_r1_string_var.set(chip_rep3_r1_button_input.name)
        current_dir_string_var.set(os.path.dirname(chip_rep3_r1_button_input.name))

def chip_rep3_r2_button_function():
    chip_rep3_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a file')
    if chip_rep3_r2_button_input:
        chip_rep3_r2_string_var.set(chip_rep3_r2_button_input.name)
        current_dir_string_var.set(os.path.dirname(chip_rep3_r2_button_input.name))

def chip_rep4_r1_button_function():
    chip_rep4_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a file')
    if chip_rep4_r1_button_input:
        chip_rep4_r1_string_var.set(chip_rep4_r1_button_input.name)
        current_dir_string_var.set(os.path.dirname(chip_rep4_r1_button_input.name))

def chip_rep4_r2_button_function():
    chip_rep4_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a file')
    if chip_rep4_r2_button_input:
        chip_rep4_r2_string_var.set(chip_rep4_r2_button_input.name)
        current_dir_string_var.set(os.path.dirname(chip_rep4_r2_button_input.name))

def chip_rep5_r1_button_function():
    chip_rep5_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a file')
    if chip_rep5_r1_button_input:
        chip_rep5_r1_string_var.set(chip_rep5_r1_button_input.name)
        current_dir_string_var.set(os.path.dirname(chip_rep5_r1_button_input.name))

def chip_rep5_r2_button_function():
    chip_rep5_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a file')
    if chip_rep5_r2_button_input:
        chip_rep5_r2_string_var.set(chip_rep5_r2_button_input.name)
        current_dir_string_var.set(os.path.dirname(chip_rep5_r2_button_input.name))

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

chip_sample_notification_label = tk.Label(root, textvariable = chip_sample_notification_string_var, width = 70, padx = 10, pady = 4)
chip_sample_notification_label.grid(sticky = "we", column = 1, row = 26, padx = 10, pady = 4, columnspan = 9)

ttk.Separator(root, style = "Line.TSeparator", orient = tk.HORIZONTAL).grid(row = 30, column = 1, columnspan = 11, sticky = "we")

########################################################################################################################

def ctrl_rep1_r1_button_function():
    ctrl_rep1_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a file')
    if ctrl_rep1_r1_button_input:
        ctrl_rep1_r1_string_var.set(ctrl_rep1_r1_button_input.name)
        current_dir_string_var.set(os.path.dirname(ctrl_rep1_r1_button_input.name))

def ctrl_rep1_r2_button_function():
    ctrl_rep1_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a file')
    if ctrl_rep1_r2_button_input:
        ctrl_rep1_r2_string_var.set(ctrl_rep1_r2_button_input.name)
        current_dir_string_var.set(os.path.dirname(ctrl_rep1_r2_button_input.name))

def ctrl_rep2_r1_button_function():
    ctrl_rep2_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a file')
    if ctrl_rep2_r1_button_input:
        ctrl_rep2_r1_string_var.set(ctrl_rep2_r1_button_input.name)
        current_dir_string_var.set(os.path.dirname(ctrl_rep2_r1_button_input.name))

def ctrl_rep2_r2_button_function():
    ctrl_rep2_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a file')
    if ctrl_rep2_r2_button_input:
        ctrl_rep2_r2_string_var.set(ctrl_rep2_r2_button_input.name)
        current_dir_string_var.set(os.path.dirname(ctrl_rep2_r2_button_input.name))

def ctrl_rep3_r1_button_function():
    ctrl_rep3_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a file')
    if ctrl_rep3_r1_button_input:
        ctrl_rep3_r1_string_var.set(ctrl_rep3_r1_button_input.name)
        current_dir_string_var.set(os.path.dirname(ctrl_rep3_r1_button_input.name))

def ctrl_rep3_r2_button_function():
    ctrl_rep3_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a file')
    if ctrl_rep3_r2_button_input:
        ctrl_rep3_r2_string_var.set(ctrl_rep3_r2_button_input.name)
        current_dir_string_var.set(os.path.dirname(ctrl_rep3_r2_button_input.name))

def ctrl_rep4_r1_button_function():
    ctrl_rep4_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a file')
    if ctrl_rep4_r1_button_input:
        ctrl_rep4_r1_string_var.set(ctrl_rep4_r1_button_input.name)
        current_dir_string_var.set(os.path.dirname(ctrl_rep4_r1_button_input.name))

def ctrl_rep4_r2_button_function():
    ctrl_rep4_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a file')
    if ctrl_rep4_r2_button_input:
        ctrl_rep4_r2_string_var.set(ctrl_rep4_r2_button_input.name)
        current_dir_string_var.set(os.path.dirname(ctrl_rep4_r2_button_input.name))

def ctrl_rep5_r1_button_function():
    ctrl_rep5_r1_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a file')
    if ctrl_rep5_r1_button_input:
        ctrl_rep5_r1_string_var.set(ctrl_rep5_r1_button_input.name)
        current_dir_string_var.set(os.path.dirname(ctrl_rep5_r1_button_input.name))

def ctrl_rep5_r2_button_function():
    ctrl_rep5_r2_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a file')
    if ctrl_rep5_r2_button_input:
        ctrl_rep5_r2_string_var.set(ctrl_rep5_r2_button_input.name)
        current_dir_string_var.set(os.path.dirname(ctrl_rep5_r2_button_input.name))

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

ctrl_sample_notification_label = tk.Label(root, textvariable = ctrl_sample_notification_string_var, width = 70, padx = 10, pady = 4)
ctrl_sample_notification_label.grid(sticky = "we", column = 11, row = 26, padx = 10, pady = 4, columnspan = 9)

ttk.Separator(root, style = "Line.TSeparator", orient = tk.HORIZONTAL).grid(row = 30, column = 11, columnspan = 10, sticky = "we")

########################################################################################################################

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

    rearrange_entry_field_function()


genome_ref_label = tk.Label(root, text = "Reference genome:")
genome_ref_label.grid(row = 35, column = 1, padx = 10, columnspan = 1, pady = (5,0))

genome_ref_drop_down = tk.OptionMenu(root, genome_ref_string_var, *genome_ref_options, command = lambda x = None : auto_assign_genome_folder_function())
CreateToolTip(genome_ref_drop_down, text = 'Select the genome assembly you want\nyour ChIP-seq reads to be aligned to.\nCurrently, ChIP-AP supports\nsix genome assemblies.\nOutside those supported by\nChIP-AP, you will need to\ngenerate the files yourself')
genome_ref_drop_down.config(takefocus = 1)
genome_ref_drop_down.grid(sticky = "we", row = 35, column = 2, columnspan = 4, pady = (5,0), padx = (0,10))

########################################################################################################################

def genome_folder_button_function():
    genome_folder_button_input = filedialog.askdirectory(initialdir = current_dir_string_var.get(), parent = root, title = 'Choose the Genome Directory')
    if genome_folder_button_input:
        genome_folder_string_var.set(genome_folder_button_input)
        current_dir_string_var.set(genome_folder_button_input)

    rearrange_entry_field_function()


genome_folder_button = tk.Button(root, text = 'Genome folder', command = lambda:genome_folder_button_function(), state = tk.DISABLED, width = 15)
CreateToolTip(genome_folder_button, text = 'Click here to browse and select\nthe directory containing your\ncustom genome reference')
genome_folder_button.grid(sticky = "we", row = 37, column = 1, padx = 10, columnspan = 1)

genome_folder_entry = tk.Entry(root, textvariable = genome_folder_string_var, width = 40, justify = tk.RIGHT, state = tk.DISABLED)
genome_folder_entry.grid(sticky = "we", row = 37, column = 2, columnspan = 4, ipady = 3, padx = (0,10))

########################################################################################################################

def known_motif_button_function():
    known_motif_button_input = filedialog.askopenfile(initialdir = current_dir_string_var.get(), parent = root, mode = 'rb', title = 'Choose a .motif file')
    if known_motif_button_input:
        known_motif_string_var.set(known_motif_button_input.name)
        current_dir_string_var.set(os.path.dirname(known_motif_button_input.name))

    rearrange_entry_field_function()


known_motif_button = tk.Button(root, text = 'Known motif file', command = lambda:known_motif_button_function(), state = tk.NORMAL, width = 15)
CreateToolTip(known_motif_button, text = 'Click here to browse and select your\n.motif file (in HOMER matrix format)')
known_motif_button.grid(sticky = "we", row = 39, column = 1, padx = 10, columnspan = 1)

known_motif_entry = tk.Entry(root, textvariable = known_motif_string_var, width = 40, justify = tk.RIGHT)
known_motif_entry.grid(sticky = "we", row = 39, column = 2, columnspan = 4, ipady = 3, padx = (0,10))

########################################################################################################################

def output_folder_button_function():
    output_folder_button_input = filedialog.askdirectory(initialdir = current_dir_string_var.get(), parent = root, title = 'Choose the Output Directory')
    if output_folder_button_input:
        output_folder_string_var.set(output_folder_button_input)
        current_dir_string_var.set(output_folder_button_input)

    rearrange_entry_field_function()


output_folder_button = tk.Button(root, text = 'Output save folder', command = lambda:output_folder_button_function(), state = tk.NORMAL, width = 15)
CreateToolTip(output_folder_button, text = 'Click here to browse and select\nthe directory to save the results\nof your ChIP-AP pipeline run.\nTo save in a new folder,\ntype in the entry field\nyour new folder name after\nthe output save directory:\nfull_path_to_output_save_directory/\nnew_folder_name')
output_folder_button.grid(sticky = "we", column = 1, row = 41, padx = 10, columnspan = 1)

output_folder_entry = tk.Entry(root, textvariable = output_folder_string_var, width = 40, justify = tk.RIGHT)
output_folder_entry.grid(sticky = "we", column = 2, row = 41, columnspan = 4, ipady = 3, padx = (0,10))

########################################################################################################################

setname_label = tk.Label(root, text = "Dataset name:", width = 15)
setname_label.grid(row = 43, column = 1, padx = 10, columnspan = 1, pady = (0,5))

setname_entry = tk.Entry(root, textvariable = setname_string_var, width = 40, justify = tk.LEFT)
CreateToolTip(setname_entry, text = 'Type in your folder name and prefix\nfor all the resulting output filenames')
setname_entry.grid(sticky = "we", row = 43, column = 2, columnspan = 4, ipady = 3, pady = (0,5), padx = (0,10))

########################################################################################################################

homer_motif_label = tk.Label(root, text = "HOMER motif enrichment:", anchor = "w")
homer_motif_label.grid(row = 35, column = 7, padx = 5, columnspan = 5, sticky = "sw", pady = (5,0))

homer_motif_drop_down = tk.OptionMenu(root, homer_motif_string_var, *homer_motif_options)
CreateToolTip(homer_motif_drop_down, text = 'Select the peak set(s) you want\nHOMER findMotifsGenome to perform\nmotif enrichment analysis on')
homer_motif_drop_down.config(takefocus = 1)
homer_motif_drop_down.grid(sticky = "we", row = 37, column = 7, padx = 5, columnspan = 5)

########################################################################################################################

meme_motif_label = tk.Label(root, text = "MEME motif enrichment:", anchor = "w")
meme_motif_label.grid(row = 39, column = 7, padx = 5, columnspan = 5, sticky = "sw")

meme_motif_drop_down = tk.OptionMenu(root, meme_motif_string_var, *meme_motif_options)
CreateToolTip(meme_motif_drop_down, text = 'Select the peak set(s) you\nwant meme-chip to perform\nmotif enrichment analysis on')
meme_motif_drop_down.config(takefocus = 1)
meme_motif_drop_down.grid(sticky = "we", row = 41, column = 7, padx = 5, columnspan = 5)

########################################################################################################################

fcmerge_checkbox = tk.Checkbutton(root, text = ' Merged fold enrichment analysis', variable = fcmerge_var, onvalue = 1, offvalue = 0)
CreateToolTip(fcmerge_checkbox, text = 'Check this box if you want\nthe fold enrichment analysis\nfrom all replicates combined as one.\nThis option will be ignored\nwhen there are unequal\nnumber of replicates between\nChIP and control samples')
fcmerge_checkbox.grid(sticky = "w", row = 35, column = 13, columnspan = 8, padx = (0,10), pady = (5,0))

########################################################################################################################

goann_checkbox = tk.Checkbutton(root, text = ' Annotate peaks with known gene ontology terms', variable = goann_var, onvalue = 1, offvalue = 0)
CreateToolTip(goann_checkbox, text = 'Check this box if you want\neach peak in the final peaks list\nto have gene ontology annotations')
goann_checkbox.grid(sticky = "w", row = 37, column = 13, columnspan = 8, padx = (0,10))

########################################################################################################################

pathann_checkbox = tk.Checkbutton(root, text = ' Annotate peaks with known pathway terms', variable = pathann_var, onvalue = 1, offvalue = 0)
CreateToolTip(pathann_checkbox, text = 'Check this box if you want\neach peak in the final peaks list\nto have pathway annotations')
pathann_checkbox.grid(sticky = "w", row = 39, column = 13, columnspan = 8, padx = (0,10))

########################################################################################################################

deltemp_checkbox = tk.Checkbutton(root, text = ' Delete large temporary files', variable = deltemp_var, onvalue = 1, offvalue = 0)
CreateToolTip(deltemp_checkbox, text = 'Check this box if you will not need\nthe large-sized intermediary files')
deltemp_checkbox.grid(sticky = "w", row = 41, column = 13, columnspan = 8, padx = (0,10))

########################################################################################################################

stdout_checkbox = tk.Checkbutton(root, text = ' Record standard outputs & errors', variable = stdout_var, onvalue = 1, offvalue = 0)
CreateToolTip(stdout_checkbox, text = 'Check this box if you want to save\npipeline standard outputs (channel 1>) and\nstandard errors (channel 2>) as text files')
stdout_checkbox.grid(sticky = "w", row = 43, column = 13, columnspan = 8, padx = (0,10), pady = (0,5))

########################################################################################################################

# stderr_checkbox = tk.Checkbutton(root, text = ' Record standard errors', variable = stderr_var, onvalue = 1, offvalue = 0)
# CreateToolTip(stderr_checkbox, text = 'Check this box if you want to save pipeline errors (channel 2>) in a text file')
# stderr_checkbox.grid(sticky = "w", row = 45, column = 11, columnspan = 10, padx = 10, pady = (5,5))

ttk.Separator(root, style = "Line.TSeparator", orient = tk.HORIZONTAL).grid(row = 46, column = 0, columnspan = 21, sticky = "we")

########################################################################################################################

cpu_count_label = tk.Label(root, text = "CPU cores to use:", width = 20, justify = tk.RIGHT)
cpu_count_label.grid(row = 62, column = 1, padx = 10, columnspan = 1, pady = 5)

cpu_count_entry = tk.Entry(root, width = 5, textvariable = cpu_count_string_var)
CreateToolTip(cpu_count_entry, text = 'Type in the number of CPU cores\nto be used by the pipeline')
cpu_count_entry.grid(sticky = "w", row = 62, column = 2, padx = (0,10), ipady = 3, pady = 5)

cpu_count_notification_label = tk.Label(root, textvariable = cpu_count_notification_string_var, width = 65, anchor = "w")
cpu_count_notification_label.grid(row = 62, column = 3, columnspan = 9, padx = 10)

########################################################################################################################

command_line_label = tk.Label(root, text = 'ChIP-AP command line:', width = 50)
command_line_label.grid(sticky = "we", column = 1, row = 47, padx = 10, columnspan = 19, pady = 5)

command_line_font = tkFont.Font(size = 12)

command_line_output_label = tk.Label(root, textvariable = command_line_output_string_var, relief = tk.SOLID, borderwidth = 1, bg = "white",
                                    width = 50, anchor = "w", justify = tk.LEFT, padx = 10, pady = 4, wraplength = 1450)

command_line_output_label['font'] = command_line_font
command_line_output_label.grid(sticky = "we", column = 1, row = 48, padx = 10, columnspan = 19, ipady = 4)

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
        fold_change_calculator_arg.get(),
        homer_findMotifsGenome_arg.get(),
        meme_chip_arg.get()]

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

    rearrange_entry_field_function()
    print_sample_table_function()
    print_setting_table_function()
    print_GS_command_line_function()
    execute_GS_command_line_function()
    cpu_count_notification_string_var.set("Scripts generated! You can close this window now")
    cpu_count_notification_label.config(fg = 'green')


def generate_scripts_and_run_function():
    if not os.path.exists(output_dir.get()):
        os.makedirs(output_dir.get())
    
    rearrange_entry_field_function()
    print_sample_table_function()
    print_setting_table_function()
    print_GSaR_command_line_function()
    cpu_count_notification_string_var.set("Pipeline started! You can close this window now")
    cpu_count_notification_label.config(fg = 'green')
    execute_GSaR_command_line_function()


generate_scripts_button = tk.Button(root, text = 'Generate scripts', command = lambda:generate_scripts_function(), state = tk.DISABLED, width = 18)
CreateToolTip(generate_scripts_button, text = 'Select this if you wish to\nrun the pipeline later\nby executing MASTER_script.sh\nwithin the output save folder')
generate_scripts_button.grid(row = 62, column = 12, columnspan = 9, sticky = "w", padx = 10, pady = 5)

generate_and_run_scripts_button = tk.Button(root, text = 'Generate and run scripts', command = lambda:generate_scripts_and_run_function(), state = tk.DISABLED, width = 25)
CreateToolTip(generate_and_run_scripts_button, text = 'Select this if you wish to\nrun the pipeline now\nNOTE: It may take up to\nseveral hours depending\non your system')
generate_and_run_scripts_button.grid(row = 62, column = 12, columnspan = 9, sticky = "e", padx = 10, pady = 5)

########################################################################################################################

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

register_sample_function()
check_required_input_function()
update_command_line_function()
rearrange_entry_field_function()
root.mainloop()
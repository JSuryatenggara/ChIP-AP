#!/usr/bin/env python3
#pyright: reportUnboundVariable=false


script_version = '1.1'

# Generates sequences for MEME-suite input based on fold-change-calculated peak list (ChIP-AP: dataset_name_all_peaks_calculated.tsv), featuring:
#   - Optional resizing of sequences length
#   - Optional resampling of number of sequences
#   - Optional masking of genomic regions with repeat sequences
#   - Optional generation of background control sequences
#   - Optional customization of background sequences GC% content
#   - Optional filtering to include only a subset of peaks called by specific peak caller(s)

# INPUT     - "dataset_name_all_peaks_calculated" - fold-change-calculated peak list
#           - Reference genome FASTA sequence
#           - Reference genome chromosome sizes

# OUTPUT    - "targets.fa" - Sequences of the peaks in the complete peak list for MEME-suite motif enrichment analysis 
#           - "background.fa" - Optional. Randomized background control sequences from the same genome for MEME-suite motif enrichment analysis 

# PATCH NOTES
#   Version 1.1     No longer needs hard-masked version of the reference genome in order to run --masked mode
#                   - Lowercase characters from the standard version of the reference genome (soft-masked) are now automatically converted to N


import pysam
import argparse
import pandas as pd
import numpy as np
import os
import subprocess
import csv
import random   
import math


subprocess.run('ulimit -n 2000', shell = True)


print('\nParsing command line flags and arguments...')

parser = argparse.ArgumentParser()

parser.add_argument('--input',
                    help = '<Required> Your fold-change-calculated peak list (ChIP-AP: dataset_name_all_peaks_calculated.tsv)',
                    required = True)

parser.add_argument('--fastadir', 
                    help = '<Required> Your reference genome FASTA file folder.', 
                    required = True)

parser.add_argument('--chrsizedir', 
                    help = '<Required> Your reference genome chromosome sizes file folder.', 
                    required = True)

parser.add_argument('--ref', 
                    help = '<Optional> Your sample organism reference genome build.', 
                    choices = ['hg19', 'hg38', 'mm9', 'mm10', 'dm6', 'sacCer3'],
                    required = True)

parser.add_argument('--filter', 
                    help = '<Optional> The peakset to process. Accepted values integers (1-4) to select all peaks with at least the given number of peak caller overlaps (4 for consensus peaks, 1 for all a.k.a. all peaks), or the name of the peak caller (case-sensitive) to select all peaks called by the given peak caller (e.g., MACS2, Genrich, etc).',
                    default = 1)

parser.add_argument('--masked',
                    help = '<Optional> Use the masked version of the reference genome FASTA (N for every base in repetitive elements instead of lowercase letter).',        
                    action = 'store_true')

parser.add_argument('--outputdir',
                    help = '<Required> Your desired output folder.',
                    required = True)

parser.add_argument('--background',
                    help = '<Optional> Generate FASTA files for background sequences.',
                    action = 'store_true')

parser.add_argument('--length',
                    help = '<Optional> Set a fixed length for your target and background sequences. When not set and background sequence are not to be generated, target sequences will have their initial varying length. When not set and background sequences are to be generated, both target and background sequences will have a fixed length equals to the mean length of all target sequences initial lengths.',
                    default = 'auto')

parser.add_argument('--gc',
                    help = '<Optional> Set a fixed GC content (within the range of +/- 0.05) of the generated background sequences. Accepted value is between 0 to 1. When not set, the GC content of the background sequences will be set to be the mean GC content of all target sequences +/- 0.05.',
                    default = 'auto')

parser.add_argument('--target_sampling',
                    help = '<Optional> Set a number of target sequences to be generated when an integer is given (e.g., 10000) as an argument. Set a multiplier (relative to the initial number of target sequences) of target sequences to be generated when a number with "x" prefix is given (e.g., x0.8) as an argument. When not set, the number of generated target sequences will be the same as the initial number of target sequences (equals to x1). Number of generated target sequences cannot be higher than initial (downsampling only).',
                    default = 'x1')

parser.add_argument('--background_sampling',
                    help = '<Optional> Set a number of background sequences to be generated when an integer is given (e.g., 10000) as an argument. Set a multiplier (relative to the number of generated target sequences) of background sequences to be generated when a number with "x" prefix is given (e.g., x2.5) as an argument. When not set, the number of generated background sequences will be the same as the number of target sequences (equals to x1).',
                    default = 'x1')

args = parser.parse_args()


# Get the full paths of the given files or directories
fasta_dir_full_path = os.path.abspath(args.fastadir)
chr_size_dir_full_path = os.path.abspath(args.chrsizedir)
input_full_path = os.path.abspath(args.input)
output_dir_full_path = os.path.abspath(args.outputdir)

fasta_full_path = os.path.abspath('{}/{}.fa'.format(fasta_dir_full_path, args.ref))
chr_size_full_path = os.path.abspath('{}/{}.chrom.sizes'.format(chr_size_dir_full_path, args.ref))


########################################################################################################################


try:
    peakset_filter = int(args.filter) # If --filter is number of peak caller overlaps
except:
    peakset_filter = str(args.filter) # If --filter is peak caller name


if 'x' in args.target_sampling: # If target sampling is a multiple of initial number
    try:
        target_sampling_multiplier = float(args.target_sampling.strip('x')) # Should be a valid multiplier value

    except:
        print('\n--target_sampling: Value error. Please key in an a valid value. See --help for more information') # Print error and exit if invalid value
        exit()

    if target_sampling_multiplier > 1:
        print('\n--target_sampling: Upsampling is not supported. Please choose a multiplier value of less than 1 (e.g., x0.2)') # Remind user to not upsample, then exit
        exit()
    


if 'x' not in args.target_sampling: # If target sampling is the desired actual number
    try:
        target_sampling_number = int(args.target_sampling) # Should be a valid integer value

    except:
        print('\n--target_sampling: Value error. Please key in an a valid value. See --help for more information') # Print error and exit if invalid value
        exit()



if 'x' in args.background_sampling: # If background sampling is a multiple of initial number
    try:
        background_sampling_multiplier = float(args.background_sampling.strip('x')) # Should be a valid multiplier value

    except:
        print('\n--background_sampling: Value error. Please key in an a valid value. See --help for more information') # Print error and exit if invalid value
        exit()



if 'x' not in args.background_sampling: # If background sampling is the desired actual number
    try:
        background_sampling_number = int(args.background_sampling) # Should be a valid integer value

    except:
        print('\n--background_sampling: Value error. Please key in an a valid value. See --help for more information') # Print error and exit if invalid value
        exit()



if args.length != 'auto': # If sequence length is set to a fixed value
    try:
        resized_target_length = int(args.length) # Should be a valid integer value
    except:
        print('\n--length: Value error. Please key in an a valid value. See --help for more information') # Print error and exit if invalid value
        exit()



if args.gc != 'auto': # If background sequences GC content is set to a fixed value +/- 0.05
    try:
        user_value_GC_content = float(args.gc) # Should be a valid multiplier value
    except:
        print('\n--gc: Value error. Please key in a valid value. See --help for more information') # Print error and exit if invalid value
        exit()

    if user_value_GC_content > 1:
        print('\n--gc: GC content value of higher than 1 (100%) is not possible. Please choose a fractional value of less than 1 (e.g., x0.5)') # "Seriously?", then exit
        exit()


########################################################################################################################


# Define pysam object for sequence reading of the reference genome FASTA file at given coordinates
genome_pysam = pysam.Fastafile(fasta_full_path)

# Read the chromosome sizes file of the respective reference genome (sets the upper limits for random generation of background sequences coordinates)
chr_size_df = pd.read_csv(chr_size_full_path, delimiter = '\t', header = None)
chr_size_dict = dict(chr_size_df.values)

# Read the peak list (dataset_name_all_peaks_calculated.tsv)
input_df = pd.read_csv(input_full_path, delimiter = '\t')
input_array = input_df.values.tolist()
input_header = input_df.columns.tolist()


# Prepare lists for the initial target sequences
filtered_input_array = []
target_info_list = []
target_length_list = []
center_column_number_list = []


print('\nReading the peak list to generate initial target sequences...')

# Parse the given peak list and find the relevant columns' index number by header values
for input_header_counter in range(len(input_header)):
    if input_header[input_header_counter] == 'Chr':
        chr_column_number = input_header_counter

    if input_header[input_header_counter] == 'Start':
        start_column_number = input_header_counter

    if input_header[input_header_counter] == 'End':
        end_column_number = input_header_counter

    if 'Peak Center' in input_header[input_header_counter]: # Parses multiple weighted peak center columns, if exist
        center_column_number_list.append(input_header_counter)

    if input_header[input_header_counter] == 'Peak Caller Combination':
        peak_caller_combination_column_number = input_header_counter

    if input_header[input_header_counter] == 'Peak Caller Overlaps':
        peak_caller_overlap_column_number = input_header_counter


# Set the number of replicates based on how many weighted peak center columns in the peak list if sequences are to be resized 
#   (Sequence padding = multi-replicated weighted peak centers = multi-replicated motif enrichment analysis)
if args.background or args.length != 'auto': # Resizing scenario
    replicate_number = len(center_column_number_list)

# Set the number of replicates to 1 if target sequences are kept at their initial length
#   (No sequence padding = multi-replicated weighted peak centers irrelevant)
else: # No resizing scenario
    replicate_number = 1


print('\nSelecting peaks based on user-determined filter...')

if type(peakset_filter) == int: # If --filter is number of peak caller overlaps
    for input_array_row in input_array:
        if input_array_row[peak_caller_overlap_column_number] >= peakset_filter: # Get only peaks with number of overlaps equal or greater than the value of the filter
            filtered_input_array.append(input_array_row)

if type(peakset_filter) == str: # If --filter is peak caller name
    for input_array_row in input_array:
        if peakset_filter in input_array_row[peak_caller_combination_column_number]: # Get only peaks containing a substring equal to the value of the filter
            filtered_input_array.append(input_array_row)

if len(filtered_input_array) == 0: # If filtering fails and ends up selecting no peak at all
    print('\nSelected peak set does not contain any peak. Please re-select your peak set with --filter or re-check the input peak list given with --input') # "Check again", then exit
    exit()


########################################################################################################################


print('\nReading the selected peaks to generate the initial target sequences...')

# Reading the filtered peak list to generate the initial target sequences
for filtered_input_array_row in filtered_input_array:
    chr = filtered_input_array_row[chr_column_number]
    start = int(filtered_input_array_row[start_column_number])
    end = int(filtered_input_array_row[end_column_number])
    center_list = filtered_input_array_row[(np.min(np.array(center_column_number_list))) : (np.max(np.array(center_column_number_list)) + 1)] # Contains multiple values if exist.

    sequence = genome_pysam.fetch(chr, start, end) # Use the pysam object defined in the beginning to rapidly access and read the sequence
    
    if args.masked: # If --masked flag is used
        if not sequence.isupper(): # If the current pulled sequence has at least one lowercase character
            masked_sequence = [] # Prepare for the conversion into hard-masked sequence

            for base in sequence: # Read the sequence one by one from the left
                if base.isupper(): # If the current base is in uppercase
                    masked_sequence.append(base) # No change, write it as is in the hard-masked sequence
                elif base.islower(): # If the current base is in lowercase
                    masked_sequence.append('N') # Mask it, write is as N in the hard-masked sequence
            
            sequence =  ''.join(masked_sequence) # Join all the base characters in the list into a string

    info = [chr, start, end, sequence, center_list] # Info refers to complete information on a peak

    target_info_list.append(info) # Append each peak info to a list

    length = len(sequence) # Get the number of characters (bases) in the sequence returned by pysam

    target_length_list.append(length) # Append each peak lenght to a list (to calculate mean below)


target_sequence_number = len(target_info_list) # Get the initial number of target sequences
print('\nInitial number of target sequences: {}'.format(target_sequence_number)) # Let user know the initial number of target sequences 

mean_target_length = int(np.mean(np.array(target_length_list))) # Get the mean initial length of target sequences 
print('\nMean length of target sequences: {}'.format(mean_target_length)) # Let user know the mean initial length of target sequences 


########################################################################################################################


print('\nSampling the initial target sequences...')

if 'x' in args.target_sampling: # If background sampling is a multiple of initial number, let the user know the effects of given target sampling value
    sampled_target_sequence_number = int(target_sequence_number * target_sampling_multiplier)

    if target_sampling_multiplier == 1:
        print('\n--target_sampling {}. Sampling multiplier is equal to 1'.format(args.target_sampling))
        print('No downsampling performed. Proceeding with the initial number of target sequences: {}'.format(sampled_target_sequence_number))

    if target_sampling_multiplier < 1:
        print('\n--target_sampling {}: Number of target sequences has been downsampled by a factor of {} to {}'.format(args.target_sampling, target_sampling_multiplier, sampled_target_sequence_number))


if 'x' not in args.target_sampling: # If background sampling is the desired actual number, let the user know the effects of given target sampling value
    if target_sampling_number >= target_sequence_number:
        target_sampling_number = target_sequence_number
        sampled_target_sequence_number = target_sampling_number
        print('\n--target_sampling {}. Sampling number is larger than or equal to the initial number of target sequences'.format(args.target_sampling))
        print('No downsampling performed. Proceeding with the initial number of target sequences: {}'.format(sampled_target_sequence_number))

    if target_sampling_number < target_sequence_number:
        sampled_target_sequence_number = target_sampling_number
        print('\n--target_sampling {}. Number of target sequences has been downsampled to {}'.format(args.target_sampling, sampled_target_sequence_number))


sampled_target_info_list = random.sample(target_info_list, sampled_target_sequence_number) # Randomly sample a number of peaks from the initial target sequences


########################################################################################################################


print('\nSetting the length for the generated target sequences...')

resized_target_info_list = [] # Prepare list for the resized target sequences

if args.length == 'auto' and not args.background: # No resizing scenario 
    print('\n--length auto: auto length mode has been chosen, and background sequences are NOT to be generated')
    print('Target sequences are maintained at their respective individual sequence length') # Let user know that the target sequences are not resized and at varying lengths

    # Assign every variable as list here regardless of having only a single value in order to make them all iterable, just like the results of resizing scenario
    for sampled_target_info in sampled_target_info_list:
        resized_chr_list = [sampled_target_info[0]]
        resized_start_list = [sampled_target_info[1]]
        resized_end_list = [sampled_target_info[2]]
        resized_sequence_list = [sampled_target_info[3]]
        resized_center_list = sampled_target_info[4]

        resized_target_info_list.append([resized_chr_list, resized_start_list, resized_end_list, resized_sequence_list, resized_center_list]) # Stack this peak info


if args.length == 'auto' and args.background: # Automatic resizing scenario 
    resized_target_length = mean_target_length
    print('\n--length auto: auto length mode has been chosen, and background sequences are to be generated')
    print('The lengths of all target sequences are adjusted to equal to their mean length') # Let user know that all target sequences are resized into a uniform length

if args.length != 'auto': # User-determined resizing scenario 
    resized_target_length = int(args.length)
    print('\n--length {}: All the resulting target and background (if generated) sequences are set to be {} bp'.format(args.length, args.length))


if (args.length == 'auto' and args.background) or args.length != 'auto': # Resizing scenario 
    print('\n{} replicates detected based on the number of peak list weighted peak center columns'.format(replicate_number)) # Let user confirm the number of replicates
    print('\nResizing the sampled target sequences...')

    # Assign every variable here, with only the weighted peak center as a list variable
    for sampled_target_info in sampled_target_info_list:
        chr = sampled_target_info[0]
        start = sampled_target_info[1]
        end = sampled_target_info[2]
        sequence = sampled_target_info[3]
        center_list = sampled_target_info[4]

        # Prepare a list for every variable regardless of having only a single replicate in order to make them all iterable regardless of replicate number
        resized_chr_list = []
        resized_start_list = []
        resized_end_list = []
        resized_sequence_list = []
        resized_center_list = []


        for center in center_list: # Replicates are defined by the number of items in the weighted peak center list

            resized_chr_list.append(chr)
            
            resized_start = center - int(math.ceil(resized_target_length / 2)) # Define start coordinate from the weighted peak center
            resized_end = resized_start + resized_target_length # Define end coordinate from the start coordinate

            if resized_start < 0: # If the start coordinate is below zero (probably will not happen, but just in case)
                resized_end += (0 - resized_start) # Slide the end coordinate to an acceptable range
                resized_start = 0 # Slide the start coordinate to an acceptable range

                resized_start_list.append(resized_start) # Append the start coordinate of this replicate
                resized_end_list.append(resized_end) # Append the end coordinate of this replicate

            if resized_end > chr_size_dict[chr]: # If the end coordinate beyond the chromosome size
                resized_start -= (resized_end - chr_size_dict[chr]) # Slide the start coordinate to an acceptable range
                resized_end = chr_size_dict[chr] # Slide the end coordinate to an acceptable range
                
                resized_start_list.append(resized_start) # Append the start coordinate of this replicate
                resized_end_list.append(resized_end) # Append the end coordinate of this replicate
        
            else: # If both coordinates are within acceptable range
                resized_start_list.append(resized_start) # Append the start coordinate of this replicate
                resized_end_list.append(resized_end) # Append the end coordinate of this replicate

            resized_sequence = genome_pysam.fetch(chr, resized_start, resized_end) # Use the pysam object defined in the beginning to rapidly access and read the sequence

            if args.masked:
                if not resized_sequence.isupper():
                    masked_sequence = []
                    
                    for base in resized_sequence:
                        if base.isupper():
                            masked_sequence.append(base)
                        elif base.islower():
                            masked_sequence.append('N')
                    
                    resized_sequence =  ''.join(masked_sequence)

            resized_sequence_list.append(resized_sequence) # Append the sequence of this replicate returned by pysam

            resized_center_list = center_list # Append the weighted peak center coordinate of this replicate (just in case)

        resized_target_info_list.append([resized_chr_list, resized_start_list, resized_end_list, resized_sequence_list, resized_center_list]) # Stack this peak info


# Prepare lists for the recorded sequence composition stats from each replicate
resized_target_total_GC_content_list = []
resized_target_total_all_N_instances_list = []
resized_target_total_half_N_instances_list = []


for replicate_counter in range(replicate_number): # Process one by one replicate
    
    target_output_list = [] # Every item appended into this output list represents one row in the generated target sequences FASTA file

    # Set every sequence composition stats to zero before beginning list iteration
    resized_target_total_GC_instances = 0 
    resized_target_total_AT_instances = 0
    resized_target_total_all_N_instances = 0
    resized_target_total_half_N_instances = 0


    for resized_target_info in resized_target_info_list: # Now process one by one row (peak)

        resized_sequence = resized_target_info[3][replicate_counter] # Asses the sequence column of the current peak of the current replicate to count the bases
        
        resized_target_total_GC_instances += (resized_sequence.count('A') + resized_sequence.count('T') + resized_sequence.count('a') + resized_sequence.count('t')) # Counter + 1 if G or C
        resized_target_total_AT_instances += (resized_sequence.count('G') + resized_sequence.count('C') + resized_sequence.count('g') + resized_sequence.count('c')) # Counter + 1 if A or T
            
        if resized_sequence.count('N') == len(resized_sequence):
            resized_target_total_all_N_instances += 1 # Counter + 1 if the sequence consist of 100% repetitive elements (N). Only relevant when using masked reference genome.

        if resized_sequence.count('N') > (0.5 * len(resized_sequence)) and resized_sequence.count('N') < len(resized_sequence):
            resized_target_total_half_N_instances += 1 # Counter + 1 if 50% < % of N < 100% in the sequence. Only relevant when using masked reference genome.

        # Register all variables as lists here regardless of having only a single replicate in order to make them all iterable regardless of replicate number
        chr_list = resized_target_info[0]
        start_list = resized_target_info[1]
        end_list = resized_target_info[2]
        sequence_list = resized_target_info[3]

        # Generic formatting of peak ID. Used later as a header for each sequence in the generated FASTA files
        peak_ID = '>{}:{}-{}'.format(chr_list[replicate_counter], start_list[replicate_counter], end_list[replicate_counter])

        target_output_list.append(peak_ID) # Append the peak ID first to the output list so that it will be written as the header of each respective sequence

        chunk_list = [sequence_list[replicate_counter][i:i+100] for i in range(0, len(sequence_list[replicate_counter]), 100)] # Slice the sequence into chunks of 100 bases

        for chunk in chunk_list:
            target_output_list.append(chunk) # Write the sequence below the peak ID, chunk by chunk, before moving on to the next peak's sequence


    # Calculate GC content of the target sequences from all peaks in this replicate
    resized_target_total_GC_content = resized_target_total_GC_instances / (resized_target_total_GC_instances + resized_target_total_AT_instances)

    # Keep all these value below for the subsequent background sequences processing
    resized_target_total_GC_content_list.append(resized_target_total_GC_content) # Append the GC content of this replicate calculated just above
    resized_target_total_all_N_instances_list.append(resized_target_total_all_N_instances) # Append the number of sequences with 100% N in this replicate
    resized_target_total_half_N_instances_list.append(resized_target_total_half_N_instances) # Append the number of sequences with 50% < N < 100% in this replicate


    if replicate_number == 1: # When there are only one replicate
        target_file_name = '{}/target.fa'.format(output_dir_full_path) # Simply name the output file target.fa

    if replicate_number > 1: # When there are multiple replicates
        target_file_name = '{}/target_rep{}.fa'.format(output_dir_full_path, (replicate_counter + 1)) # Put the replicate number behind the output file name

    # Write down the target output list under their designated file name, replicate-wise
    target_file = open(target_file_name, 'w')
    target_writer = csv.writer(target_file, delimiter='\n')
    target_writer.writerow(target_output_list)
    target_file.close()


    # Let the user know the stats of the generated target sequences, replicate-wise
    print('\nTarget sequences (replicate {}) has been generated: {}'.format((replicate_counter + 1), target_file_name))
    print('Target sequences (replicate {}) stats:'.format(replicate_counter + 1))
    print('Number of sequences: {}'.format(sampled_target_sequence_number))
    print('Length of sequences: {}'.format(resized_target_length))
    print('Total G or C characters: {}'.format(resized_target_total_GC_instances))
    print('Total A or T characters: {}'.format(resized_target_total_AT_instances))
    print('Overall GC content: {}'.format(round(resized_target_total_GC_content, 4)))
    # print('Number of all-N sequences: {}'.format(resized_target_total_all_N_instances))
    # print('Number of half-N sequences: {}'.format(resized_target_total_half_N_instances))


########################################################################################################################


print('\nSetting the number of background sequences...')

if args.background:

    if 'x' in args.background_sampling:
        background_sequence_number = int(sampled_target_sequence_number * background_sampling_multiplier)
        
        if background_sampling_multiplier > 1:
            print('\n--background_sampling {}. Number of target sequences has been upsampled by a factor of {} to {}'.format(args.background_sampling, background_sampling_multiplier, background_sequence_number))

        if background_sampling_multiplier == 1:
            print('\n--background_sampling {}.  Sampling value is larger than or equal to the number of target sequences'.format(args.background_sampling))
            print('No resampling performed. Proceeding with the number of target sequences: {}'.format(background_sequence_number))

        if background_sampling_multiplier < 1:
            print('\n--background_sampling {}. Number of target sequences has been downsampled by a factor of {} to {}'.format(args.background_sampling, background_sampling_multiplier, background_sequence_number))

    if 'x' not in args.background_sampling:
        background_sequence_number = background_sampling_number
        print('\n--background_sampling {}: Number of target sequences has been resampled to {}'.format(args.background_sampling, background_sequence_number))


    print('\nSetting the length of background sequences (always the same as resized target sequences length)...') # Duh, but maybe helps user understand better

    background_length = resized_target_length


    for replicate_counter in range(replicate_number): # Process one by one replicate

        # Prepare empty lists for the generated background info (may be necessary for module expansion) and output (to write the to the generated FASTA file(s))
        # Replicate wise, gets reset to empty in the beginning of each replicate
        background_info_list = []
        background_output_list = []

        # Set every sequence composition stats to zero before beginning list iteration
        # Replicate wise, gets reset to empty in the beginning of each replicate
        background_total_GC_instances = 0
        background_total_AT_instances = 0
        background_total_GC_content = 0
        background_total_all_N_instances = 0
        background_total_half_N_instances = 0   


        print('\nSetting the GC content of the background sequences (replicate {})...'.format(replicate_counter + 1))

        if args.gc == 'auto': # When GC content of background sequences follows the overall GC content of target sequences
            GC_content_lower_tolerance = resized_target_total_GC_content_list[replicate_counter] - 0.05 # Set the lowest accepted GC content for every generated sequence
            GC_content_upper_tolerance = resized_target_total_GC_content_list[replicate_counter] + 0.05 # Set the highest accepted GC content for every generated sequence
            print('\n--gc {}: GC content of the background sequences is set to be equal to the GC content of the target sequences ({}) +/- 0.05'.format(args.gc, resized_target_total_GC_content_list[replicate_counter]))

        elif args.gc != 'auto': # When GC content of background sequences is set to user-determined value
            GC_content_lower_tolerance = user_value_GC_content - 0.05 # Set the lowest accepted GC content for every generated sequence
            GC_content_upper_tolerance = user_value_GC_content + 0.05 # Set the highest accepted GC content for every generated sequence
            print('\n--gc {}: GC content of the background sequences is set to be equal to the inputted user determined value ({}) +/- 0.05'.format(args.gc, user_value_GC_content))

        resized_target_total_all_N_instances = resized_target_total_all_N_instances_list[replicate_counter] # Get the number of 100% N target sequences of current replicate
        # Sets the limit so the number of generated 100% N background sequences in current replicate does not exceed those in target sequences
        background_total_all_N_instances_limit = int(resized_target_total_all_N_instances * (background_sequence_number / sampled_target_sequence_number))

        resized_target_total_half_N_instances = resized_target_total_half_N_instances_list[replicate_counter] # Get the number of 50% < % of N < 100% target sequences of current replicate
        # Sets the limit so the number of generated 50% < % of N < 100% background sequences in current replicate does not exceed those in target sequences
        background_total_half_N_instances_limit = int(resized_target_total_half_N_instances * (background_sequence_number / sampled_target_sequence_number))


        print('\nGenerating background sequences for target sequences replicate {}...'.format(replicate_counter + 1))

        for sequence_counter in range(background_sequence_number): # Keep processing until the set number of background sequences is reached
            while True: # Keep looping until acceptable sequence composition is fulfilled. If the process does not meet a break, the loop is restarted to look for a new candidate.
                background_chr = random.choice([*chr_size_dict]) # Get sequence from a randomized chromosome (based on available chromosomes in the given chromosome sizes file)
                background_start = random.randint(0, (chr_size_dict[background_chr] - background_length)) # Set a randomized sequence start coordinate within acceptable range
                background_end = background_start + background_length # Set an end coordinate based on the randomized start coordinate
                background_sequence = genome_pysam.fetch(background_chr, background_start, background_end) # Get the sequence based on the randomized chr, start, and end coordinate above
                
                if args.masked:
                    if not background_sequence.isupper():
                        masked_sequence = []

                        for base in background_sequence:
                            if base.isupper():
                                masked_sequence.append(base)
                            elif base.islower():
                                masked_sequence.append('N')

                        background_sequence = ''.join(masked_sequence)

                # This sequence's G or C counter + 1 if G or C
                background_GC_instances = (background_sequence.count('A') + background_sequence.count('T') + background_sequence.count('a') + background_sequence.count('t')) 

                # This sequence's A or T counter + 1 if A or T
                background_AT_instances = (background_sequence.count('G') + background_sequence.count('C') + background_sequence.count('g') + background_sequence.count('c'))

                if (background_GC_instances + background_AT_instances) > 0: # If generated sequence is not 100% N...
                    background_GC_content = background_GC_instances / (background_GC_instances + background_AT_instances) # Calculate that sequence's GC content

                    # AND if the sequence composition is within the set acceptable GC content, it passes the test...
                    if (background_GC_content >= GC_content_lower_tolerance and background_GC_content <= GC_content_upper_tolerance):
                        background_total_GC_instances += background_GC_instances # If it passes, add the number of G or C bases to this replicate's counter
                        background_total_AT_instances += background_AT_instances # If it passes, add the number of A or T bases to this replicate's counter

                        # Unless the sequence composition is 50% < % of N < 100%
                        if background_sequence.count('N') > (0.5 * background_length) and background_sequence.count('N') < background_length:
                            
                            # Then the sequence only passes half N sequences limit has not been reached yet...
                            if background_total_half_N_instances < background_total_half_N_instances_limit:
                                background_total_half_N_instances += 1 # If it passes, count it towards the accepted half N sequences limit 
                            
                                break # The sequence passed all tests and will be generated as a background sequence

                        else: # However, if the sequence composition has less than 50% N...

                            break  # The sequence passed all tests and will be generated as a background sequence
                
                # Otherwise, if the sequence contains no A, C, T, or G (most likely because the randomized coordinate hits a repetitive element entirely a.k.a. 100% N)...
                elif (background_GC_instances + background_AT_instances) == 0:
                    
                    # The sequence only passes all N sequences limit has not been reached yet...
                    if background_total_all_N_instances < background_total_all_N_instances_limit: 

                        if background_sequence.count('N') == background_length: # Just to double check if 100% N is really the case
                            background_total_all_N_instances += 1 # If it passes, count it towards the accepted all N sequences limit 

                        break # The sequence passed all tests and will be generated as a background sequence


            # When the loop is broken (current background sequence candidate is accepted), all the processes below are carried out 

            background_info_list.append([background_chr, background_start, background_end, background_sequence]) # Append the accepted candidate info into the list


            # Generic formatting of peak ID. Used later as a header for each sequence in the generated FASTA files
            peak_ID = '>{}:{}-{}'.format(background_chr, background_start, background_end)
            background_output_list.append(peak_ID) # Append the peak ID first to the output list so that it will be written as the header of each respective sequence

            chunk_list = [background_sequence[i:i+100] for i in range(0, len(background_sequence), 100)] # Slice the sequence into chunks of 100 bases

            for chunk in chunk_list:
                background_output_list.append(chunk) # Write the sequence below the peak ID, chunk by chunk, before moving on to the next peak's sequence


        # Calculate GC content of the background sequences from all peaks in this replicate
        background_total_GC_content = background_total_GC_instances / (background_total_GC_instances + background_total_AT_instances)


        if replicate_number == 1: # When there are only one replicate
            background_file_name = '{}/background.fa'.format(output_dir_full_path) # Simply name the output file background.fa

        if replicate_number > 1: # When there are multiple replicates
            background_file_name = '{}/background_rep{}.fa'.format(output_dir_full_path, (replicate_counter + 1)) # Put the replicate number behind the output file name

        # Write down the target output list under their designated file name, replicate-wise
        background_file = open(background_file_name, 'w')
        background_writer = csv.writer(background_file, delimiter='\n')
        background_writer.writerow(background_output_list)
        background_file.close()


        # Let the user know the stats of the generated background sequences, replicate-wise
        print('\nBackground sequences (replicate {}) has been generated: {}'.format((replicate_counter + 1), background_file_name))
        print('Background sequences (replicate {}) stats:'.format(replicate_counter + 1))
        print('Number of sequences: {}'.format(background_sequence_number))
        print('Length of sequences: {}'.format(background_length))
        print('Total G or C characters: {}'.format(background_total_GC_instances))
        print('Total A or T characters: {}'.format(background_total_AT_instances))
        print('Individual background sequence GC content: {}'.format(round(background_total_GC_content, 4)))
        # print('Number of all-N sequences: {}'.format(background_total_all_N_instances))
        # print('Number of half-N sequences: {}'.format(background_total_half_N_instances))
        print('')
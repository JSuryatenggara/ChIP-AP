#!/usr/bin/env python3
#pyright: reportUnboundVariable=false


script_version = '1.0'


# Processes the peak list file "dataset_name_all_peaks_calculated" by:
#   - Adding the -log10 IDR values of the peak entries between two compared peak lists
#       - Specifically for implementation in ChIP-AP:
#           - Multiple peak IDR values of pairwise union peak list vs individual peak callers list (MACS/GEM/SICER2/HOMER/Genrich)
#           - Number of columns as many as individual peak callers list being compared to the union peak list (columns position: after GC%) 
#           - These columns will only be present in the output peak list directly resulting from this script run
#           - These columns will be removed from any subsequent resulting peak lists of downstream process(es)
#   - Adding the sum of -log10 IDR values (column position: after Number of Motifs).
#   - Adding the converted IDR values (column position: after Number of Motifs).


# INPUT     - "dataset_name_all_peaks_calculated" - gene annotated, fold change calculated, complete peaks list 
#           - Multiple peak IDR-calculated peak lists of pairwise union peak list vs individual peak callers list (MACS/GEM/SICER2/HOMER/Genrich)

# OUTPUT    - "dataset_name_all_peaks_calculated" - gene annotated, fold change calculated, IDR-integrated, complete peaks list.
#                                                 - saved under the same file name as the input peak list, replacing it.



print('Importing required modules')
import pandas as pd
import subprocess
import multiprocessing
import os
import argparse


print('Setting argument parser')

parser = argparse.ArgumentParser()

parser.add_argument('--input_tsv', 
                    help = '<Required> Input annotated peak list file, output of ChIP-AP fold_change_calculator.py script, .tsv extension.', 
                    required = True)

parser.add_argument('--idr_tsv', 
                    nargs = '+', help = '<Required> Input IDR-calculated peak list file(s), output of IDR script, .tsv extension, separated by space.', 
                    required = True)

parser.add_argument('--output_tsv', 
                    help = '<Required> Output IDR-integrated annotated peak list file, .tsv extension', 
                    required = True)

parser.add_argument('--thread', 
                    help = '<Optional> Maximum number of processes to use. Default is half the maximum available.', 
                    type = int, 
                    choices = range(1, (multiprocessing.cpu_count() + 1), 1),
                    metavar = "[1-{}]".format(multiprocessing.cpu_count()),
                    default = int(multiprocessing.cpu_count() / 2))

args = parser.parse_args()

subprocess.run('ulimit -n 2000', shell = True)


print('Parsing arguments')

cpu_count = args.thread

input_tsv_full_path = os.path.abspath(args.input_tsv) # ChIP-AP: dataset_name_all_peaks_calculated.tsv (output of fold_change_calculator.py)
output_tsv_full_path = os.path.abspath(args.output_tsv) # ChIP-AP: dataset_name_all_peaks_calculated.tsv (same name, peaks are now with IDR values)

idr_tsv_full_path_list = []

idr_tsv_list = args.idr_tsv # List of IDR module outputs from running dataset_name_all_peaks_calculated.tsv (union peak set) against individual peak caller sets (.narrowPeak or .broadPeak)

# Listing all the full paths to all the IDR output files
for idr_tsv in idr_tsv_list:
    idr_tsv_full_path = os.path.abspath(idr_tsv)
    idr_tsv_full_path_list.append(idr_tsv_full_path) 


print('Opening file: {}'.format(input_tsv_full_path)) # Reading the input fold change calculated peak list file
input_peak_df = pd.read_csv(input_tsv_full_path, delimiter = '\t')
input_peak_df.sort_values(by = ['Chr', 'Start'], inplace = True) # Sorting the peaks list based on chromosomal position
input_peak_array = input_peak_df.values.tolist()
input_peak_header = input_peak_df.columns.tolist()

print('Parsing header names for columns containing target data')
for input_peak_header_counter in range(len(input_peak_header)):

    # Parsing for the column containing the chromosome in which the peak is located in the peak list
    if input_peak_header[input_peak_header_counter] == 'Chr':
        chr_column_number = input_peak_header_counter

    # Parsing for the column containing the starting base coordinate of the peak in the peak list
    if input_peak_header[input_peak_header_counter] == 'Start':
        start_column_number = input_peak_header_counter

    # Parsing for the column containing the ending base coordinate of the peak in the peak list
    if input_peak_header[input_peak_header_counter] == 'End':
        end_column_number = input_peak_header_counter

    # Parsing for the column containing number of motifs that overlaps the peaks in the peak list
    #   Would be all zero if -motif flag was not used during annotatePeaks.pl run
    if input_peak_header[input_peak_header_counter] == 'Number of Motifs':
        motif_count_column_number = input_peak_header_counter


peak_idr_array_list = []

for idr_tsv_full_path in idr_tsv_full_path_list:
    print('Opening file: {}'.format(idr_tsv_full_path)) # Reading the input IDR output peak list files
    peak_idr_df = pd.read_csv(idr_tsv_full_path, delimiter = '\t', header = None)
    peak_idr_df.sort_values(by = [0, 1], inplace = True) # Sorting the peaks list based on chromosomal position
    peak_idr_array = peak_idr_df.values.tolist()
    peak_idr_array_list.append(peak_idr_array)

    output_peak_array = []


# Function to assign the correct -log10IDR values from the input IDR output peak list files into the correct peak entries
def peak_IDR_integration_function(input_peak_array_row):

    idr_value_list = []

    for idr_counter in range(len(peak_idr_array_list)): # Do this one by one for every input IDR output peak list files
    
        for peak_idr_array_row in peak_idr_array_list[idr_counter]: # Do this one by one for every peak entry in the input fold change calculated peak list file
            
            # If the peak coordinates match
            if input_peak_array_row[chr_column_number] == peak_idr_array_row[0]:
                if input_peak_array_row[start_column_number] == peak_idr_array_row[1]:
                    if input_peak_array_row[end_column_number] == peak_idr_array_row[2]:
                        
                        # Do this if the current input IDR output peak list file has .narrowPeak extension
                        if idr_tsv_full_path_list[idr_counter].endswith('.narrowPeak'):
                            idr_value = round(float(peak_idr_array_row[11]), 8)
                            idr_value_list.append(idr_value)
                            break

                        else: # If the current input IDR output peak list file does not have .narrowPeak extension
                            
                            # Do the same as above if the the number of columns match that of a .narrowPeak formatted file
                            if len(peak_idr_array_row) == 20:
                                idr_value = round(float(peak_idr_array_row[11]), 8)
                                idr_value_list.append(idr_value)
                                break
                            
                            else:
                                pass # Do nothing if conditions are not satisfied

                        # Do this if the current input IDR output peak list file has .broadPeak extension
                        if idr_tsv_full_path_list[idr_counter].endswith('.broadPeak'):
                            idr_value = round(float(peak_idr_array_row[10]), 8)
                            idr_value_list.append(idr_value)
                            break

                        else: # If the current input IDR output peak list file does not have .broadPeak extension
                            
                            # Do the same as above if the the number of columns match that of a .broadPeak formatted file
                            if len(peak_idr_array_row) == 17:
                                idr_value = round(float(peak_idr_array_row[10]), 8)
                                idr_value_list.append(idr_value)
                                break

                            else:
                                pass # Do nothing if conditions are not satisfied


    output_peak_array_row = input_peak_array_row
    
    summed_idr_value = round(sum(idr_value_list), 8) # Sum all -log10IDR values

    converted_idr_value = round(10 ** (-summed_idr_value), 8) # Convert the sum of -log10IDR values into IDR value (fraction with range from 0 to 1)

    # Insert summed_idr_value and converted_idr_value behind the column "Number of Motifs" in the input fold change calculated peak list file
    output_peak_array_row[(motif_count_column_number + 1):(motif_count_column_number + 1)] = [summed_idr_value, converted_idr_value]
    
    # Append the -log10IDR value columns (from each input IDR output peak list file) at the right end of the the input fold change calculated peak list file
    for idr_value in idr_value_list:
        output_peak_array_row.append(idr_value)

    return output_peak_array_row


# Multiprocessing the IDR value assignment and calculation
print('Integrating IDR values into the peak list')
pool_IDR = multiprocessing.Pool(processes = cpu_count)
output_peak_array = pool_IDR.map(peak_IDR_integration_function, input_peak_array)
pool_IDR.close()
pool_IDR.join()


output_peak_header = input_peak_header

# Insert summed_idr_value and converted_idr_value behind the column "Number of Motifs" in the input fold change calculated peak list file header
output_peak_header[(motif_count_column_number + 1):(motif_count_column_number + 1)] = ['negLog10_IDR', 'IDR']

# Column header names when dealing with narrow peak datasets (SICER2 instead of GEM)
if 'SICER2' in idr_tsv_full_path_list[1]:
    output_peak_header = output_peak_header + ['negLog10_IDR_MACS2_vs_union', 'negLog10_IDR_SICER2_vs_union', 'negLog10_IDR_HOMER_vs_union', 'negLog10_IDR_Genrich_vs_union']

# Column header names when dealing with broad peak datasets (GEM instead of SICER2)
elif 'GEM' in idr_tsv_full_path_list[1]:
    output_peak_header = output_peak_header + ['negLog10_IDR_MACS2_vs_union', 'negLog10_IDR_GEM_vs_union', 'negLog10_IDR_HOMER_vs_union', 'negLog10_IDR_Genrich_vs_union']


# Writing the output file (which is a modified version of the input file)
# ChIP-AP: dataset_name_all_peaks_calculated.tsv (same name, peaks are now with IDR values)
print('Writing the result to file {}'.format(output_tsv_full_path))
output_peak_df = pd.DataFrame(data = output_peak_array, columns = output_peak_header)
output_peak_df.to_csv(output_tsv_full_path, sep = '\t', index = False)
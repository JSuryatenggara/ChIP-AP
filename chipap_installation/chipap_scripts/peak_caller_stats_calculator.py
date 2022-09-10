#!/usr/bin/env python3
#pyright: reportUnboundVariable=false


script_version = '1.0'


print('Running statistics analysis on each peak caller combination')

print('Importing required modules')
import sys
import csv
import os
import numpy as np

def combi_stats(input_tsv, output_tsv):
    # Assigning full paths to the input (HOMER annotated + fold_change_calculator processed peak list) 
    #   and the resulting output (peak caller combination stats table)
    input_tsv_full_path = os.path.abspath(input_tsv)
    output_tsv_full_path = os.path.abspath(output_tsv)
    print('Opening file: {}'.format(input_tsv_full_path))

    # Opening the HOMER annotated + fold_change_calculator processed peak list 
    with open(input_tsv_full_path) as peak_file:
        peak_read = [peak_read_line for peak_read_line in peak_file.readlines()]
        peak_header = peak_read[0].strip().split('\t')
        peak_data = [peak_read[peak_read_counter].strip() for peak_read_counter in range(1,len(peak_read),1)]
        peak_array = []
        
        chip_tag_count_column_number_list = []
        fold_change_column_number_list = []
        


        print('Parsing header names for columns containing target data')
        for peak_header_counter in range(len(peak_header)):

            # Parsing for the column containing information of which peak callers called the peaks in the peak list
            if peak_header[peak_header_counter] == 'Peak Caller Combination':
                peak_caller_ID_column_number = peak_header_counter

            # Parsing for the column containing number of motifs that overlaps the peaks in the peak list. 
            #   Would be all zero if -motif flag was not used during annotatePeaks.pl run
            if peak_header[peak_header_counter] == 'Number of Motifs':
                motif_count_column_number = peak_header_counter

            # Parsing all columns containing number of ChIP reads within the peak. 
            #   Number of columns = number of replicates as different .bam file gives different number of reads
            if 'ChIP Tag Count' in peak_header[peak_header_counter]:
                chip_tag_count_column_number_list.append(peak_header_counter)

            # Parsing all columns containing ChIP vs control reads fold change value at the weighted peak center. 
            #   Number of columns = number of replicates for the same reason
            if 'Fold Change' in peak_header[peak_header_counter]:
                fold_change_column_number_list.append(peak_header_counter)
        


        print('Reading the opened file and loading the contents into arrays')
        # Split all read lines (tab-separated) into a list of values			
        for peak_data_line in peak_data:
            peak_array.append(peak_data_line.split('\t'))

        print('Storing target data in form of list variables')
        # Generate a list of peak caller combinations
        peak_caller_ID_list = [peak_array_row[peak_caller_ID_column_number] for peak_array_row in peak_array]

        # Generate a list of motif counts
        motif_count_list = [int(peak_array_row[motif_count_column_number]) for peak_array_row in peak_array]
        
        # Generate a list that contains the average value of the ChIP read counts in all replicates 
        chip_tag_count_list = [np.mean(np.array([peak_array_row[chip_tag_count_column_number] for chip_tag_count_column_number in chip_tag_count_column_number_list]).astype(np.float)) for peak_array_row in peak_array]

        # Generate a list that contains the average value of the ChIP vs control fold change value in all replicates 
        fold_change_list = [np.mean(np.array([peak_array_row[fold_change_column_number] for fold_change_column_number in fold_change_column_number_list]).astype(np.float)) for peak_array_row in peak_array]



        print('Iterating through each peak. Stats of peaks with the same "Peak Caller Combination" are combined and summarized')
        # Getting all possible combinations of peak callers from the HOMER annotated + fold_change_calculator processed peak list
        peak_caller_dict = {}
        for peak_caller_ID_row in peak_caller_ID_list:
            if peak_caller_ID_row not in peak_caller_dict:
                peak_caller_dict.update({peak_caller_ID_row:0})

        # Prepare all the stats counters for all possible peak caller combinations
        for peak_caller_dict_key, peak_caller_dict_value in peak_caller_dict.items():
            motif_count_exclusive       = 0
            chip_read_sum_exclusive     = 0
            fold_change_sum_exclusive   = 0
            peak_count_exclusive        = 0
            pos_peak_count_exclusive    = 0
            motif_count_inclusive       = 0
            chip_read_sum_inclusive     = 0
            fold_change_sum_inclusive   = 0
            peak_count_inclusive        = 0
            pos_peak_count_inclusive    = 0

            # Original peak caller combination is pipe-separated ('|'). 
            #   This needs to be splitted into a list for inclusive peak caller combination stats
            peak_caller_dict_key_split = peak_caller_dict_key.split('|')



            # Exclusive meaning: peaks called by A and B is NOT included in peaks called by A. 
            #   Peak sets from different combinations are exclusive to each other. AND boolean rule.

            for peak_caller_ID_list_counter in range(len(peak_caller_ID_list)):
                # When the iterated peak caller combination completely matches with 
                #   the one currently accessed in the dictionary, do the things below:
                if peak_caller_ID_list[peak_caller_ID_list_counter] == peak_caller_dict_key:

                    # Add the number of motifs of the iterated peak entry to the motif count counter
                    motif_count_exclusive += motif_count_list[peak_caller_ID_list_counter]

                    # Add the chip reads count of the iterated peak entry to the chip tag count counter
                    chip_read_sum_exclusive += chip_tag_count_list[peak_caller_ID_list_counter]

                    # Add the fold change value of the iterated peak entry to the fold change counter
                    fold_change_sum_exclusive += fold_change_list[peak_caller_ID_list_counter]

                    # Add one to the peak count counter
                    peak_count_exclusive += 1
                    
                    if motif_count_list[peak_caller_ID_list_counter] != 0:

                        # If the iterated peak entry's number of motif is not zero, add one to the positive peak counter
                        pos_peak_count_exclusive += 1

                    # Hit rate and average values below are basically the values above divided by 
                    #   the total number of peak entry as recorded in the peak count counter
                    motif_hit_rate_exclusive        = round(motif_count_exclusive/peak_count_exclusive, 2)
                    pos_peak_hit_rate_exclusive     = round(pos_peak_count_exclusive/peak_count_exclusive, 2)
                    chip_read_average_exclusive     = round(chip_read_sum_exclusive/peak_count_exclusive, 2)
                    fold_change_average_exclusive   = round(fold_change_sum_exclusive/peak_count_exclusive, 2)



                # Inclusive meaning: peaks called by A and B is included in peaks called by A. 
                #   Peak sets from different combinations are shared among other related sets. OR boolean rule.

                # Splits the iterated peak peak caller combination into a list
                peak_caller_ID_split = peak_caller_ID_list[peak_caller_ID_list_counter].split('|')

                # When all elements of the iterated peak caller combination are present in 
                #   the list currently accessed in the dictionary, do the things below:
                if all(element in peak_caller_ID_split for element in peak_caller_dict_key_split):

                    # Add the number of motifs of the iterated peak entry to the motif count counter
                    motif_count_inclusive += motif_count_list[peak_caller_ID_list_counter] 

                    # Add the chip reads count of the iterated peak entry to the chip tag count counter
                    chip_read_sum_inclusive += chip_tag_count_list[peak_caller_ID_list_counter]

                    # Add the fold change value of the iterated peak entry to the fold change counter
                    fold_change_sum_inclusive += fold_change_list[peak_caller_ID_list_counter]
                    
                    # Add one to the peak count counter
                    peak_count_inclusive += 1 

                    if motif_count_list[peak_caller_ID_list_counter] != 0:

                        # If the iterated peak entry's number of motif is not zero, add one to the positive peak counter
                        pos_peak_count_inclusive += 1

                    # Hit rate and average values below are basically the values above divided by 
                    #   the total number of peak entry as recorded in the peak count counter
                    motif_hit_rate_inclusive        = round(motif_count_inclusive/peak_count_inclusive, 2)
                    pos_peak_hit_rate_inclusive     = round(pos_peak_count_inclusive/peak_count_inclusive, 2)
                    chip_read_average_inclusive     = round(chip_read_sum_inclusive/peak_count_inclusive, 2)
                    fold_change_average_inclusive   = round(fold_change_sum_inclusive/peak_count_inclusive, 2)



                # Update the original values (zeroes) of the keys in the dictionary with 
                #   the statistical counters obtained from the iterations above
                peak_caller_dict.update({peak_caller_dict_key:[peak_count_exclusive, 
                                                            pos_peak_count_exclusive, 
                                                            motif_count_exclusive, 
                                                            pos_peak_hit_rate_exclusive, 
                                                            motif_hit_rate_exclusive, 
                                                            chip_read_average_exclusive, 
                                                            fold_change_average_exclusive, 
                                                            peak_count_inclusive, 
                                                            pos_peak_count_inclusive, 
                                                            motif_count_inclusive, 
                                                            pos_peak_hit_rate_inclusive, 
                                                            motif_hit_rate_inclusive, 
                                                            chip_read_average_inclusive, 
                                                            fold_change_average_inclusive]})



    print('Writing a summary table: {}'.format(output_tsv_full_path))
    with open(output_tsv_full_path, 'w') as output_tsv:
    	# Writing the output: Summary table containing peak caller performance-related of the called peaks, 
        #   grouped by the combinations of peak callers that called each of them.
        w = csv.writer(output_tsv, delimiter='\t')
        w.writerow(['peak_callers', 
                    'peak_count_exclusive', 
                    'pos_peak_count_exclusive', 
                    'motif_count_exclusive', 
                    'pos_peak_hit_rate_exclusive', 
                    'motif_hit_rate_exclusive', 
                    'chip_read_average_exclusive', 
                    'fold_change_average_exclusive', 
                    'peak_count_inclusive', 
                    'pos_peak_count_inclusive', 
                    'motif_count_inclusive', 
                    'pos_peak_hit_rate_inclusive', 
                    'motif_hit_rate_inclusive', 
                    'chip_read_average_inclusive', 
                    'fold_change_average_inclusive'])
        
        for peak_caller_dict_key, peak_caller_dict_value in peak_caller_dict.items():
            # Write the dictionary key (peak caller combination) in front all peak statistics values. 
            #   Will be the leftmost column in the resulting table.
            peak_caller_dict_value.insert(0, peak_caller_dict_key)
            w.writerow(peak_caller_dict_value)
    
    print('All processes finished! The output file is: {}'.format(output_tsv_full_path))

combi_stats(sys.argv[1], sys.argv[2])
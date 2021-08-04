#!/usr/bin/env python3
#pyright: reportUnboundVariable=false


script_version = '1.2'


# PATCH NOTES
#   Version 1.1     Written due to some term lists exceeding characters limit per cell when viewed in spreadsheet application
#                       Bracketed string removed from each term to save character space in each term list
#                       If a term list still exceeds the character limit, anything beyond the 32000th character is removed
#   Version 1.2     Added a function to remove 'IDR vs Peakset' column(s) if present



print('Running GO annotation script')

print('Importing required modules')
import csv
import os
import multiprocessing
import argparse
import subprocess

print('Setting argument parser')

parser = argparse.ArgumentParser()

parser.add_argument('--thread', 
                    help = '<Optional> Maximum number of processes to use. Default is half the maximum available.', 
                    type = int, 
                    choices = range(1, (multiprocessing.cpu_count() + 1), 1),
                    metavar = "[1-{}]".format(multiprocessing.cpu_count()),
                    default = int(multiprocessing.cpu_count() / 2))

parser.add_argument('--input_tsv', 
                    help = '<Required> Input peak list file, in HOMER annotatePeaks format, .tsv extension.', 
                    required = True)

parser.add_argument('--output_tsv', 
                    help = '<Required> Output peak list file, in HOMER annotatePeaks format (modified), .tsv extension', 
                    required = True)

parser.add_argument('--go_folder', 
                    help = '<Required> The "gene_ontology" folder, resulting from HOMER annotatePeaks with -go function on', 
                    required = True)

args = parser.parse_args()



subprocess.run('ulimit -n 2000', shell = True)

print('Parsing arguments')
cpu_count = args.thread

print('Establishing paths to the HOMER gene ontology files')
input_tsv_full_path 	= os.path.abspath(args.input_tsv)
output_tsv_full_path 	= os.path.abspath(args.output_tsv)
gene_ontology_dir 		= os.path.abspath(args.go_folder)

# Assigning full paths to the databases. These "filename.txt" are direct output of 
# 	HOMER annotatePeaks.pl with -go switch toggled on
biological_process_full_path = gene_ontology_dir + "/biological_process.txt"
molecular_function_full_path = gene_ontology_dir + "/molecular_function.txt"
cellular_component_full_path = gene_ontology_dir + "/cellular_component.txt"


########################################################################################################################


print('Defining required functions for multiprocessing')

def biological_process_annotation_function(biological_process_annotation_function_arg):
    peak_current_entrez_ID = peak_array[biological_process_annotation_function_arg][peak_entrez_ID_column_number]
    biological_process_term_field_list = []
    for biological_process_array_filtered_counter in range(len(biological_process_array_filtered)):
        # Get all the entrez IDs (comma-separated string) associated to the biological process term in the current iteration
        biological_process_current_entrez_ID = biological_process_array_filtered[biological_process_array_filtered_counter][biological_process_entrez_ID_column_number].split(',')
        for biological_process_current_entrez_ID_elem in biological_process_current_entrez_ID:
            if biological_process_current_entrez_ID_elem == peak_current_entrez_ID:
                biological_process_term = biological_process_array_filtered[biological_process_array_filtered_counter][biological_process_term_column_number].split('(')[0].strip()
                # Append all biological_process annotations associated to the entrez ID of the currently iterated peak into a list
                if biological_process_term not in biological_process_term_field_list:
                    biological_process_term_field_list.append(biological_process_term)
    
    # Join the list into a comma-separated string
    biological_process_term_field_list_joined = ', '.join(biological_process_term_field_list)
    if len(biological_process_term_field_list_joined) > 32000:
        biological_process_term_field_list_joined = biological_process_term_field_list_joined[:32001]
        biological_process_term_field_list_joined = biological_process_term_field_list_joined[:biological_process_term_field_list_joined.rfind(',')]
    
    if biological_process_term_field_list:
        # Return the string if not empty
        return biological_process_term_field_list_joined
    else:
        # Return '-' (dash) if the list is empty. Helps with visualization in spreadsheet: 
        # 	prevents adjacent cells with long value from overflowing towards an empty cell.
        return '-'

def molecular_function_annotation_function(molecular_function_annotation_function_arg):
    peak_current_entrez_ID = peak_array[molecular_function_annotation_function_arg][peak_entrez_ID_column_number]
    molecular_function_term_field_list = []
    for molecular_function_array_filtered_counter in range(len(molecular_function_array_filtered)):
        # Get all the entrez IDs (comma-separated string) associated to the molecular function term in the current iteration
        molecular_function_current_entrez_ID = molecular_function_array_filtered[molecular_function_array_filtered_counter][molecular_function_entrez_ID_column_number].split(',')
        for molecular_function_current_entrez_ID_elem in molecular_function_current_entrez_ID:
            if molecular_function_current_entrez_ID_elem == peak_current_entrez_ID:
                molecular_function_term = molecular_function_array_filtered[molecular_function_array_filtered_counter][molecular_function_term_column_number].split('(')[0].strip()
                # Append all molecular_function annotations associated to the entrez ID of the currently iterated peak into a list
                if molecular_function_term not in molecular_function_term_field_list:
                    molecular_function_term_field_list.append(molecular_function_term)
    
    # Join the list into a comma-separated string
    molecular_function_term_field_list_joined = ', '.join(molecular_function_term_field_list)
    if len(molecular_function_term_field_list_joined) > 32000:
        molecular_function_term_field_list_joined = molecular_function_term_field_list_joined[:32001]
        molecular_function_term_field_list_joined = molecular_function_term_field_list_joined[:molecular_function_term_field_list_joined.rfind(',')]
    
    if molecular_function_term_field_list:
        # Return the string if not empty
        return molecular_function_term_field_list_joined
    else:
        # Return '-' (dash) if the list is empty. Helps with visualization in spreadsheet: 
        # 	prevents adjacent cells with long value from overflowing towards an empty cell.
        return '-'

def cellular_component_annotation_function(cellular_component_annotation_function_arg):
    peak_current_entrez_ID = peak_array[cellular_component_annotation_function_arg][peak_entrez_ID_column_number]
    cellular_component_term_field_list = []
    for cellular_component_array_filtered_counter in range(len(cellular_component_array_filtered)):
        # Get all the entrez IDs (comma-separated string) associated to the cellular component term in the current iteration
        cellular_component_current_entrez_ID = cellular_component_array_filtered[cellular_component_array_filtered_counter][cellular_component_entrez_ID_column_number].split(',')
        for cellular_component_current_entrez_ID_elem in cellular_component_current_entrez_ID:
            if cellular_component_current_entrez_ID_elem == peak_current_entrez_ID:
                cellular_component_term = cellular_component_array_filtered[cellular_component_array_filtered_counter][cellular_component_term_column_number].split('(')[0].strip()
                # Append all cellular_component annotations associated to the entrez ID of the currently iterated peak into a list
                if cellular_component_term not in cellular_component_term_field_list:
                    cellular_component_term_field_list.append(cellular_component_term)
    
    # Join the list into a comma-separated string
    cellular_component_term_field_list_joined = ', '.join(cellular_component_term_field_list)
    if len(cellular_component_term_field_list_joined) > 32000:
        cellular_component_term_field_list_joined = cellular_component_term_field_list_joined[:32001]
        cellular_component_term_field_list_joined = cellular_component_term_field_list_joined[:cellular_component_term_field_list_joined.rfind(',')]
    
    if cellular_component_term_field_list:
        # Return the string if not empty
        return cellular_component_term_field_list_joined
    else:
        # Return '-' (dash) if the list is empty. Helps with visualization in spreadsheet: 
        # 	prevents adjacent cells with long value from overflowing towards an empty cell.
        return '-'


########################################################################################################################


if __name__ == '__main__':
    print('Opening file: {}'.format(input_tsv_full_path))
    # Opening the annotated peak list
    with open(input_tsv_full_path) as peak_file:
        peak_read 			= [peak_read_line for peak_read_line in peak_file.readlines()]
        peak_header 		= peak_read[0].strip().split('\t')
        peak_header_insert 	= []
        peak_data 			= [peak_read[peak_read_counter].strip() for peak_read_counter in range(1,len(peak_read),1)]
        peak_array 			= []

        print('Parsing header names for columns containing target data')
        # Parsing for the column containing Entrez ID in the annotated peak list

        individual_IDR_columns_present = 'no'

        for peak_header_counter in range(len(peak_header)):
            if 'Entrez' in peak_header[peak_header_counter]:
                peak_entrez_ID_column_number = peak_header_counter

            if 'IDR vs Peakset' in peak_header[peak_header_counter]:
                first_peak_IDR_column_number = peak_header_counter
                individual_IDR_columns_present = 'yes'
                peak_header = peak_header[:first_peak_IDR_column_number]
                break
                
        print('Reading the opened file and loading the contents into arrays')
        for peak_data_line in peak_data:
            if individual_IDR_columns_present == 'yes':
                peak_array.append((peak_data_line.split('\t')[:first_peak_IDR_column_number]))
            if individual_IDR_columns_present == 'no':
                peak_array.append(peak_data_line.split('\t'))


########################################################################################################################


        print('Opening file: biological_process.txt')
        # Opening the HOMER biological process database file located in the path assigned by --go_folder flag
        with open(biological_process_full_path) as biological_process_file:
            biological_process_read 			= [biological_process_line for biological_process_line in biological_process_file.readlines()]
            biological_process_header 			= biological_process_read[0].strip().split('\t')
            biological_process_data 			= [biological_process_read[biological_process_read_counter].strip() for biological_process_read_counter in range(1,len(biological_process_read),1)]
            biological_process_array 			= []
            biological_process_array_filtered 	= []

            print('Parsing header names for columns containing Biological Process terms')
            # Parsing for column containing Entrez ID in the database file
            for biological_process_header_counter in range(len(biological_process_header)):
                if 'Entrez' in biological_process_header[biological_process_header_counter]:
                    biological_process_entrez_ID_column_number = biological_process_header_counter
            
            # Parsing for column containing target genes in the database file
            for biological_process_header_counter in range(len(biological_process_header)):
                if biological_process_header[biological_process_header_counter] == 'Target Genes in Term':
                    biological_process_target_genes_count_column_number = biological_process_header_counter

            # Parsing for column containing biological process terms in the database file
            for biological_process_header_counter in range(len(biological_process_header)):
                if biological_process_header[biological_process_header_counter] == 'Term':
                    biological_process_term_column_number = biological_process_header_counter

            # Split all read lines (tab-separated) into a list of values			
            for biological_process_data_line in biological_process_data:
                biological_process_array.append(biological_process_data_line.split('\t'))

            # Filter out all biological process terms that has no registered genes (and thus, no Entrez ID)
            for biological_process_array_row in biological_process_array:
                if biological_process_array_row[biological_process_target_genes_count_column_number] != '0':
                    if biological_process_array_row[biological_process_term_column_number] != 'NA':
                        biological_process_array_filtered.append(biological_process_array_row)

            print('Matching, extracting, and concatenating relevant Biological Process terms for each peak, then storing them in form of list variable')
            # Get the biological process terms from the database relevant to the genes' Entrez ID in each annotated peak
            pool_biological_process = multiprocessing.Pool(processes = cpu_count)
            biological_process_term_array = pool_biological_process.map(biological_process_annotation_function, range(len(peak_array)))
            pool_biological_process.close()
            pool_biological_process.join()

        print('Appending a new column (Biological Process) into the peak arrays')
        peak_header_insert.append('Biological Process')
        for peak_array_counter in range(len(peak_array)):
            # Append the biological process terms as the last column of the peak list
            peak_array[peak_array_counter].append(biological_process_term_array[peak_array_counter])


########################################################################################################################


        print('Opening file: molecular_function.txt')
        # Opening the HOMER molecular function database file located in the path assigned by --go_folder flag
        with open(molecular_function_full_path) as molecular_function_file:
            molecular_function_read 			= [molecular_function_line for molecular_function_line in molecular_function_file.readlines()]
            molecular_function_header 			= molecular_function_read[0].strip().split('\t')
            molecular_function_data 			= [molecular_function_read[molecular_function_read_counter].strip() for molecular_function_read_counter in range(1,len(molecular_function_read),1)]
            molecular_function_array 			= []
            molecular_function_array_filtered 	= []

            print('Parsing header names for columns containing Molecular Function terms')
            # Parsing for column containing Entrez ID in the database file
            for molecular_function_header_counter in range(len(molecular_function_header)):
                if 'Entrez' in molecular_function_header[molecular_function_header_counter]:
                    molecular_function_entrez_ID_column_number = molecular_function_header_counter

            # Parsing for column containing target genes in the database file
            for molecular_function_header_counter in range(len(molecular_function_header)):
                if molecular_function_header[molecular_function_header_counter] == 'Target Genes in Term':
                    molecular_function_target_genes_count_column_number = molecular_function_header_counter

            # Parsing for column containing molecular function terms in the database file
            for molecular_function_header_counter in range(len(molecular_function_header)):
                if molecular_function_header[molecular_function_header_counter] == 'Term':
                    molecular_function_term_column_number = molecular_function_header_counter

            # Split all read lines (tab-separated) into a list of values			
            for molecular_function_data_line in molecular_function_data:
                molecular_function_array.append(molecular_function_data_line.split('\t'))

            # Filter out all molecular function terms that has no registered genes (and thus, no Entrez ID)
            for molecular_function_array_row in molecular_function_array:
                if molecular_function_array_row[molecular_function_target_genes_count_column_number] != '0':
                    if molecular_function_array_row[molecular_function_term_column_number] != 'NA':
                        molecular_function_array_filtered.append(molecular_function_array_row)

            print('Matching, extracting, and concatenating relevant Molecular Function terms for each peak, then storing them in form of list variable')
            # Get the molecular function terms from the database relevant to the genes' Entrez ID in each annotated peak
            pool_molecular_function = multiprocessing.Pool(processes = cpu_count)
            molecular_function_term_array = pool_molecular_function.map(molecular_function_annotation_function, range(len(peak_array)))
            pool_molecular_function.close()
            pool_molecular_function.join()

        print('Appending a new column (Molecular Function) into the peak arrays')
        peak_header_insert.append('Molecular Function')
        for peak_array_counter in range(len(peak_array)):
            # Append the molecular function terms as the last column of the peak list
            peak_array[peak_array_counter].append(molecular_function_term_array[peak_array_counter])


########################################################################################################################


        print('Opening file: cellular_component.txt')
        # Opening the HOMER cellular component database file located in the path assigned by --go_folder flag
        with open(cellular_component_full_path) as cellular_component_file:
            cellular_component_read 			= [cellular_component_line for cellular_component_line in cellular_component_file.readlines()]
            cellular_component_header 			= cellular_component_read[0].strip().split('\t')
            cellular_component_data 			= [cellular_component_read[cellular_component_read_counter].strip() for cellular_component_read_counter in range(1,len(cellular_component_read),1)]
            cellular_component_array 			= []
            cellular_component_array_filtered 	= []

            print('Parsing header names for columns containing Cellular Component terms')
            # Parsing for column containing Entrez ID in the database file
            for cellular_component_header_counter in range(len(cellular_component_header)):
                if 'Entrez' in cellular_component_header[cellular_component_header_counter]:
                    cellular_component_entrez_ID_column_number = cellular_component_header_counter

            # Parsing for column containing target genes in the database file
            for cellular_component_header_counter in range(len(cellular_component_header)):
                if cellular_component_header[cellular_component_header_counter] == 'Target Genes in Term':
                    cellular_component_target_genes_count_column_number = cellular_component_header_counter

            # Parsing for column containing cellular component terms in the database file
            for cellular_component_header_counter in range(len(cellular_component_header)):
                if cellular_component_header[cellular_component_header_counter] == 'Term':
                    cellular_component_term_column_number = cellular_component_header_counter
            
            # Split all read lines (tab-separated) into a list of values			
            for cellular_component_data_line in cellular_component_data:
                cellular_component_array.append(cellular_component_data_line.split('\t'))

            # Filter out all cellular component terms that has no registered genes (and thus, no Entrez ID)
            for cellular_component_array_row in cellular_component_array:
                if cellular_component_array_row[cellular_component_target_genes_count_column_number] != '0':
                    if cellular_component_array_row[cellular_component_term_column_number] != 'NA':
                        cellular_component_array_filtered.append(cellular_component_array_row)

            print('Matching, extracting, and concatenating relevant Cellular Component terms for each peak, then storing them in form of list variable')
            # Get the cellular component terms from the database relevant to the genes' Entrez ID in each annotated peak
            pool_cellular_component = multiprocessing.Pool(processes = cpu_count)
            cellular_component_term_array = pool_cellular_component.map(cellular_component_annotation_function, range(len(peak_array)))
            pool_cellular_component.close()
            pool_cellular_component.join()

        print('Appending a new column (Cellular Component) into the peak arrays')
        peak_header_insert.append('Cellular Component')
        for peak_array_counter in range(len(peak_array)):
            # Append the cellular component terms as the last column of the peak list
            peak_array[peak_array_counter].append(cellular_component_term_array[peak_array_counter])


########################################################################################################################


    print('Writing the updated peak arrays into a new file: {}'.format(output_tsv_full_path))
    # Writing the output: HOMER annotated peak list, appended with terms from gene ontology databases
    with open(output_tsv_full_path, 'w') as output_tsv:
        w = csv.writer(output_tsv, delimiter='\t')
        peak_header_new = peak_header + peak_header_insert
        w.writerow(peak_header_new)
        for peak_array_row in peak_array:
            w.writerow(peak_array_row)

    print('All processes finished! The output file is: {}'.format(output_tsv_full_path))
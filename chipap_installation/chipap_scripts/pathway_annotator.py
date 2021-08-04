#!/usr/bin/env python3
#pyright: reportUnboundVariable=false


script_version = '1.2'


# PATCH NOTES
#   Version 1.1     Written due to some term lists exceeding characters limit per cell when viewed in spreadsheet application
#                       Bracketed string removed from each term to save character space in each term list
#                       If a term list still exceeds the character limit, anything beyond the 32000th character is removed
#   Version 1.2     Added a function to remove 'IDR vs Peakset' column(s) if present



print('Running pathway annotation script')

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

print('Establishing paths to the HOMER pathway files')
input_tsv_full_path 	= os.path.abspath(args.input_tsv)
output_tsv_full_path 	= os.path.abspath(args.output_tsv)
gene_ontology_dir 		= os.path.abspath(args.go_folder)

# Assigning full paths to the databases. These "filename.txt" are direct output of 
# 	HOMER annotatePeaks.pl with -go switch toggled on
common_interactions_full_path 	= gene_ontology_dir + "/interactions.txt"
cosmic_full_path 				= gene_ontology_dir + "/cosmic.txt"
kegg_full_path 					= gene_ontology_dir + "/kegg.txt"
biocyc_full_path 				= gene_ontology_dir + "/biocyc.txt"
pathwayInteractionDB_full_path 	= gene_ontology_dir + "/pathwayInteractionDB.txt"
reactome_full_path 				= gene_ontology_dir + "/reactome.txt"
smpdb_full_path 				= gene_ontology_dir + "/smpdb.txt"
wikipathways_full_path 			= gene_ontology_dir + "/wikipathways.txt"


########################################################################################################################


print('Defining required functions for multiprocessing')

def common_interactions_annotation_function(common_interactions_annotation_function_arg):
    peak_current_entrez_ID = peak_array[common_interactions_annotation_function_arg][peak_entrez_ID_column_number]
    common_interactions_term_field_list = []
    for common_interactions_array_filtered_counter in range(len(common_interactions_array_filtered)):
        # Get all the entrez IDs (comma-separated string) associated to the common interactions term in the current iteration
        common_interactions_current_entrez_ID = common_interactions_array_filtered[common_interactions_array_filtered_counter][common_interactions_entrez_ID_column_number].split(',')
        for common_interactions_current_entrez_ID_elem in common_interactions_current_entrez_ID:
            if common_interactions_current_entrez_ID_elem == peak_current_entrez_ID:
                common_interactions_term = common_interactions_array_filtered[common_interactions_array_filtered_counter][common_interactions_term_column_number].split('(')[0].strip()
                # Append all common_interactions annotations associated to the entrez ID of the currently iterated peak into a list
                if common_interactions_term not in common_interactions_term_field_list:
                    common_interactions_term_field_list.append(common_interactions_term)
    
    # Join the list into a comma-separated string
    common_interactions_term_field_list_joined = ', '.join(common_interactions_term_field_list)
    if len(common_interactions_term_field_list_joined) > 32000:
        common_interactions_term_field_list_joined = common_interactions_term_field_list_joined[:32001]
        common_interactions_term_field_list_joined = common_interactions_term_field_list_joined[:common_interactions_term_field_list_joined.rfind(',')]
    
    if common_interactions_term_field_list:
        # Return the string if not empty
        return common_interactions_term_field_list_joined
    else:
        # Return '-' (dash) if the list is empty. Helps with visualization in spreadsheet: 
        # 	prevents adjacent cells with long value from overflowing towards an empty cell.
        return '-'

def cosmic_annotation_function(cosmic_annotation_function_arg):
    peak_current_entrez_ID = peak_array[cosmic_annotation_function_arg][peak_entrez_ID_column_number]
    cosmic_term_field_list = []
    for cosmic_array_filtered_counter in range(len(cosmic_array_filtered)):
        # Get all the entrez IDs (comma-separated string) associated to the cosmic term in the current iteration
        cosmic_current_entrez_ID = cosmic_array_filtered[cosmic_array_filtered_counter][cosmic_entrez_ID_column_number].split(',')
        for cosmic_current_entrez_ID_elem in cosmic_current_entrez_ID:
            if cosmic_current_entrez_ID_elem == peak_current_entrez_ID:
                cosmic_term = cosmic_array_filtered[cosmic_array_filtered_counter][cosmic_term_column_number].split('(')[0].strip()
                # Append all cosmic annotations associated to the entrez ID of the currently iterated peak into a list
                if cosmic_term not in cosmic_term_field_list:
                    cosmic_term_field_list.append(cosmic_term)
    
    # Join the list into a comma-separated string
    cosmic_term_field_list_joined = ', '.join(cosmic_term_field_list)
    if len(cosmic_term_field_list_joined) > 32000:
        cosmic_term_field_list_joined = cosmic_term_field_list_joined[:32001]
        cosmic_term_field_list_joined = cosmic_term_field_list_joined[:cosmic_term_field_list_joined.rfind(',')]
    
    if cosmic_term_field_list:
        # Return the string if not empty
        return cosmic_term_field_list_joined
    else:
        # Return '-' (dash) if the list is empty. Helps with visualization in spreadsheet: 
        # 	prevents adjacent cells with long value from overflowing towards an empty cell.
        return '-'

def kegg_annotation_function(kegg_annotation_function_arg):
    peak_current_entrez_ID = peak_array[kegg_annotation_function_arg][peak_entrez_ID_column_number]
    kegg_term_field_list = []
    for kegg_array_filtered_counter in range(len(kegg_array_filtered)):
        # Get all the entrez IDs (comma-separated string) associated to the kegg term in the current iteration
        kegg_current_entrez_ID = kegg_array_filtered[kegg_array_filtered_counter][kegg_entrez_ID_column_number].split(',')
        for kegg_current_entrez_ID_elem in kegg_current_entrez_ID:
            if kegg_current_entrez_ID_elem == peak_current_entrez_ID:
                kegg_term = kegg_array_filtered[kegg_array_filtered_counter][kegg_term_column_number].split('(')[0].strip()
                # Append all kegg annotations associated to the entrez ID of the currently iterated peak into a list
                if kegg_term not in kegg_term_field_list:
                    kegg_term_field_list.append(kegg_term)
    
    # Join the list into a comma-separated string
    kegg_term_field_list_joined = ', '.join(kegg_term_field_list)
    if len(kegg_term_field_list_joined) > 32000:
        kegg_term_field_list_joined = kegg_term_field_list_joined[:32001]
        kegg_term_field_list_joined = kegg_term_field_list_joined[:kegg_term_field_list_joined.rfind(',')]
    
    if kegg_term_field_list:
        # Return the string if not empty
        return kegg_term_field_list_joined
    else:
        # Return '-' (dash) if the list is empty. Helps with visualization in spreadsheet: 
        # 	prevents adjacent cells with long value from overflowing towards an empty cell.
        return '-'

def biocyc_annotation_function(biocyc_annotation_function_arg):
    peak_current_entrez_ID = peak_array[biocyc_annotation_function_arg][peak_entrez_ID_column_number]
    biocyc_term_field_list = []
    for biocyc_array_filtered_counter in range(len(biocyc_array_filtered)):
        # Get all the entrez IDs (comma-separated string) associated to the biocyc term in the current iteration
        biocyc_current_entrez_ID = biocyc_array_filtered[biocyc_array_filtered_counter][biocyc_entrez_ID_column_number].split(',')
        for biocyc_current_entrez_ID_elem in biocyc_current_entrez_ID:
            if biocyc_current_entrez_ID_elem == peak_current_entrez_ID:
                biocyc_term = biocyc_array_filtered[biocyc_array_filtered_counter][biocyc_term_column_number].split('(')[0].strip()
                # Append all biocyc annotations associated to the entrez ID of the currently iterated peak into a list
                if biocyc_term not in biocyc_term_field_list:
                    biocyc_term_field_list.append(biocyc_term)
    
    # Join the list into a comma-separated string
    biocyc_term_field_list_joined = ', '.join(biocyc_term_field_list)
    if len(biocyc_term_field_list_joined) > 32000:
        biocyc_term_field_list_joined = biocyc_term_field_list_joined[:32001]
        biocyc_term_field_list_joined = biocyc_term_field_list_joined[:biocyc_term_field_list_joined.rfind(',')]
    
    if biocyc_term_field_list:
        # Return the string if not empty
        return biocyc_term_field_list_joined
    else:
        # Return '-' (dash) if the list is empty. Helps with visualization in spreadsheet: 
        # 	prevents adjacent cells with long value from overflowing towards an empty cell.
        return '-'

def pathwayInteractionDB_annotation_function(pathwayInteractionDB_annotation_function_arg):
    peak_current_entrez_ID = peak_array[pathwayInteractionDB_annotation_function_arg][peak_entrez_ID_column_number]
    pathwayInteractionDB_term_field_list = []
    for pathwayInteractionDB_array_filtered_counter in range(len(pathwayInteractionDB_array_filtered)):
        # Get all the entrez IDs (comma-separated string) associated to the pathwayInteractionDB term in the current iteration
        pathwayInteractionDB_current_entrez_ID = pathwayInteractionDB_array_filtered[pathwayInteractionDB_array_filtered_counter][pathwayInteractionDB_entrez_ID_column_number].split(',')
        for pathwayInteractionDB_current_entrez_ID_elem in pathwayInteractionDB_current_entrez_ID:
            if pathwayInteractionDB_current_entrez_ID_elem == peak_current_entrez_ID:
                pathwayInteractionDB_term = pathwayInteractionDB_array_filtered[pathwayInteractionDB_array_filtered_counter][pathwayInteractionDB_term_column_number].split('(')[0].strip()
                # Append all pathwayInteractionDB annotations associated to the entrez ID of the currently iterated peak into a list
                if pathwayInteractionDB_term not in pathwayInteractionDB_term_field_list:
                    pathwayInteractionDB_term_field_list.append(pathwayInteractionDB_term)
    
    # Join the list into a comma-separated string
    pathwayInteractionDB_term_field_list_joined = ', '.join(pathwayInteractionDB_term_field_list)
    if len(pathwayInteractionDB_term_field_list_joined) > 32000:
        pathwayInteractionDB_term_field_list_joined = pathwayInteractionDB_term_field_list_joined[:32001]
        pathwayInteractionDB_term_field_list_joined = pathwayInteractionDB_term_field_list_joined[:pathwayInteractionDB_term_field_list_joined.rfind(',')]
    
    if pathwayInteractionDB_term_field_list:
        # Return the string if not empty
        return pathwayInteractionDB_term_field_list_joined
    else:
        # Return '-' (dash) if the list is empty. Helps with visualization in spreadsheet: 
        # 	prevents adjacent cells with long value from overflowing towards an empty cell.
        return '-'

def reactome_annotation_function(reactome_annotation_function_arg):
    peak_current_entrez_ID = peak_array[reactome_annotation_function_arg][peak_entrez_ID_column_number]
    reactome_term_field_list = []
    for reactome_array_filtered_counter in range(len(reactome_array_filtered)):
        # Get all the entrez IDs (comma-separated string) associated to the reactome term in the current iteration
        reactome_current_entrez_ID = reactome_array_filtered[reactome_array_filtered_counter][reactome_entrez_ID_column_number].split(',')
        for reactome_current_entrez_ID_elem in reactome_current_entrez_ID:
            if reactome_current_entrez_ID_elem == peak_current_entrez_ID:
                reactome_term = reactome_array_filtered[reactome_array_filtered_counter][reactome_term_column_number].split('(')[0].strip()
                # Append all reactome annotations associated to the entrez ID of the currently iterated peak into a list
                if reactome_term not in reactome_term_field_list:
                    reactome_term_field_list.append(reactome_term)
    
    # Join the list into a comma-separated string
    reactome_term_field_list_joined = ', '.join(reactome_term_field_list)
    if len(reactome_term_field_list_joined) > 32000:
        reactome_term_field_list_joined = reactome_term_field_list_joined[:32001]
        reactome_term_field_list_joined = reactome_term_field_list_joined[:reactome_term_field_list_joined.rfind(',')]
    
    if reactome_term_field_list:
        # Return the string if not empty
        return reactome_term_field_list_joined
    else:
        # Return '-' (dash) if the list is empty. Helps with visualization in spreadsheet: 
        # 	prevents adjacent cells with long value from overflowing towards an empty cell.
        return '-'

def smpdb_annotation_function(smpdb_annotation_function_arg):
    peak_current_entrez_ID = peak_array[smpdb_annotation_function_arg][peak_entrez_ID_column_number]
    smpdb_term_field_list = []
    for smpdb_array_filtered_counter in range(len(smpdb_array_filtered)):
        # Get all the entrez IDs (comma-separated string) associated to the smpdb term in the current iteration
        smpdb_current_entrez_ID = smpdb_array_filtered[smpdb_array_filtered_counter][smpdb_entrez_ID_column_number].split(',')
        for smpdb_current_entrez_ID_elem in smpdb_current_entrez_ID:
            if smpdb_current_entrez_ID_elem == peak_current_entrez_ID:
                smpdb_term = smpdb_array_filtered[smpdb_array_filtered_counter][smpdb_term_column_number].split('(')[0].strip()
                # Append all smpdb annotations associated to the entrez ID of the currently iterated peak into a list
                if smpdb_term not in smpdb_term_field_list:
                    smpdb_term_field_list.append(smpdb_term)
    
    # Join the list into a comma-separated string
    smpdb_term_field_list_joined = ', '.join(smpdb_term_field_list)
    if len(smpdb_term_field_list_joined) > 32000:
        smpdb_term_field_list_joined = smpdb_term_field_list_joined[:32001]
        smpdb_term_field_list_joined = smpdb_term_field_list_joined[:smpdb_term_field_list_joined.rfind(',')]
    
    if smpdb_term_field_list:
        # Return the string if not empty
        return smpdb_term_field_list_joined
    else:
        # Return '-' (dash) if the list is empty. Helps with visualization in spreadsheet: 
        # 	prevents adjacent cells with long value from overflowing towards an empty cell.
        return '-'

def wikipathways_annotation_function(wikipathways_annotation_function_arg):
    peak_current_entrez_ID = peak_array[wikipathways_annotation_function_arg][peak_entrez_ID_column_number]
    wikipathways_term_field_list = []
    for wikipathways_array_filtered_counter in range(len(wikipathways_array_filtered)):
        # Get all the entrez IDs (comma-separated string) associated to the wikipathways term in the current iteration
        wikipathways_current_entrez_ID = wikipathways_array_filtered[wikipathways_array_filtered_counter][wikipathways_entrez_ID_column_number].split(',')
        for wikipathways_current_entrez_ID_elem in wikipathways_current_entrez_ID:
            if wikipathways_current_entrez_ID_elem == peak_current_entrez_ID:
                wikipathways_term = wikipathways_array_filtered[wikipathways_array_filtered_counter][wikipathways_term_column_number].split('(')[0].strip()
                # Append all wikipathways annotations associated to the entrez ID of the currently iterated peak into a list
                if wikipathways_term not in wikipathways_term_field_list:
                    wikipathways_term_field_list.append(wikipathways_term)
    
    # Join the list into a comma-separated string
    wikipathways_term_field_list_joined = ', '.join(wikipathways_term_field_list)
    if len(wikipathways_term_field_list_joined) > 32000:
        wikipathways_term_field_list_joined = wikipathways_term_field_list_joined[:32001]
        wikipathways_term_field_list_joined = wikipathways_term_field_list_joined[:wikipathways_term_field_list_joined.rfind(',')]
    
    if wikipathways_term_field_list:
        # Return the string if not empty
        return wikipathways_term_field_list_joined
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


        print('Opening file: interactions.txt')
        # Opening the HOMER common interactions database file located in the path assigned by --go_folder flag
        with open(common_interactions_full_path) as common_interactions_file:
            common_interactions_read 			= [common_interactions_line for common_interactions_line in common_interactions_file.readlines()]
            common_interactions_header 			= common_interactions_read[0].strip().split('\t')
            common_interactions_data 			= [common_interactions_read[common_interactions_read_counter].strip() for common_interactions_read_counter in range(1,len(common_interactions_read),1)]
            common_interactions_array 			= []
            common_interactions_array_filtered 	= []

            print('Parsing header names for columns containing Interaction with Common Protein terms')
            # Parsing for column containing Entrez ID in the database file
            for common_interactions_header_counter in range(len(common_interactions_header)):
                if 'Entrez' in common_interactions_header[common_interactions_header_counter]:
                    common_interactions_entrez_ID_column_number = common_interactions_header_counter
            
            # Parsing for column containing target genes in the database file
            for common_interactions_header_counter in range(len(common_interactions_header)):
                if common_interactions_header[common_interactions_header_counter] == 'Target Genes in Term':
                    common_interactions_target_genes_count_column_number = common_interactions_header_counter

            # Parsing for column containing common interaction terms in the database file
            for common_interactions_header_counter in range(len(common_interactions_header)):
                if common_interactions_header[common_interactions_header_counter] == 'Term':
                    common_interactions_term_column_number = common_interactions_header_counter

            # Split all read lines (tab-separated) into a list of values			
            for common_interactions_data_line in common_interactions_data:
                common_interactions_array.append(common_interactions_data_line.split('\t'))

            # Filter out all common interactions terms that has no registered genes (and thus, no Entrez ID)
            for common_interactions_array_row in common_interactions_array:
                if common_interactions_array_row[common_interactions_target_genes_count_column_number] != '0':
                    common_interactions_array_filtered.append(common_interactions_array_row)

            print('Matching, extracting, and concatenating relevant Interaction with Common Protein terms for each peak, then storing them in form of list variable')
            # Get the common interactions terms from the database relevant to the genes' Entrez ID in each annotated peak
            pool_common_interactions = multiprocessing.Pool(processes = cpu_count)
            common_interactions_term_array = pool_common_interactions.map(common_interactions_annotation_function, range(len(peak_array)))
            pool_common_interactions.close()
            pool_common_interactions.join()

        print('Appending a new column (Interaction with Common Protein) into the peak arrays')
        peak_header_insert.append('Interaction with Common Protein')
        for peak_array_counter in range(len(peak_array)):
            # Append the common interactions terms as the last column of the peak list
            peak_array[peak_array_counter].append(common_interactions_term_array[peak_array_counter])


########################################################################################################################


        print('Opening file: cosmic.txt')
        # Opening the HOMER cosmic database file located in the path assigned by --go_folder flag
        with open(cosmic_full_path) as cosmic_file:
            cosmic_read 			= [cosmic_line for cosmic_line in cosmic_file.readlines()]
            cosmic_header 			= cosmic_read[0].strip().split('\t')
            cosmic_data 			= [cosmic_read[cosmic_read_counter].strip() for cosmic_read_counter in range(1,len(cosmic_read),1)]
            cosmic_array 			= []
            cosmic_array_filtered 	= []

            print('Parsing header names for columns containing Somatic Mutations (COSMIC) terms')
            # Parsing for column containing Entrez ID in the database file
            for cosmic_header_counter in range(len(cosmic_header)):
                if 'Entrez' in cosmic_header[cosmic_header_counter]:
                    cosmic_entrez_ID_column_number = cosmic_header_counter
                    
            # Parsing for column containing target genes in the database file
            for cosmic_header_counter in range(len(cosmic_header)):
                if cosmic_header[cosmic_header_counter] == 'Target Genes in Term':
                    cosmic_target_genes_count_column_number = cosmic_header_counter

            # Parsing for column containing cosmic terms in the database file
            for cosmic_header_counter in range(len(cosmic_header)):
                if cosmic_header[cosmic_header_counter] == 'Term':
                    cosmic_term_column_number = cosmic_header_counter

            # Split all read lines (tab-separated) into a list of values			
            for cosmic_data_line in cosmic_data:
                cosmic_array.append(cosmic_data_line.split('\t'))

            # Filter out all cosmic terms that has no registered genes (and thus, no Entrez ID)
            for cosmic_array_row in cosmic_array:
                if cosmic_array_row[cosmic_target_genes_count_column_number] != '0':
                    cosmic_array_filtered.append(cosmic_array_row)

            print('Matching, extracting, and concatenating relevant Somatic Mutations (COSMIC) terms for each peak, then storing them in form of list variable')
            # Get the cosmic terms from the database relevant to the genes' Entrez ID in each annotated peak
            pool_cosmic = multiprocessing.Pool(processes = cpu_count)
            cosmic_term_array = pool_cosmic.map(cosmic_annotation_function, range(len(peak_array)))
            pool_cosmic.close()
            pool_cosmic.join()
        
        print('Appending a new column (Somatic Mutations (COSMIC)) into the peak arrays')
        peak_header_insert.append('Somatic Mutations (COSMIC)')
        for peak_array_counter in range(len(peak_array)):
            # Append the cosmic terms as the last column of the peak list
            peak_array[peak_array_counter].append(cosmic_term_array[peak_array_counter])


########################################################################################################################


        print('Opening file: kegg.txt')
        # Opening the HOMER kegg database file located in the path assigned by --go_folder flag
        with open(kegg_full_path) as kegg_file:
            kegg_read 			= [kegg_line for kegg_line in kegg_file.readlines()]
            kegg_header 		= kegg_read[0].strip().split('\t')
            kegg_data 			= [kegg_read[kegg_read_counter].strip() for kegg_read_counter in range(1,len(kegg_read),1)]
            kegg_array 			= []
            kegg_array_filtered = []

            print('Parsing header names for columns containing Pathway (KEGG) terms')
            # Parsing for column containing Entrez ID in the database file
            for kegg_header_counter in range(len(kegg_header)):
                if 'Entrez' in kegg_header[kegg_header_counter]:
                    kegg_entrez_ID_column_number = kegg_header_counter

            # Parsing for column containing target genes in the database file
            for kegg_header_counter in range(len(kegg_header)):
                if kegg_header[kegg_header_counter] == 'Target Genes in Term':
                    kegg_target_genes_count_column_number = kegg_header_counter

            # Parsing for column containing kegg terms in the database file
            for kegg_header_counter in range(len(kegg_header)):
                if kegg_header[kegg_header_counter] == 'Term':
                    kegg_term_column_number = kegg_header_counter

            # Split all read lines (tab-separated) into a list of values			
            for kegg_data_line in kegg_data:
                kegg_array.append(kegg_data_line.split('\t'))

            # Filter out all kegg terms that has no registered genes (and thus, no Entrez ID)
            for kegg_array_row in kegg_array:
                if kegg_array_row[kegg_target_genes_count_column_number] != '0':
                    kegg_array_filtered.append(kegg_array_row)

            print('Matching, extracting, and concatenating relevant Pathway (KEGG) terms for each peak, then storing them in form of list variable')
            # Get the kegg terms from the database relevant to the genes' Entrez ID in each annotated peak
            pool_kegg = multiprocessing.Pool(processes = cpu_count)
            kegg_term_array = pool_kegg.map(kegg_annotation_function, range(len(peak_array)))
            pool_kegg.close()
            pool_kegg.join()

        print('Appending a new column (Pathway (KEGG)) into the peak arrays')
        peak_header_insert.append('Pathway (KEGG)')
        for peak_array_counter in range(len(peak_array)):
            # Append the kegg terms as the last column of the peak list
            peak_array[peak_array_counter].append(kegg_term_array[peak_array_counter])


########################################################################################################################


        print('Opening file: biocyc.txt')
        # Opening the HOMER biocyc database file located in the path assigned by --go_folder flag
        with open(biocyc_full_path) as biocyc_file:
            biocyc_read 			= [biocyc_line for biocyc_line in biocyc_file.readlines()]
            biocyc_header 			= biocyc_read[0].strip().split('\t')
            biocyc_data 			= [biocyc_read[biocyc_read_counter].strip() for biocyc_read_counter in range(1,len(biocyc_read),1)]
            biocyc_array 			= []
            biocyc_array_filtered 	= []

            print('Parsing header names for columns containing Pathway (BIOCYC) terms')
            # Parsing for column containing Entrez ID in the database file
            for biocyc_header_counter in range(len(biocyc_header)):
                if 'Entrez' in biocyc_header[biocyc_header_counter]:
                    biocyc_entrez_ID_column_number = biocyc_header_counter
            
            # Parsing for column containing target genes in the database file
            for biocyc_header_counter in range(len(biocyc_header)):
                if biocyc_header[biocyc_header_counter] == 'Target Genes in Term':
                    biocyc_target_genes_count_column_number = biocyc_header_counter

            # Parsing for column containing biocyc terms in the database file
            for biocyc_header_counter in range(len(biocyc_header)):
                if biocyc_header[biocyc_header_counter] == 'Term':
                    biocyc_term_column_number = biocyc_header_counter
            
            # Split all read lines (tab-separated) into a list of values			
            for biocyc_data_line in biocyc_data:
                biocyc_array.append(biocyc_data_line.split('\t'))

            # Filter out all biocyc terms that has no registered genes (and thus, no Entrez ID)
            for biocyc_array_row in biocyc_array:
                if biocyc_array_row[biocyc_target_genes_count_column_number] != '0':
                    biocyc_array_filtered.append(biocyc_array_row)

            print('Matching, extracting, and concatenating relevant Pathway (BIOCYC) terms for each peak, then storing them in form of list variable')
            # Get the biocyc terms from the database relevant to the genes' Entrez ID in each annotated peak
            pool_biocyc = multiprocessing.Pool(processes = cpu_count)
            biocyc_term_array = pool_biocyc.map(biocyc_annotation_function, range(len(peak_array)))
            pool_biocyc.close()
            pool_biocyc.join()

        print('Appending a new column (Pathway (BIOCYC)) into the peak arrays')
        peak_header_insert.append('Pathway (BIOCYC)')
        for peak_array_counter in range(len(peak_array)):
            # Append the biocyc terms as the last column of the peak list
            peak_array[peak_array_counter].append(biocyc_term_array[peak_array_counter])


########################################################################################################################


        print('Opening file: pathwayInteractionDB.txt')
        # Opening the HOMER pathwayInteractionDB database file located in the path assigned by --go_folder flag
        with open(pathwayInteractionDB_full_path) as pathwayInteractionDB_file:
            pathwayInteractionDB_read 			= [pathwayInteractionDB_line for pathwayInteractionDB_line in pathwayInteractionDB_file.readlines()]
            pathwayInteractionDB_header 		= pathwayInteractionDB_read[0].strip().split('\t')
            pathwayInteractionDB_data 			= [pathwayInteractionDB_read[pathwayInteractionDB_read_counter].strip() for pathwayInteractionDB_read_counter in range(1,len(pathwayInteractionDB_read),1)]
            pathwayInteractionDB_array 			= []
            pathwayInteractionDB_array_filtered = []

            print('Parsing header names for columns containing Pathway (pathwayInteractionDB) terms')
            # Parsing for column containing Entrez ID in the database file
            for pathwayInteractionDB_header_counter in range(len(pathwayInteractionDB_header)):
                if 'Entrez' in pathwayInteractionDB_header[pathwayInteractionDB_header_counter]:
                    pathwayInteractionDB_entrez_ID_column_number = pathwayInteractionDB_header_counter
            
            # Parsing for column containing target genes in the database file
            for pathwayInteractionDB_header_counter in range(len(pathwayInteractionDB_header)):
                if pathwayInteractionDB_header[pathwayInteractionDB_header_counter] == 'Target Genes in Term':
                    pathwayInteractionDB_target_genes_count_column_number = pathwayInteractionDB_header_counter

            # Parsing for column containing pathwayInteractionDB terms in the database file
            for pathwayInteractionDB_header_counter in range(len(pathwayInteractionDB_header)):
                if pathwayInteractionDB_header[pathwayInteractionDB_header_counter] == 'Term':
                    pathwayInteractionDB_term_column_number = pathwayInteractionDB_header_counter

            # Split all read lines (tab-separated) into a list of values			
            for pathwayInteractionDB_data_line in pathwayInteractionDB_data:
                pathwayInteractionDB_array.append(pathwayInteractionDB_data_line.split('\t'))

            # Filter out all pathwayInteractionDB terms that has no registered genes (and thus, no Entrez ID)
            for pathwayInteractionDB_array_row in pathwayInteractionDB_array:
                if pathwayInteractionDB_array_row[pathwayInteractionDB_target_genes_count_column_number] != '0':
                    pathwayInteractionDB_array_filtered.append(pathwayInteractionDB_array_row)

            print('Matching, extracting, and concatenating relevant Pathway (pathwayInteractionDB) terms for each peak, then storing them in form of list variable')
            # Get the pathwayInteractionDB terms from the database relevant to the genes' Entrez ID in each annotated peak
            pool_pathwayInteractionDB = multiprocessing.Pool(processes = cpu_count)
            pathwayInteractionDB_term_array = pool_pathwayInteractionDB.map(pathwayInteractionDB_annotation_function, range(len(peak_array)))
            pool_pathwayInteractionDB.close()
            pool_pathwayInteractionDB.join()

        print('Appending a new column (Pathway (pathwayInteractionDB)) into the peak arrays')
        peak_header_insert.append('Pathway (pathwayInteractionDB)')
        for peak_array_counter in range(len(peak_array)):
            # Append the pathwayInteractionDB terms as the last column of the peak list
            peak_array[peak_array_counter].append(pathwayInteractionDB_term_array[peak_array_counter])


########################################################################################################################


        print('Opening file: reactome.txt')
        # Opening the HOMER reactome database file located in the path assigned by --go_folder flag
        with open(reactome_full_path, 'rb') as reactome_file:
            reactome_read 			= [reactome_line for reactome_line in reactome_file.readlines()]
            reactome_header 		= reactome_read[0].decode('utf8').strip().split('\t')
            reactome_data 			= [reactome_read[reactome_read_counter].strip() for reactome_read_counter in range(1,len(reactome_read),1)]
            reactome_array 			= []
            reactome_array_filtered = []

            print('Parsing header names for columns containing Pathway (REACTOME) terms')
            # Parsing for column containing Entrez ID in the database file
            for reactome_header_counter in range(len(reactome_header)):
                if 'Entrez' in reactome_header[reactome_header_counter]:
                    reactome_entrez_ID_column_number = reactome_header_counter

            # Parsing for column containing target genes in the database file
            for reactome_header_counter in range(len(reactome_header)):
                if reactome_header[reactome_header_counter] == 'Target Genes in Term':
                    reactome_target_genes_count_column_number = reactome_header_counter

            # Parsing for column containing reactome terms in the database file
            for reactome_header_counter in range(len(reactome_header)):
                if reactome_header[reactome_header_counter] == 'Term':
                    reactome_term_column_number = reactome_header_counter

            # Split all read lines (tab-separated) into a list of values			
            for reactome_data_line in reactome_data:
                reactome_array.append(reactome_data_line.decode(encoding='utf-8', errors='ignore').split('\t'))

            # Filter out all reactome terms that has no registered genes (and thus, no Entrez ID)
            for reactome_array_row in reactome_array:
                if reactome_array_row[reactome_target_genes_count_column_number] != '0':
                    reactome_array_filtered.append(reactome_array_row)

            print('Matching, extracting, and concatenating relevant Pathway (REACTOME) terms for each peak, then storing them in form of list variable')
            # Get the reactome terms from the database relevant to the genes' Entrez ID in each annotated peak
            pool_reactome = multiprocessing.Pool(processes = cpu_count)
            reactome_term_array = pool_reactome.map(reactome_annotation_function, range(len(peak_array)))
            pool_reactome.close()
            pool_reactome.join()

        print('Appending a new column (Pathway (REACTOME)) into the peak arrays')
        peak_header_insert.append('Pathway (REACTOME)')
        for peak_array_counter in range(len(peak_array)):
            # Append the reactome terms as the last column of the peak list
            peak_array[peak_array_counter].append(reactome_term_array[peak_array_counter])


########################################################################################################################


        print('Opening file: smpdb.txt')
        # Opening the HOMER smpdb database file located in the path assigned by --go_folder flag
        with open(smpdb_full_path) as smpdb_file:
            smpdb_read 				= [smpdb_line for smpdb_line in smpdb_file.readlines()]
            smpdb_header 			= smpdb_read[0].strip().split('\t')
            smpdb_data 				= [smpdb_read[smpdb_read_counter].strip() for smpdb_read_counter in range(1,len(smpdb_read),1)]
            smpdb_array 			= []
            smpdb_array_filtered 	= []

            print('Parsing header names for columns containing Pathway (SMPDB) terms')
            # Parsing for column containing Entrez ID in the database file
            for smpdb_header_counter in range(len(smpdb_header)):
                if 'Entrez' in smpdb_header[smpdb_header_counter]:
                    smpdb_entrez_ID_column_number = smpdb_header_counter

            # Parsing for column containing target genes in the database file
            for smpdb_header_counter in range(len(smpdb_header)):
                if smpdb_header[smpdb_header_counter] == 'Target Genes in Term':
                    smpdb_target_genes_count_column_number = smpdb_header_counter

            # Parsing for column containing smpdb terms in the database file
            for smpdb_header_counter in range(len(smpdb_header)):
                if smpdb_header[smpdb_header_counter] == 'Term':
                    smpdb_term_column_number = smpdb_header_counter

            # Split all read lines (tab-separated) into a list of values			
            for smpdb_data_line in smpdb_data:
                smpdb_array.append(smpdb_data_line.split('\t'))

            # Filter out all smpdb terms that has no registered genes (and thus, no Entrez ID)
            for smpdb_array_row in smpdb_array:
                if smpdb_array_row[smpdb_target_genes_count_column_number] != '0':
                    smpdb_array_filtered.append(smpdb_array_row)

            print('Matching, extracting, and concatenating relevant Pathway (SMPDB) terms for each peak, then storing them in form of list variable')
            # Get the smpdb terms from the database relevant to the genes' Entrez ID in each annotated peak
            pool_smpdb = multiprocessing.Pool(processes = cpu_count)
            smpdb_term_array = pool_smpdb.map(smpdb_annotation_function, range(len(peak_array)))
            pool_smpdb.close()
            pool_smpdb.join()

        print('Appending a new column (Pathway (SMPDB)) into the peak arrays')
        peak_header_insert.append('Pathway (SMPDB)')
        for peak_array_counter in range(len(peak_array)):
            # Append the smpdb terms as the last column of the peak list
            peak_array[peak_array_counter].append(smpdb_term_array[peak_array_counter])


########################################################################################################################


        print('Opening file: wikipathways.txt')
        # Opening the HOMER wikipathways database file located in the path assigned by --go_folder flag
        with open(wikipathways_full_path) as wikipathways_file:
            wikipathways_read 			= [wikipathways_line for wikipathways_line in wikipathways_file.readlines()]
            wikipathways_header 		= wikipathways_read[0].strip().split('\t')
            wikipathways_data 			= [wikipathways_read[wikipathways_read_counter].strip() for wikipathways_read_counter in range(1,len(wikipathways_read),1)]
            wikipathways_array 			= []
            wikipathways_array_filtered = []

            print('Parsing header names for columns containing Pathway (Wikipathways) terms')
            # Parsing for column containing Entrez ID in the database file
            for wikipathways_header_counter in range(len(wikipathways_header)):
                if 'Entrez' in wikipathways_header[wikipathways_header_counter]:
                    wikipathways_entrez_ID_column_number = wikipathways_header_counter
            
            # Parsing for column containing target genes in the database file
            for wikipathways_header_counter in range(len(wikipathways_header)):
                if wikipathways_header[wikipathways_header_counter] == 'Target Genes in Term':
                    wikipathways_target_genes_count_column_number = wikipathways_header_counter

            # Parsing for column containing wikipathways terms in the database file
            for wikipathways_header_counter in range(len(wikipathways_header)):
                if wikipathways_header[wikipathways_header_counter] == 'Term':
                    wikipathways_term_column_number = wikipathways_header_counter

            # Split all read lines (tab-separated) into a list of values			
            for wikipathways_data_line in wikipathways_data:
                wikipathways_array.append(wikipathways_data_line.split('\t'))

            # Filter out all wikipathways terms that has no registered genes (and thus, no Entrez ID)
            for wikipathways_array_row in wikipathways_array:
                if wikipathways_array_row[wikipathways_target_genes_count_column_number] != '0':
                    wikipathways_array_filtered.append(wikipathways_array_row)

            print('Matching, extracting, and concatenating relevant Pathway (Wikipathways) terms for each peak, then storing them in form of list variable')
            # Get the wikipathways terms from the database relevant to the genes' Entrez ID in each annotated peak
            pool_wikipathways = multiprocessing.Pool(processes = cpu_count)
            wikipathways_term_array = pool_wikipathways.map(wikipathways_annotation_function, range(len(peak_array)))
            pool_wikipathways.close()
            pool_wikipathways.join()

        print('Appending a new column (Pathway (Wikipathways)) into the peak arrays')
        peak_header_insert.append('Pathway (Wikipathways)')
        for peak_array_counter in range(len(peak_array)):
            # Append the wikipathways terms as the last column of the peak list
            peak_array[peak_array_counter].append(wikipathways_term_array[peak_array_counter])


########################################################################################################################


    print('Writing the updated peak arrays into a new file: {}'.format(output_tsv_full_path))
    # Writing the output: HOMER annotated peak list, appended with terms from 
    # 	protein interactions, cancer-related somatic mutations and pathway databases
    with open(output_tsv_full_path, 'w') as output_tsv:
        w = csv.writer(output_tsv, delimiter='\t')
        peak_header_new = peak_header + peak_header_insert
        w.writerow(peak_header_new)
        for peak_array_row in peak_array:
            w.writerow(peak_array_row)

    print('All processes finished! The output file is: {}'.format(output_tsv_full_path))
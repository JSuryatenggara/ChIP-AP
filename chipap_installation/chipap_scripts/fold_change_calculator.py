#!/usr/bin/env python3
#pyright: reportUnboundVariable=false


script_version = '3.1'


# Processes the peak list file "dataset_name_all_peaks_annotated" by:
#   - Adding the ChIP read depth (a.k.a. ChIP tag count) at each weighted peak center
#   - Adding the background control read depth (a.k.a. Control tag count) at each weighted peak center
#   - Adding the fold enrichment (a.k.a Fold Change) at each weighted peak center
#   - Adding the number of peak caller overlaps of each peak
#   - Adding the the number of overlapped DNA-protein binding motif (if provided by user)
#   - Renaming the number of motifs column
#   - Renaming the peak caller combinations column

# INPUT     - "dataset_name_all_peaks_annotated" - gene annotated complete peaks list 
#           - Chromatin IP BAM files (+ replicates if available) - to calculate the peak depths at each peak coordinate
#           - Background control BAM files (+ replicates if available) - to calculate the peak depths at each peak coordinate

# OUTPUT    - "dataset_name_all_peaks_calculated" - peak stats calculated gene annotated complete peaks list



# PATCH NOTES
#   Version 2.1     Added .fillna(0) to the line peak_df['Entrez ID'] = peak_df['Entrez ID'].astype('Int64')
#                       so that all Entrez ID values can be coverted into integers
#
#   Version 3.0     Now appends weighted peak center coordinates into the output peak list
#                       Just like tag counts and fold change, this weighted peak center is replicate-wise
#                       Due to the absence of weighted peak center in broad peaks, weighted peak center is replaced by 
#                           simply the midpoint between the broad peak's start and end coordinates
#
#   Version 3.1     Generates narrowPeak- or broadPeak-formatted output for downstream peak IDR calculation
#                   "Strand" column now contains all "." values instead of all "+" values to better represent the strandless nature of the peaks


print('Importing required modules')
import pandas as pd
import subprocess
import multiprocessing
import os
import argparse



print('Defining required functions')

# Function to parse the absolute path, basename, and extension of the assigned argument: 

# INPUT     - file_relative_path_list - a list of relative paths to n files

# OUTPUT    - absolute_path - a list of the absolute paths of file_relative_path_list elements
#           - original_name - a list of the file basenames of file_relative_path_list elements
#           - original_extension a list of the file extensions of file_relative_path_list elements

def file_basename_parsing(file_relative_path_list):
    absolute_path = [os.path.abspath(relative_path_element) for relative_path_element in file_relative_path_list]
        
    original_name = []
    original_extension = []

    for absolute_path_element in absolute_path:
        current_extension = '.bam'
        current_file_name = absolute_path_element.split('/')[-1].strip(current_extension)
        original_name.append(current_file_name)
        original_extension.append(current_extension)
        
    return absolute_path, original_name, original_extension



# Function to measure the number of reads at ChIP and control weighted peak center;
#   also calculates the ChIP vs control fold change

# INPUT     - chip_calculator_arg - the current peak coordinate in the format of chr:start-end, as listed in variable Peak_ID
#           - current_chip_bam - the absolute path to the Chromatin IP BAM file currently being calculated
#           - current_ctrl_bam - the absolute path to the Background control BAM file currently being calculated

# OUTPUT    - chip_weighted_peak_center_read_count - the read depth at the current weighted peak center within the current ChIP BAM file
#           - ctrl_weighted_peak_center_read_count - the read depth at the current weighted peak center within the current Control BAM file
#           - weighted_peak_center_fold_change - the fold enrichment value at the current weighted peak center within the current ChIP and Control BAM files

def fold_change_calculator_function_narrow(chip_calculator_arg):

    # Samtools depth run on a .bam file returns a list of chromosomal locations (single base resolution) 
    #   paired with the read depth (how many mapped reads are overlapping, similar to read counts) at those locations, 
    #   then temporarily store the list in form of text file named after its chromosomal location code (chr:start-end)
    subprocess.run('samtools depth -aa -r {} {} > temp_{}/{}.reads'.format(chip_calculator_arg, current_chip_bam, chip_bam_name[chip_replicate_counter], chip_calculator_arg), shell = True)

    # Read the list in samtools depth output file above
    chip_depth_df           = pd.read_csv('temp_{}/{}.reads'.format(chip_bam_name[chip_replicate_counter], chip_calculator_arg), delimiter = '\t', header = None)
    chip_depth_array        = chip_depth_df.values.tolist()

    # Generate a list of read depths that spans from the start to the end of the peak
    chip_depth_count_list   = [chip_depth_array_row[2] for chip_depth_array_row in chip_depth_array]

    # Sum up the number of read depths within the peak
    chip_depth_total_count  = sum(chip_depth_count_list) 

    # Determine the median of the read depths sum, which points out to the single base location of the weighted peak center
    chip_depth_median_count = int(0.5 * chip_depth_total_count)
    accumulated_read_count  = 0



    # Iterates the list starting from the leftmost base
    for chip_depth_array_row in chip_depth_array:

        # While keep adding up the read depths from each iteration to the accumulated_read_count
        accumulated_read_count += chip_depth_array_row[2]

        # When accumulated_read_count went over the median of the read depths sum, it has reached the weighted peak center position
        if accumulated_read_count >= chip_depth_median_count:

            # Record the weighted peak center position later to be used to get the read depth of the control sample
            chip_weighted_peak_center_chromosome = chip_depth_array_row[0]

            # Record the weighted peak center position later to be used to get the read depth of the control sample 
            chip_weighted_peak_center_coordinate = chip_depth_array_row[1]

            # Record the ChIP read depth at weighted peak center
            chip_weighted_peak_center_read_count = chip_depth_array_row[2]

            # Normalize ChIP read depth using the assigned normalization factor
            chip_weighted_peak_center_read_count = chip_weighted_peak_center_read_count / chip_normalization_factor_value

            # No use iterating the half remaining peak, break out of the loop
            break



    # Assign the chromosomal location code for getting the read depth by means of samtools view
    chip_weighted_peak_center_ID = '{}:{}-{}'.format(chip_weighted_peak_center_chromosome, chip_weighted_peak_center_coordinate, chip_weighted_peak_center_coordinate)
    
    # Call samtools depth on the control .bam file using one base chromosomal location (start = end);
    #   will return a line containing the position and the read depth at that exact one base position
    popen_ctrl = subprocess.Popen('samtools depth -aa -r {} {}'.format(chip_weighted_peak_center_ID, current_ctrl_bam), shell = True, stdout = subprocess.PIPE)

    # Convert the byte-type parsed output of samtools depth
    ctrl_out = popen_ctrl.communicate()[0].decode("utf-8").split()

    # Get the read depth number we want, which is at the end of the list
    if len(ctrl_out) == 3:
        ctrl_weighted_peak_center_read_count = int(ctrl_out[-1])
    else:
        ctrl_weighted_peak_center_read_count = 1
            
    # Normalize control read depth using the assigned normalization factor
    ctrl_weighted_peak_center_read_count = ctrl_weighted_peak_center_read_count / ctrl_normalization_factor_value



    # Calculate the fold change by simply dividing the read depth of the ChIP sample by the read depth of the control sample
    # If control read depth less than 1, assume that it is 1 (to avoid extreme inflation of 
    #   fold change value due to extremely small (or zero) denominator)
    
    if ctrl_weighted_peak_center_read_count >= 1: 
        weighted_peak_center_fold_change = chip_weighted_peak_center_read_count / ctrl_weighted_peak_center_read_count

    else: 
        weighted_peak_center_fold_change = chip_weighted_peak_center_read_count

    # Return the values of interest
    return round(chip_weighted_peak_center_read_count, 2), round(ctrl_weighted_peak_center_read_count, 2), round(weighted_peak_center_fold_change, 2), round(chip_weighted_peak_center_coordinate)




def fold_change_calculator_function_broad(chip_calculator_arg):

    # Call samtools view on the ChIP .bam file using the Peak ID (chr:start-end);
    #   will return the number of all reads found between the start and end positions
    popen_chip = subprocess.Popen('samtools view -c {} {}'.format(current_chip_bam, chip_calculator_arg), shell = True, stdout = subprocess.PIPE)
    chip_out = popen_chip.communicate()[0] # Get the number of all reads
    chip_tag_count = int(chip_out.strip()) # Remove all whitespaces, so can be converted into integer
    
    # Call samtools view on the control .bam file using the Peak ID;
    #   will return the number of all reads found between the start and end positions
    popen_ctrl = subprocess.Popen('samtools view -c {} {}'.format(current_ctrl_bam, chip_calculator_arg), shell = True, stdout = subprocess.PIPE)
    ctrl_out = popen_ctrl.communicate()[0] # Get the number of all reads
    ctrl_tag_count = int(ctrl_out.strip()) # Remove all whitespaces, so can be converted into integer
    
    # Normalize ChIP read depth using the assigned normalization factor
    chip_tag_count = chip_tag_count / chip_normalization_factor_value

    # Normalize control read depth using the assigned normalization factor
    ctrl_tag_count = ctrl_tag_count / ctrl_normalization_factor_value



    # Calculate the fold change by simply dividing the read depth of the ChIP sample by the read depth of the control sample
    # If control read depth less than 1, assume that it is 1 (to avoid extreme inflation of 
    #   fold change value due to extremely small (or zero) denominator)
    
    if ctrl_tag_count >= 1: 
        average_fold_change = chip_tag_count / ctrl_tag_count

    else: 
        average_fold_change = chip_tag_count

    # Return the values of interest
    return round(chip_tag_count, 2), round(ctrl_tag_count, 2), round(average_fold_change, 2)




# Function to generate Peak_ID out of chromosomal coordinates of the current peak
# INPUT     - chr_value     - the chromosome number at which the current peak is located
#           - start_value   - the smallest base number at which the current peak is located
#           - end_value     - the largest base number at which the current peak is located
# OUTPUT    - the current Peak_ID, made out of the current peak coordinates, in the format of chr:start-end
def peak_ID_generator_function(chr_value, start_value, end_value):
    return '{}:{}-{}'.format(chr_value, start_value, end_value)



# Function to split a pipe-separated string, converting it into a list, then count the number of elements
# INPUT     - peak_caller_combination - pipe-separated string containing name of the peak callers that called for the current peak
# OUTPUT    - the number of the peak callers names contained in the peak_caller_combination string
def overlap_number_calculator_function(peak_caller_combination):
    return len(peak_caller_combination.split('|'))



# Function to only output the input values before the brackets
# INPUT     - annotation - the annotation value of the current peak in the default HOMER format
# OUTPUT    - simplified annotation: without the brackets and bracketed values that comes after
def annotation_simplifier_function(annotation):
    if isinstance(annotation, float):
        return 'None'
    else:
        return annotation.split('(')[0]



print('Setting argument parser')

parser = argparse.ArgumentParser()

parser.add_argument('--thread', 
                    help = '<Optional> Maximum number of processes to use. Default is half the maximum available.', 
                    type = int, 
                    choices = range(1, (multiprocessing.cpu_count() + 1), 1),
                    metavar = "[1-{}]".format(multiprocessing.cpu_count()),
                    default = int(multiprocessing.cpu_count() / 2))

parser.add_argument('--peak', 
                    help = '<Optional> Peak type. Narrow peaks for transcription factors (default). Broad peaks for histone modifiers.',
                    choices = ['narrow', 'broad'],
                    default = 'narrow')

parser.add_argument('--input_tsv', 
                    help = '<Required> Input peak list file, in HOMER annotatePeaks format, .tsv extension.', 
                    required = True)

parser.add_argument('--output_tsv', 
                    help = '<Required> Output peak list file, in HOMER annotatePeaks format (modified), .tsv extension', 
                    required = True)

parser.add_argument('--chip_bam', 
                    nargs = '+', help = '<Required> The ChIP dataset aligned reads (.bam) file: ordered by replicate number, separated by space.', 
                    required = True)

parser.add_argument('--ctrl_bam', 
                    nargs = '+', help = '<Required> The control dataset aligned reads (.bam) file: ordered by replicate number, separated by space.', 
                    required = True)

parser.add_argument('--normfactor', 
                    help = '<Optional> Assign the factor to use for read counts normalization', 
                    choices = ['mapped', 'uniquely_mapped', 'user_value'])

parser.add_argument('--chip_norm', 
                    nargs = '+', 
                    help = '<Required if using --normfactor user_value>, Assign the custom user values for the ChIP replicates: order by replicate number, separated by space', 
                    type = int)

parser.add_argument('--ctrl_norm', 
                    nargs = '+', 
                    help = '<Required if using --normfactor user_value>, Assign the custom user values for the control replicates: order by replicate number, separated by space', 
                    type = int)

args = parser.parse_args()



subprocess.run('ulimit -n 2000', shell = True)

print('Parsing arguments')
cpu_count = args.thread

peak_type = args.peak # Experiment-dependent peak type. Narrow for TF peaks. Broad for histone modifier peaks.

input_tsv_full_path = os.path.abspath(args.input_tsv)
output_tsv_full_path = os.path.abspath(args.output_tsv)

# Parsing the absolute path, basename, and extension of the ChIP and control .bam files:
chip_bam_absolute_path, chip_bam_name, chip_bam_extension = file_basename_parsing(args.chip_bam)
ctrl_bam_absolute_path, ctrl_bam_name, ctrl_bam_extension = file_basename_parsing(args.ctrl_bam)



print('Opening file: {}'.format(input_tsv_full_path)) # Reading the HOMER annotated peak file
peak_df = pd.read_csv(input_tsv_full_path, delimiter = '\t')
peak_df.sort_values(by = ['Chr', 'Start'], inplace = True) # Sorting the peaks list based on chromosomal position
peak_header = peak_df.columns.tolist()

peak_df['Peak ID'] = peak_df.apply(lambda peak_df: peak_ID_generator_function(peak_df['Chr'], peak_df['Start'], peak_df['End']), axis = 1)

# Basically renaming the column into something that makes more sense than what's given by HOMER
peak_df['Peak Caller Combination']  = peak_df['Focus Ratio/Region Size'] 

# Generating a list of number of peak caller overlaps based on the list length of
#   pipe-symbol-splitted-string of 'Peak Caller Combination' value
peak_df['Peak Caller Overlaps']     = peak_df.apply(lambda peak_df: overlap_number_calculator_function(peak_df['Peak Caller Combination']), axis = 1)

# Generating a list of simplified information regarding the region of the peak location, relative to the nearest gene
# The function is to get rid of the irrelevant values in the bracket
peak_df['Annotation']               = peak_df.apply(lambda peak_df: annotation_simplifier_function(peak_df['Annotation']), axis = 1)

peak_df['Strand']                   = '.'



if len(chip_bam_absolute_path) == 1 and len(ctrl_bam_absolute_path) == 1:
    # Assign column names for fold_change_calculator_function outputs in scenario where there are no multiple replicates
    calculator_output_column_name = ['ChIP Tag Count', 'Control Tag Count', 'Fold Change', 'Peak Center']

if len(chip_bam_absolute_path) > 1 and len(ctrl_bam_absolute_path) > 1:
    # Assign column names for fold_change_calculator_function outputs in scenario where there are multiple replicates
    calculator_output_column_name = ['ChIP Tag Count Rep {}'.format(chip_replicate_counter + 1) for chip_replicate_counter in range(len(chip_bam_absolute_path))] + ['Control Tag Count Rep {}'.format(chip_replicate_counter + 1) for chip_replicate_counter in range(len(chip_bam_absolute_path))] + ['Fold Change Rep {}'.format(chip_replicate_counter + 1) for chip_replicate_counter in range(len(chip_bam_absolute_path))] + ['Peak Center Rep {}'.format(chip_replicate_counter + 1) for chip_replicate_counter in range(len(chip_bam_absolute_path))]



for chip_replicate_counter in range(len(chip_bam_absolute_path)):
    print('Running the calculator on ChIP replicate {} and {}'.format(chip_bam_name[chip_replicate_counter], ctrl_bam_name[chip_replicate_counter]))

    current_chip_bam = chip_bam_absolute_path[chip_replicate_counter]
    current_ctrl_bam = ctrl_bam_absolute_path[chip_replicate_counter]

    # Default normalization factor
    chip_normalization_factor_value = 1
    # Default normalization factor
    ctrl_normalization_factor_value = 1



    if args.normfactor == 'mapped':
        popen_chip_mapped = subprocess.Popen('samtools view -F4 -c {}'.format(current_chip_bam), shell = True, stdout = subprocess.PIPE)
        # Get the number of mapped reads in the ChIP .bam file
        chip_mapped_out = popen_chip_mapped.communicate()[0]
        chip_mapped_read_count = int(chip_mapped_out.strip())
        print('chip_mapped_read_count:', chip_mapped_read_count)

        popen_ctrl_mapped = subprocess.Popen('samtools view -F4 -c {}'.format(current_ctrl_bam), shell = True, stdout = subprocess.PIPE)
        # Get the number of mapped reads in the control .bam file
        ctrl_mapped_out = popen_ctrl_mapped.communicate()[0]
        ctrl_mapped_read_count = int(ctrl_mapped_out.strip())
        print('ctrl_mapped_read_count:', ctrl_mapped_read_count)
        
        # Calculate the normalization factor for the ChIP sample read depths
        chip_normalization_factor_value = chip_mapped_read_count / ctrl_mapped_read_count
        ctrl_normalization_factor_value = 1



    if args.normfactor == 'uniquely_mapped':
        popen_chip_uniquely_mapped = subprocess.Popen('samtools view -F256 -c {}'.format(current_chip_bam), shell = True, stdout = subprocess.PIPE)
        # Get the number of uniquely mapped reads in the ChIP .bam file
        chip_uniquely_mapped_out = popen_chip_uniquely_mapped.communicate()[0]
        chip_uniquely_mapped_read_count = int(chip_uniquely_mapped_out.strip())
        print('chip_uniquely_mapped_read_count:', chip_uniquely_mapped_read_count)

        popen_ctrl_uniquely_mapped = subprocess.Popen('samtools view -F256 -c {}'.format(current_ctrl_bam), shell = True, stdout = subprocess.PIPE)
        # Get the number of uniquely mapped reads in the control .bam file
        ctrl_uniquely_mapped_out = popen_ctrl_uniquely_mapped.communicate()[0]
        ctrl_uniquely_mapped_read_count = int(ctrl_uniquely_mapped_out.strip())
        print('ctrl_uniquely_mapped_read_count:', ctrl_uniquely_mapped_read_count)

        # Calculate the normalization factor for the ChIP sample read depths
        chip_normalization_factor_value = chip_uniquely_mapped_read_count / ctrl_uniquely_mapped_read_count
        ctrl_normalization_factor_value = 1
    


    if args.normfactor == 'user_value':
        # Get the ChIP normalization factor directly from the user inputted argument
        chip_norm_factor = args.chip_norm
        # Get the control normalization factor directly from the user inputted argument
        ctrl_norm_factor = args.ctrl_norm

        chip_normalization_factor_value = chip_norm_factor[chip_replicate_counter]
        ctrl_normalization_factor_value = ctrl_norm_factor[chip_replicate_counter]

    print('chip_normalization_factor_value:', chip_normalization_factor_value)
    print('ctrl_normalization_factor_value:', ctrl_normalization_factor_value)

    # Getting the folder ready for the temporary samtools depth files generated in fold_change_calculator_function
    subprocess.run('mkdir -p temp_{}'.format(chip_bam_name[chip_replicate_counter]), shell = True)



    pool_calculator = multiprocessing.Pool(processes = cpu_count)
    
    # Takes in a list of each peak's chromosomal location codes (chr:start-end)
    # Reads the current_chip_bam and current_ctrl_bam assigned above at those locations
    # Returns ChIP read depth, control read depth, and fold change value at the weighted peak center of each peak
    # Only for transcription factor target protein sample (narrow peaks)
    if peak_type == 'narrow':
        chip_tag_count_list, ctrl_tag_count_list, fold_change_list, peak_center_list = zip(*pool_calculator.map(fold_change_calculator_function_narrow, peak_df['Peak ID']))
    

    # Takes in a list of each peak's chromosomal location codes (chr:start-end)
    # Reads the current_chip_bam and current_ctrl_bam assigned above at those locations
    # Returns ChIP read count, control read count, and average fold change value along the whole range of each peak
    # Only for histone modifier target protein sample (broad peaks)
    if peak_type == 'broad':
        chip_tag_count_list, ctrl_tag_count_list, fold_change_list = zip(*pool_calculator.map(fold_change_calculator_function_broad, peak_df['Peak ID']))
        peak_center_df = (0.5 * (peak_df['Start'] + peak_df['End'])).round(0)
        peak_center_list = peak_center_df.values.tolist()

    pool_calculator.close()
    pool_calculator.join()

    # Remove the samtools depth output temporary folder
    subprocess.run('rm -rf temp_{}'.format(chip_bam_name[chip_replicate_counter]), shell = True)



    if len(chip_bam_absolute_path) == 1 and len(ctrl_bam_absolute_path) == 1:
        # Save the final fold_change_calculator_function outputs in their appropriate pandas columns
        peak_df['ChIP Tag Count']       = chip_tag_count_list
        peak_df['Control Tag Count']    = ctrl_tag_count_list
        peak_df['Fold Change']          = fold_change_list
        peak_df['Peak Center']          = peak_center_list

    if len(chip_bam_absolute_path) > 1 and len(ctrl_bam_absolute_path) > 1:
        # Save the current fold_change_calculator_function outputs in their 
        #   appropriate pandas columns before going for another iteration
        peak_df['ChIP Tag Count Rep {}'.format(chip_replicate_counter + 1)]     = chip_tag_count_list
        peak_df['Control Tag Count Rep {}'.format(chip_replicate_counter + 1)]  = ctrl_tag_count_list
        peak_df['Fold Change Rep {}'.format(chip_replicate_counter + 1)]        = fold_change_list
        peak_df['Peak Center Rep {}'.format(chip_replicate_counter + 1)]        = peak_center_list



motif_number_column_exist = 0

# Parsing for the existence of number of motifs column. The column might be absent if 
#   -motif flag was not used during HOMER annotatePeaks.pl run upstream
for peak_header_counter in range(len(peak_header)):
    if 'Distance From' in peak_header[peak_header_counter]:
        peak_df['Number of Motifs'] = peak_df.iloc[:, peak_header_counter] # Rename the column into something that makes more sense
        motif_number_column_exist = 1
        break

if motif_number_column_exist == 0:
    peak_df['Number of Motifs'] = 0 # If the column does not exist, then the resulting output values will be all zeroes.



# Determining all the columns wanted in the output table, and their orders
peak_df = peak_df[['Peak ID', 
                    'Chr', 
                    'Start', 
                    'End', 
                    'Strand', 
                    'Peak Caller Combination', 
                    'Peak Caller Overlaps'
                    ] + calculator_output_column_name + [
                    'Number of Motifs', 
                    'Annotation', 
                    'Detailed Annotation', 
                    'Distance to TSS', 
                    'Nearest PromoterID', 
                    'Entrez ID', 
                    'Nearest Unigene', 
                    'Nearest Refseq', 
                    'Nearest Ensembl', 
                    'Gene Name', 
                    'Gene Alias', 
                    'Gene Description', 
                    'Gene Type', 
                    'CpG%', 
                    'GC%'
                    ]]

# Fix the bug where pandas automatically adds one decimal point to all the numbers in the 'Entrez ID' column
# Some entries in the 'Entrez ID' column are NaN. When pandas parse column containing numbers 
#   with one or more NaNs, it automatically assigns all values into float type
peak_df['Entrez ID'] = peak_df['Entrez ID'].fillna(0).astype('Int64')

print('Writing the result to file {}'.format(output_tsv_full_path))
# HOMER annotated peak list augmented with essential information regarding the 
#   'quality' of each peak (fold change, read depth, motif count)
peak_df.to_csv(output_tsv_full_path, sep = '\t', index = False)



# Creating "pseudo"-columns required to convert the fold change calculated peak list into .narrowPeak or .broadPeak
# Necessary so that the fold change calculated peak list (union peak set) can be processed by the IDR module
peak_df['score'] = 0 # Creating a new column "score" in peak_df, which values are all 0
peak_df['p-value'] = -1 # Creating a new column "p-value" in peak_df, which values are all -1
peak_df['q-value'] = -1 # Creating a new column "q-value" in peak_df, which values are all -1
peak_df['peak'] = -1 # Creating a new column "peak" in peak_df, which values are all -1 (only needed for .narrowPeak)

if len(chip_bam_absolute_path) == 1 and len(ctrl_bam_absolute_path) == 1: # If there are only a single "Fold Change" column in the fold change calculated peak list
    # Creating a new column "signalValue" in peak_df, which values are obtained from the "Fold Change" column values
    peak_df['signalValue'] = peak_df[['Fold Change']]

if len(chip_bam_absolute_path) > 1 and len(ctrl_bam_absolute_path) > 1: # If there are multi-replicated "Fold Change" columns in the fold change calculated peak list
    # Creating a new column "signalValue" in peak_df, which values are obtained from the average of "Fold Change" columns values
    peak_df['signalValue'] = peak_df[['Fold Change Rep {}'.format(chip_replicate_counter + 1) for chip_replicate_counter in range(len(chip_bam_absolute_path))]].mean(axis = 1)

# Prepare the peak_df to be a proper input for IDR module processing: ranking the peaks according to the number of peak caller overlaps and signalValue (fold change)
peak_df.sort_values(by = ['Peak Caller Overlaps', 'signalValue'], inplace = True, ascending = False)



# Make a .narrowPeak formatted peak list by getting all the necessary columns from the peak_df
if peak_type == 'narrow':
    narrowpeak_df = peak_df[['Chr',
                            'Start',
                            'End',
                            'Peak ID',
                            'score',
                            'Strand',
                            'signalValue',
                            'p-value',
                            'q-value',
                            'peak'
                            ]]

    # Save the .narrowPeak formatted peak list (the union peak set for IDR input) under .narrowPeak extension
    output_tsv_extension = '.' + output_tsv_full_path.split('.')[-1] if len(output_tsv_full_path.split('.')) > 1 else ''
    output_narrowpeak_full_path = output_tsv_full_path.strip(output_tsv_extension) + '.narrowPeak'
    print('Writing narrowPeak-formatted output for IDR calculation {}'.format(output_narrowpeak_full_path))
    narrowpeak_df.to_csv(output_narrowpeak_full_path, sep = '\t', index = False, header = None)



# Make a .broadPeak formatted peak list by getting all the necessary columns from the peak_df
if peak_type == 'broad':
    broadpeak_df = peak_df[['Chr',
                            'Start',
                            'End',
                            'Peak ID',
                            'score',
                            'Strand',
                            'signalValue',
                            'p-value',
                            'q-value'
                            ]]

    # Save the .broadPeak formatted peak list (the union peak set for IDR input) under .broadPeak extension
    output_tsv_extension = '.' + output_tsv_full_path.split('.')[-1] if len(output_tsv_full_path.split('.')) > 1 else ''
    output_broadpeak_full_path = output_tsv_full_path.strip(output_tsv_extension) + '.broadPeak'
    print('Writing broadPeak-formatted output for IDR calculation {}'.format(output_broadpeak_full_path))
    broadpeak_df.to_csv(output_broadpeak_full_path, sep = '\t', index = False, header = None)
#!/usr/bin/env python3

import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--fa', required = True)
parser.add_argument('--output', required = False)
args = parser.parse_args()

genome_ref_fa_full_path = os.path.abspath(args.fa)

if args.output:
    genome_ref_fa_dir = os.path.abspath(args.output) + '/'

else:
    genome_ref_fa_filename = genome_ref_fa_full_path.split('/')[-1]
    genome_ref_fa_dir = genome_ref_fa_full_path.strip(genome_ref_fa_filename)

with open(genome_ref_fa_full_path, 'r') as fa_file:

    while True:
        file_row = fa_file.readline()

        chr_value_list = []

        if file_row == '':
            break

        if file_row != '':
            
            if '>' in file_row.strip():
                
                if len(chr_value_list) == 0:
                    chr_value = file_row.strip().strip('>')
                    chromfa_file = open('{}{}.fa'.format(genome_ref_fa_dir, chr_value), 'w')
                    continue
                
                if len(chr_value_list) != 0:
                    chromfa_file.close()
                    chr_value = file_row.strip().strip('>')
                    chromfa_file = open('{}{}.fa'.format(genome_ref_fa_dir, chr_value), 'w')
                    continue

            if '>' not in file_row.strip():
                chromfa_file.write(file_row)

    chromfa_file.close()
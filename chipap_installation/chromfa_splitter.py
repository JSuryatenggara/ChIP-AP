#!/usr/bin/env python3

import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--fa', required = True)
parser.add_argument('--output', required = False)
parser.add_argument('--gem', required = False, action = 'store_true')
args = parser.parse_args()

genome_ref_fa_full_path = os.path.abspath(args.fa)

if args.output:
    genome_ref_fa_dir = os.path.abspath(args.output)

else:
    genome_ref_fa_filename = genome_ref_fa_full_path.split('/')[-1]
    genome_ref_fa_dir = '/'.join(genome_ref_fa_full_path.split('/')[:-1])

with open(genome_ref_fa_full_path, 'r') as fa_file:

    while True:
        file_row = fa_file.readline()

        if file_row == '':
            break

        if file_row != '':
            
            if '>' in file_row.strip():
                chr_value = file_row.strip().strip('>')
                if args.gem and not chr_value.lower().startswith('chr'):
                    chromfa_file = open('{}/chr{}.fa'.format(genome_ref_fa_dir, chr_value), 'w')
                else:
                    chromfa_file = open('{}/{}.fa'.format(genome_ref_fa_dir, chr_value), 'w')
                chromfa_file.write(file_row)

            if '>' not in file_row.strip():
                chromfa_file.write(file_row)

    chromfa_file.close()
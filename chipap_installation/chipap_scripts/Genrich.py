#!/usr/bin/env python3
#pyright: reportUnboundVariable=false


script_version = '1.0'


# Usage: ./Genrich  -t <file>  -o <file>  [optional arguments]
# Required arguments:
#   -t  <file>       Input SAM/BAM file(s) for experimental sample(s)
#   -o  <file>       Output peak file (in ENCODE narrowPeak format)
# Optional I/O arguments:
#   -c  <file>       Input SAM/BAM file(s) for control sample(s)
#   -f  <file>       Output bedgraph-ish file for p/q values
#   -k  <file>       Output bedgraph-ish file for pileups and p-values
#   -b  <file>       Output BED file for reads/fragments/intervals
#   -R  <file>       Output file for PCR duplicates (only with -r)
# Filtering options:
#   -r               Remove PCR duplicates
#   -e  <arg>        Comma-separated list of chromosomes to exclude
#   -E  <file>       Input BED file(s) of genomic regions to exclude
#   -m  <int>        Minimum MAPQ to keep an alignment (default = 0)
#   -s  <float>      Keep sec alns with AS >= bestAS - <float> (default = 0)
#   -y               Keep unpaired alignments (default = false)
#   -w  <int>        Keep unpaired alns, lengths changed to <int>
#   -x               Keep unpaired alns, lengths changed to paired avg
# Options for ATAC-seq:
#   -j               Use ATAC-seq mode (default = false)
#   -d  <int>        Expand cut sites to <int> bp (default = 100)
#   -D               Skip Tn5 adjustments of cut sites (default = false)
# Options for peak-calling (ignored when --adjustp is used):
#   -p  <float>      Maximum p-value (default = 0.01)
#   -q  <float>      Maximum q-value (FDR-adjusted p-value; default = 1)
#   -a  <float>      Minimum AUC for a peak (default = 200)
#   -l  <int>        Minimum length of a peak (default = 0)
#   -g  <int>        Maximum distance between signif. sites (default = 100)
# Other options:
#   -X               Skip peak-calling
#   -P               Call peaks directly from a log file (-f)
#   -z               Option to gzip-compress output(s)
#   -v               Option to print status updates/counts to stderr

################################################################################

import subprocess
import os
import multiprocessing
import argparse
import math
import numpy as np
import pandas as pd

subprocess.run('ulimit -n 2000', shell = True)

parser = argparse.ArgumentParser()

parser.add_argument('--mode', 
                    help = '<Required> Single-end or paired-end sequencing read.', 
                    required = True, 
                    choices = ['single', 'paired'])
parser.add_argument('--thread', 
                    help = '<Optional> Maximum number of threads/processes to use. Default is half the maximum available.', 
                    type = int, 
                    choices = range(1,61,1),
                    metavar = "[0-100]")
parser.add_argument('--adjustp', 
                    help = '<Optional> Alleviate Genrich behaviour in low reads datasets by auto-adjusting p-value threshold',
                    action = 'store_true')

parser.add_argument('-t', 
                    help = '<Required> <file> ChIP name-sorted BAM file: ordered by replicate number, separated by space.',
                    nargs = '+',
                    required = True)
parser.add_argument('-o', 
                    help = '<Required> <file> Output peak file (in ENCODE narrowPeak format).', 
                    required = True)

parser.add_argument('-c', 
                    help = '<Optional> <file> Control name-sorted BAM file: ordered by replicate number, separated by space.',
                    nargs = '+')
parser.add_argument('-f', 
                    help = '<Optional> <file> Output bedgraph-ish file for p/q values.')
parser.add_argument('-k', 
                    help = '<Optional> <file> Output bedgraph-ish file for pileups and p-values.')
parser.add_argument('-b', 
                    help = '<Optional> <file> Output BED file for reads/fragments/intervals.')
parser.add_argument('-R', 
                    help = '<Optional> <file> Output file for PCR duplicates (use only with -r).')

parser.add_argument('-r', 
                    help = '<Optional> Remove PCR duplicates.',
                    action = 'store_true')
parser.add_argument('-e', 
                    help = '<Optional> <arg> Comma-separated list of chromosomes to exclude.')
parser.add_argument('-E', 
                    help = '<Optional> <file> Input BED file(s) of genomic regions to exclude.')
parser.add_argument('-m', 
                    help = '<Optional> <int> Minimum MAPQ to keep an alignment (default = 0).',
                    type = int)
parser.add_argument('-s', 
                    help = '<Optional> <float> Keep sec alignments with AS >= bestAS - <float> (default = 0).',
                    type = float)
parser.add_argument('-y', 
                    help = '<Optional> Keep unpaired alignments (default = false).',
                    action = 'store_true')
parser.add_argument('-w', 
                    help = '<Optional> <int> Keep unpaired alignments, lengths changed to <int>.')
parser.add_argument('-x', 
                    help = '<Optional> Keep unpaired alignments, lengths changed to paired average.',
                    action = 'store_true')

parser.add_argument('-j', 
                    help = '<Optional> Use ATAC-seq mode (default = false).',
                    action = 'store_true')
parser.add_argument('-d', 
                    help = '<Optional> <int> Expand cut sites to <int> bp (default = 100).')
parser.add_argument('-D', 
                    help = '<Optional> Skip Tn5 adjustments of cut sites (default = false).',
                    action = 'store_true')

parser.add_argument('-p', 
                    help = '<Optional> <float> Maximum p-value (default = 0.01).')
parser.add_argument('-q', 
                    help = '<Optional> <float> Maximum q-value (FDR-adjusted p-value; default = 1).')
parser.add_argument('-a', 
                    help = '<Optional> <float> Minimum AUC for a peak (default = 200).')
parser.add_argument('-l', 
                    help = '<Optional> <int> Minimum length of a peak (default = 0).')
parser.add_argument('-g', 
                    help = '<Optional> <int> Maximum distance between signif. sites (default = 100).')

parser.add_argument('-X', 
                    help = '<Optional> Skip peak-calling.',
                    action = 'store_true')
parser.add_argument('-P', 
                    help = '<Optional> Call peaks directly from a log file (-f).',
                    action = 'store_true')
parser.add_argument('-z', 
                    help = '<Optional> Option to gzip-compress output(s).',
                    action = 'store_true')
parser.add_argument('-v', 
                    help = '<Optional> Option to print status updates/counts to stderr.',
                    action = 'store_true')

parser.add_argument('--stdout', 
                    help = '<Optional> <file> Standard output channel log file (stdout / 1>).')
parser.add_argument('--stderr', 
                    help = '<Optional> <file> Standard error channel log file (stderr / 2>).')

args = parser.parse_args()

################################################################################

read_mode = args.mode

chip_bam_list = [os.path.abspath(chip_bam) for chip_bam in args.t]
genrich_chip_string = ','.join(chip_bam_list)

if args.c:
    ctrl_bam_list = [os.path.abspath(ctrl_bam) for ctrl_bam in args.c]
    genrich_ctrl_string = ','.join(ctrl_bam_list)

output_dir = os.path.abspath(args.o)

if args.thread:
    cpu_count = args.thread
elif not args.thread:    
    cpu_count = multiprocessing.cpu_count() / 2

if args.adjustp:

    mapped_read_count_list = []

    for list_counter in range(len(chip_bam_list)):

        current_chip_bam = '{}'.format(chip_bam_list[list_counter])
        
        popen_mapped = subprocess.Popen('samtools view -F4 -c {} -@ {}'.format(current_chip_bam, cpu_count), shell = True, stdout = subprocess.PIPE)
        # Get the number of mapped reads in the ChIP .bam file
        mapped_out = popen_mapped.communicate()[0]
        mapped_read_count = int(mapped_out.strip())
        mapped_read_count_list.append(mapped_read_count)

        if read_mode == 'single':
            mapped_read_count_average = np.mean(mapped_read_count_list)
        if read_mode == 'paired':
            mapped_read_count_average = np.mean(mapped_read_count_list) * 2

        print('ChIP (replicate #{}) mapped read count: {} reads'.format(list_counter + 1, mapped_read_count))

    print('ChIP average mapped read count: {} reads'.format(mapped_read_count_average))

    genrich_neg_log_p_threshold = 2 ** (((-0.016 * (math.log(mapped_read_count_average, 2))) + 0.5) * math.log(mapped_read_count_average, 2))
    adjusted_genrich_neg_log_p_threshold = 0.5 * genrich_neg_log_p_threshold
    genrich_p_threshold = 10 ** (-1 * adjusted_genrich_neg_log_p_threshold)

    print('Genrich peak-calling calculated minimum -logp value threshold: {}'.format(genrich_neg_log_p_threshold))
    print('Genrich peak-calling adjusted minimum -logp value threshold: {}'.format(adjusted_genrich_neg_log_p_threshold))
    print('Genrich peak-calling calculated maximum p-value threshold: {}'.format(genrich_p_threshold))

################################################################################

genrich_arg_list = []
genrich_arg_list.append('-t {}'.format(genrich_chip_string))
genrich_arg_list.append('-o {}'.format(args.o))

if args.c:
    genrich_arg_list.append('-c {}'.format(genrich_ctrl_string))
if args.f:
    genrich_arg_list.append('-f {}'.format(args.f))
if args.k:
    genrich_arg_list.append('-k {}'.format(args.k))
if args.b:
    genrich_arg_list.append('-b {}'.format(args.b))
if args.R:
    genrich_arg_list.append('-R {}'.format(args.R))

if args.r:
    genrich_arg_list.append('-r')
if args.e:
    genrich_arg_list.append('-e {}'.format(args.e))
if args.E:
    genrich_arg_list.append('-E {}'.format(args.E))
if args.m:
    genrich_arg_list.append('-m {}'.format(args.m))
if args.s:
    genrich_arg_list.append('-s {}'.format(args.s))
if args.y or read_mode == 'single':
    genrich_arg_list.append('-y')
if args.w:
    genrich_arg_list.append('-w {}'.format(args.w))
if args.x:
    genrich_arg_list.append('-x')

if args.j:
    genrich_arg_list.append('-j')
if args.d:
    genrich_arg_list.append('-d {}'.format(args.d))
if args.D:
    genrich_arg_list.append('-D')

if not args.adjustp:
    if args.p:
        genrich_arg_list.append('-p {}'.format(args.p))
    if args.q:
        genrich_arg_list.append('-q {}'.format(args.q))
    if args.a:
        genrich_arg_list.append('-a {}'.format(args.a))
    if args.l:
        genrich_arg_list.append('-l {}'.format(args.l))
    if args.g:
        genrich_arg_list.append('-g {}'.format(args.g))

if args.adjustp:
    genrich_arg_list.append('-p {}'.format(genrich_p_threshold))

if args.X:
    genrich_arg_list.append('-X')
if args.P:
    genrich_arg_list.append('-P')
if args.z:
    genrich_arg_list.append('-z')
if args.v:
    genrich_arg_list.append('-v')

if args.stdout:
    genrich_arg_list.append('1>> {}'.format(args.stdout))
if args.stderr:
    genrich_arg_list.append('2>> {}'.format(args.stderr))


genrich_arg_string = ' '.join(genrich_arg_list)

################################################################################

subprocess.run('Genrich {}'.format(genrich_arg_string), shell = True)

################################################################################
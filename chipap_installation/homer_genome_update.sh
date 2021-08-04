#!/bin/bash

set -euxo pipefail

perl /home/chip-ap_user/anaconda3/envs/python3_chipap/share/homer/.//configureHomer.pl -install hg38
perl /home/chip-ap_user/anaconda3/envs/python3_chipap/share/homer/.//configureHomer.pl -install hg19
perl /home/chip-ap_user/anaconda3/envs/python3_chipap/share/homer/.//configureHomer.pl -install mm9
perl /home/chip-ap_user/anaconda3/envs/python3_chipap/share/homer/.//configureHomer.pl -install mm10
perl /home/chip-ap_user/anaconda3/envs/python3_chipap/share/homer/.//configureHomer.pl -install sacCer3
perl /home/chip-ap_user/anaconda3/envs/python3_chipap/share/homer/.//configureHomer.pl -install dm6
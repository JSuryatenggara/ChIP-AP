#!/bin/bash

set -euxo pipefail

perl ~/anaconda3/share/homer/.//configureHomer.pl -install hg38
perl ~/anaconda3/share/homer/.//configureHomer.pl -install hg19
perl ~/anaconda3/share/homer/.//configureHomer.pl -install mm9
perl ~/anaconda3/share/homer/.//configureHomer.pl -install mm10
perl ~/anaconda3/share/homer/.//configureHomer.pl -install sacCer3
perl ~/anaconda3/share/homer/.//configureHomer.pl -install dm6
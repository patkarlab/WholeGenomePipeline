!/usr/bin/bash

bamfile=$1
output_prefix=$2


source activate new_base

mosdepth -n -t 3 --fast-mode -m --by 500 ${output_prefix} ${bamfile}


conda deactivate

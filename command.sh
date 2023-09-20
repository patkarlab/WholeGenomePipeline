#!/usr/bin/env bash 

source activate new_base
./wgs.nf -entry WGS \
--sequences /home/diagnostics/pipelines/WholeGenomePipeline/sequences/ \
--input /home/diagnostics/pipelines/WholeGenomePipeline/sample_list.csv \
-resume -bg
conda deactivate 

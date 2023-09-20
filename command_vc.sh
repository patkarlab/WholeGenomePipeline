#!/usr/bin/env bash

source activate new_base
./wgs_vc.nf -entry WGS \
--sequences /home/diagnostics/pipelines/WholeGenomePipeline/sequences/ \
--input /home/diagnostics/pipelines/WholeGenomePipeline/sample_list.csv \
-resume -bg

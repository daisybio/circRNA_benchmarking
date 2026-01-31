#!/bin/bash
nextflow run /nfs/data3/CIRCEST/pipeline_new_benchmarking -resume -profile apptainer,keep_work -with-tower

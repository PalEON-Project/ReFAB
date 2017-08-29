#!/bin/bash -l

#$ -q "geo*"
#$ -pe omp 5
#$ -l h_rt=24:00:00
#$ -j y
#$ -t 1-3100

Rscript master_cross_validation.R $SGE_TASK_ID
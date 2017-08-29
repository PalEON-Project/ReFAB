#!/bin/bash -l

#$ -q "geo*"
#$ -pe smp 1
#$ -l h_rt=24:00:00
#$ -j y
#$ -o logs/
#$ -t 1-3100

Rscript master_cross_validation.R $SGE_TASK_ID
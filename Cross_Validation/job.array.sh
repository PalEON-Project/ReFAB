#!/bin/bash -l

#$ -q "geo*"
#$ -l h_rt=10:00:00
#$ -j y
#$ -o logs/
#$ -t 8785-8846

Rscript master_cross_validation.R $SGE_TASK_ID

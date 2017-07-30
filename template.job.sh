#!/bin/bash -l

#$ -l h_rt=24:00:00
#$ -pe smp 1
#$ -N @SITE@.raiho
#$ -q "geo*"

module load R/3.4.0

R CMD BATCH  master.@SITE@.R
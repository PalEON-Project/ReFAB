#!/bin/bash -l

#$ -l h_rt=08:00:00
#$ -pe smp 1
#$ -N raiho.master.run.@SITE@

module load R/3.4.0

R CMD BATCH  master.@SITE@.R
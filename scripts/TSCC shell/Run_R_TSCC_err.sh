#!/bin/bash
#PBS -q hotel
#PBS -N 16S_dada_err2
#PBS -l nodes=1:ppn=2
#PBS -l walltime=12:00:00
#PBS -o /projects/ps-shurinlab/users/cbwall/Yos_water_16S/output/outs/Rscript_out
#PBS -e /projects/ps-shurinlab/users/cbwall/Yos_water_16S/output/errs/Rscript_err
#PBS -V
#PBS -M cbwall@ucsd.edu
#PBS -m ae
#PBS -A shurin-group

export MODULEPATH=/projects/builder-group/jpg/modulefiles/applications:$MODULEPATH
module load R/4.1.2
Rscript --no-save /projects/ps-shurinlab/users/cbwall/scripts/16S_dada_err2.R


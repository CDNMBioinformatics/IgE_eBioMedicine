#!/bin/bash
#$ -S /bin/bash
#$ -cwd
if [ -e /modules/tcl/init/bash ]
then
  source /modules/tcl/init/bash
fi
module load R/4.0.3
Rscript /udd/reprk/projects/TOPMed/scripts/CRA_850k_meffil_qc.R

Rscript --version

#R CMD BATCH /udd/reprk/projects/TOPMed/scripts/robust_regression_adapted_lme.R
#Rscript /udd/reprk/projects/TOPMed/scripts/robust_regression_adapted_lme.R
# qsub -l h_vmem=100111M -l lx6 combine_csv.sh


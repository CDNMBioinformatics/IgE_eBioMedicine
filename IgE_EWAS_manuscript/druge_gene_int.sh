#!/bin/bash
#$ -S /bin/bash
#$ -cwd
if [ -e /modules/tcl/init/bash ]
then
  source /modules/tcl/init/bash
fi
module load R/4.0.3
Rscript /udd/reprk/projects/TOPMed/IgE_EWAS_manuscript/druge_gene_int.R

Rscript --version

# qsub -l h_vmem=100111M -l lx7 test.sh

#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-12:00
#SBATCH -p short
#SBATCH -o Plasma.out
#SBATCH -e Plasma.err
#SBATCH --mem=200G
#SBATCH --mail-type=ALL
#SBATCH --job-name=Plasma_fine_typing

module load gcc/9.2.0 R/4.1.2 python/3.10.11 geos/3.10.2 cmake/3.22.2 gdal/3.1.4 udunits/2.2.28 
source /home/mup728/jupytervenv/bin/activate

Rscript fine_typing_all_Plasma_cells.r


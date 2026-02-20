#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-24:00
#SBATCH -p priority
#SBATCH -o Epi.out
#SBATCH -e Epi.err
#SBATCH --mem=200G
#SBATCH --mail-type=ALL
#SBATCH --job-name=Epi_fine_typing

module load gcc/9.2.0 R/4.1.2 python/3.10.11 geos/3.10.2 cmake/3.22.2 gdal/3.1.4 udunits/2.2.28 
source /home/mup728/jupytervenv/bin/activate

Rscript fine_typing_all_Epi_cells.r


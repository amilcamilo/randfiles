#!/bin/bash
#SBATCH --job-name=loggle_stock_market
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=10G
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amilcar.camilo@upf.edu

module load R/4.0.0-foss-2020a
Rscript covidstonks.r

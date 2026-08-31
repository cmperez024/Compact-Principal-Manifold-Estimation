#!/bin/bash

#SBATCH -n 12
#SBATCH --mem-per-cpu 56G
#SBATCH --job-name="Sphere 12n 56G N/2"
#SBATCH -A genacc_q
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL

module load R/4.4.0
R CMD BATCH RScripts/sphere_job.R
#!/bin/bash

#SBATCH --time 3-00:00

#SBATCH --output="log/%A-%x-std-output.out"
#SBATCH --error="log/%A-%x-err-output.out"

#SBATCH --mail-user=daiki.tagami@hertford.ox.ac.uk
#SBATCH --mail-type=ALL

#SBATCH --mem=10G
#SBATCH --cpus-per-task=1

#SBATCH --job-name="underdominant-selection"

#SBATCH --cluster=swan
#SBATCH --partition=standard-statgen-cpu
#SBATCH --nodelist=smew01.cpu.stats.ox.ac.uk

source /homes/tagami/miniconda3/bin/activate simulation

slim -d seed=10 underdominant-selection.slim
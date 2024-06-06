#!/bin/bash
#SBATCH --job-name=dec_tND
#SBATCH --output=disp_dec_tND.out
#SBATCH --error=disp_dec_tND.err
#SBATCH --time=30-00:00:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=u6562250@anu.edu.au
#SBATCH --mail-type=ALL

# path to your home directory
HOME_DIR=/mnt/data/dayhoff/home/u6562250

# activate conda environment
source /opt/conda/bin/activate $HOME_DIR/.conda/envs/iaa-birds-dispersals

Rscript /mnt/data/dayhoff/home/u6562250/dispersals/bgb_DEC_tND.R

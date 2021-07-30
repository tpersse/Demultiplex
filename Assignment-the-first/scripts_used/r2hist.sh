#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --output=histgen2%j.out
#SBATCH --error=histgen2%j.err
#SBATCH --mail-user='tpersse@uoregon.edu'
#SBATCH --mail-type=END,FAIL

conda activate bgmp_py39

/usr/bin/time -v python demux_hist.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -o /projects/bgmp/tpersse/bioinformatics/Bi622/demux/hist_read3.png
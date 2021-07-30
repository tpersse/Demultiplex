#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --output=length2%j.out
#SBATCH --error=length%j.err
#SBATCH --mail-user='tpersse@uoregon.edu'
#SBATCH --mail-type=END,FAIL

zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz | wc -l >> /projects/bgmp/tpersse/bioinformatics/Bi622/demux/length.txt
#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --output=demuxtest%j.out
#SBATCH --error=demuxtest%j.err
#SBATCH --mail-user='tpersse@uoregon.edu'
#SBATCH --mail-type=END,FAIL

conda activate bgmp_py39

/usr/bin/time -v python demux.py -Q 30 --index1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz --index2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz --read1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz --read2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz
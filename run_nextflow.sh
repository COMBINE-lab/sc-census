#!/bin/bash
#SBATCH --job-name=nextflow        # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=2        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=2G         # memory per cpu-core (4G is default)
#SBATCH --partition=cbcb       # total run time limit (HH:MM:SS)
#SBATCH --qos=highmem       # total run time limit (HH:MM:SS)
#SBATCH --account=cbcb       # total run time limit (HH:MM:SS)
#SBATCH --time=21-00:00:00       # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=dhe17@umd.edu

PROJDIR="/fs/nexus-projects/sc_read_census/nextflow"
cd $PROJDIR

$PROJDIR/nextflow $PROJDIR/main.nf -resume

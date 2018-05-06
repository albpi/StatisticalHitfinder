#!/bin/sh
#This is an example script for running SLURM jobs
#
# Set the name of the job
#SBATCH -J cxidb-78
#   
#
##SBATCH --exclusive --nodes=2
#SBATCH --nodes=4
# Ask for 2 tasks(processes)
#SBATCH --ntasks-per-node=8
##SBATCH --ntasks=48
#
#
# Choose what resources each task will use.
# We will ask for 1 CPU per task
#SBATCH --cpus-per-task=1
#
# The partition to use.

##SBATCH -p regular --mem-per-cpu=4000M

# 'A' partition
#SBATCH -p regular -x a015 -C GTX680 --mem-per-cpu=3000M

# 'C' partition
##SBATCH -p regular -C TitanBlack --mem-per-cpu=4000M 

##--mem-per-cpu=5000M

# Run a script using MPI with. The number
# of MPI processes started will match the
# requested number of tasks

#mpirun parallel.py
#mpirun --mca btl ^openib /home/alberto/StatisticalHitfinder/src/main_analysis_back.py
#mpirun --mca btl ^openib /home/alberto/StatisticalHitfinder/src/main_analysis_front_1.py
mpirun -v -tag-output /home/alberto/Statisticalhitfinder/src/hitfinding.py

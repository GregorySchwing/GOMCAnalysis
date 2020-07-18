#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jpotoff@wayne.edu
#SBATCH --time=100:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --exclude=potoff7,potoff40,potoff43
Running on host `hostname`
echo time is `date`
cd /home4/jpotoff/GOMC/SPCE/1r5
echo Directory is `pwd`
./GOMC_CPU_GCMC in.conf >out5a.dat

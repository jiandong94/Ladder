#!/bin/bash
#PBS -N U40_fourth
#PBS -l nodes=1:ppn=8
#PBS -l walltime=200:00:00
#PBS -l mem=20G
#PBS -j oe
#PBS -k oe
#PBS -V

APP=ladder_d4
INPUTFILE=IF_RU10U40_2


# Number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes
#export OMP_NUM_THREADS=8

# Creat workdir
ID=`echo "$PBS_JOBID" | grep -o "^[0-9]*"`
WORKDIR=$PBS_O_WORKDIR
echo Working directory is $WORKDIR
cd $WORKDIR

WN=`sed -n "/WriteNum\s*=\s*[0-9]*/p" $INPUTFILE | grep -o "[0-9]*"`

date
./$APP $INPUTFILE  2>&1 | tee o${ID}.${PBS_JOBNAME}_$WN
date


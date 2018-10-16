#!/bin/bash
#PBS -N test
#PBS -l nodes=1:ppn=8
#PBS -l walltime=200:00:00
#PBS -l mem=10G
#PBS -j oe
#PBS -k oe
#PBS -V

APP=ladder_d4
INPUTFILE=inputfile_ladder_test


# Number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo $PBS_NODEFILE
echo This job has allocated $NPROCS nodes
#export OMP_NUM_THREADS=8

# Creat workdir
ID=`echo "$PBS_JOBID" | grep -o "^[0-9]*"`
WORKDIR=$PBS_O_WORKDIR/$ID\_$PBS_JOBNAME
echo Working directory is $WORKDIR
mkdir $WORKDIR
cd $WORKDIR

# Copy inputfile to workdir
cp ../$INPUTFILE ./

date
../$APP $INPUTFILE 2>&1 | tee o$ID\.$PBS_JOBNAME
date

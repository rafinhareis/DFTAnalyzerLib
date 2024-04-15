#!/bin/bash
#PBS -S /bin/bash
#PBS -j oe
#PBS -l walltime=72:00:00
#PBS -l nodes=n01.fisica.ufmg.br:ppn=10+n02.fisica.ufmg.br:ppn=10
#PBS -N test
#PBS -q opmp
###

source /home/codes/intel/oneapi/setvars.sh

cd $PBS_O_WORKDIR
echo "PWD: " $PWD
#### run
mpirun -np 20 /home/codes/intel/programs/vasp.5.4.4.pl2/bin/vasp_std > test.out

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "Seu calculo terminou em:" $dt

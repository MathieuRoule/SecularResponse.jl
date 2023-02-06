#!/bin/bash
#PBS -S /bin/sh
#PBS -N Sellwood
#PBS -o ./log/Sellwood.log
#PBS -e ./log/Sellwood.err
#PBS -l nodes=1:ppn=128,walltime=24:00:00

module purge
module load julia/1.7.2

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export JULIA_NUM_THREADS=120
export JULIA_CPU_THREADS=1

JULIA=julia
PROJECT=/home/roule/Disc_2D/SecularResponse
CODE=./runsecular.jl

PREFIX=/home/roule/Disc_2D/SecularResponse/examples/Sellwood
cd ${PREFIX}

${JULIA} --project=${PROJECT} --threads ${JULIA_NUM_THREADS} ${CODE}

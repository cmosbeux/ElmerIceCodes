#!/bin/bash
#SBATCH -J MESH_OPTIM

#SBATCH --nodes=1
#SBATCH --constraint=HSW24

#SBATCH --ntasks=24

#SBATCH --ntasks-per-node=24
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1

#SBATCH --time=05:00:00

#SBATCH --output MESH_OPTIM.%j.output
#SBATCH --error  MESH_OPTIM.%j.error

#SBATCH --exclusive

ulimit -s unlimited

echo MESH_OPTIM_Antarctica.sif > ELMERSOLVER_STARTINFO


srun --mpi=pmi2 -K1 -n 1 ElmerSolver 



#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem-per-cpu=3850
#SBATCH --time=23:59:00
#SBATCH --account=su006-030

module purge
module load GCC/10.3.0 OpenMPI/4.1.1
module load OpenFOAM/v2106
source $FOAM_BASH

#./Allrun
decomposePar > log.decomposePar -allRegions 2>&1
mpirun  icoBoilingFoam -parallel > log.icoBoilingFoam 2>&1
reconstructPar > log.reconstructPar -allRegions 2>&1

#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128
#SBATCH --mem-per-cpu=3850
#SBATCH --time=47:59:00
#SBATCH --account=su006-014

module purge
module load GCC/10.3.0 OpenMPI/4.1.1
module load OpenFOAM/v2106
source $FOAM_BASH

#./Allclean
./Allrun
decomposePar > log.decomposePar -allRegions 2>&1
mpirun  icoBoilingFoam -parallel > log.icoBoilingFoam 2>&1
reconstructPar > log.reconstructPar -allRegions 2>&1

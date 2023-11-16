#!/bin/bash

#SBATCH --job-name=chtOneBubble
#SBATCH --nodes=2
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=47:00:00

# Replace [budget code] below with your project code (e.g. t01)
#SBATCH --account=e702 
#SBATCH --partition=standard
#SBATCH --qos=long

# Setup the job environment (this module needs to be loaded before any other modules)
module load epcc-job-env

# Load the appropriate modules and source the OpenFOAM bashrc file
# The first line makes PrgEnv-gnu available on the back end nodes.

module -s restore /etc/cray-pe.d/PrgEnv-gnu
module load openfoam/com/v2106
export HOME='/work/e702/e702/mmagnini'

source ${FOAM_INSTALL_DIR}/etc/bashrc

# Run OpenFOAM work
./Allrun
decomposePar > log.decomposePar -allRegions 2>&1
srun --distribution=block:block --hint=nomultithread icoBoilingFoam -parallel
reconstructPar > log.reconstructPar -allRegions -latestTime 2>&1 
#rm -r processor*

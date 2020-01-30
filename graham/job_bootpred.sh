#!/bin/bash
#SBATCH --account=def-bayne     # replace this with your own account
#SBATCH --nodes=1               # number of whole nodes
#SBATCH --ntasks-per-node=32    # 32 cores on each node
#SBATCH --mem-per-cpu=8G        # use --mem=0 for all ~3.9G mem per core
#SBATCH --time=12:00:00         # time (HH:MM:SS)
#SBATCH --job-name=bootpred
#SBATCH --output=%x-%j.out
#SBATCH --mail-user=solymos@ualberta.ca
#SBATCH --mail-type=ALL

# Load modules
module nixpkgs/16.09
module load gcc/7.3.0
module load gdal/3.0.1
module load proj4-fortran/1.0
module load openmpi/3.1.2
module load r/3.5.1

# Export the nodes names.
# If all processes are allocated on the same node,
# NODESLIST contains : node1 node1 node1 node1
export NODESLIST=$(echo $(srun hostname))

# Run R script
Rscript --vanilla bootpred.R

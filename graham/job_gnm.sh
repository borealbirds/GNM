#!/bin/bash
#SBATCH --account=def-psolymos  # replace this with your own account
#SBATCH --nodes=2               # number of whole nodes
#SBATCH --ntasks-per-node=32    # 32 cores on each node
#SBATCH --mem=0                 # use all ~3.9G mem per core
#SBATCH --time=12:00:00         # time (HH:MM:SS)
#SBATCH --job-name=north
#SBATCH --output=%x-%j.out
#SBATCH --mail-user=solymos@ualberta.ca
#SBATCH --mail-type=ALL

# Load modules
module nixpkgs/16.09
module load gcc/7.3.0
module load openmpi/3.1.2
module load r/3.5.1

# Export the nodes names.
# If all processes are allocated on the same node,
# NODESLIST contains : node1 node1 node1 node1
export NODESLIST=$(echo $(srun hostname))

# Run R script
Rscript --vanilla north.R

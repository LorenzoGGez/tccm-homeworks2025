#!/bin/bash


##SBATCH --job-name=multi_gpu_job
##SBATCH --time=08:00:00
##SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
##SBATCH --cpus-per-task=24
##SBATCH --gres=gpu:1)
##SBATCH --partition=g100_usr_interactive
##SBATCH --output=multiGPUJob.out
##SBATCH --error=multiGPUJob.err
##SBATCH --account=tra25_tccm

# Load necessary modules
module load nvhpc 

# Optional: Set environment variables for performance tuning
#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK   # Set OpenMP threads per task
#export NCCL_DEBUG=INFO                        # Enable NCCL debugging (for multi-GPU communication)

# Launch the distributed GPU application
# Replace with your actual command (e.g., mpirun or srun)
#srun --mpi=pmix ./my_distributed_gpu_app --config config.yaml

srun ./main < input > output

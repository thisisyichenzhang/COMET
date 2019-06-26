#!/bin/bash -l

#SBATCH --account=def-saram # adjust this to match the accounting group you are using to submit jobs
#SBATCH --time=0-5:00           # time (DD-HH:MM)
#SBATCH --mem=3G      			# memory; default unit is megabytes
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --output=%x-%j.out
#SBATCH --cpus-per-task=1      # adjust this if you are using parallel commands
#SBATCH --mail-user=zhangyichen93@gmail.com # adjust this to match your email address
#SBATCH --mail-type=ALL

# Choose a version of MATLAB by loading a module:
module load matlab/2018b
# Remove -singleCompThread below if you are using parallel commands:
srun matlab -nodisplay -nodesktop -nosplash -r "simmaster(1000, 100, 20, 10, 0, 1), exit"
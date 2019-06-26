#!/bin/bash

for n in 100
do
	for m in 15000;
	do  
		for d in 500;
		do
			for S_noise in 0 0.5 1 10;
			do
				for X_noise in 10;
				do 
					job_file="./n${n}_m${m}_d${d}_S${S_noise}_X${X_noise}.job"

			    echo "#!/bin/bash -l
#SBATCH --account=def-saram   # replace this with your own account
#SBATCH --time=02-00:00           # time (DD-HH:MM)
#SBATCH --mem=10G      			# memory; default unit is megabytes
#SBATCH --output=%x-%j.out
#SBATCH --mail-user=zhangyichen93@gmail.com # Send email updates to you or someone else
#SBATCH --mail-type=ALL          # send an email in all cases (job started, job ended, job aborted)
module load matlab/2018b   
srun matlab -nodisplay -nodesktop -nosplash -r 'simmaster($n, $m, $d, $X_noise, $S_noise, 100); exit' " > $job_file
    			sbatch $job_file
				done 


			done
	
		done

	done
    
done


for n in 1000
do
	for m in 15000;
	do  
		for d in 500;
		do
			for S_noise in 0 0.5 1 10;
			do
				for X_noise in 10;
				do 
					job_file="./n${n}_m${m}_d${d}_S${S_noise}_X${X_noise}.job"

			    echo "#!/bin/bash -l
#SBATCH --account=def-saram   # replace this with your own account
#SBATCH --time=02-00:00           # time (DD-HH:MM)
#SBATCH --mem=10G      			# memory; default unit is megabytes
#SBATCH --output=%x-%j.out
#SBATCH --mail-user=zhangyichen93@gmail.com # Send email updates to you or someone else
#SBATCH --mail-type=ALL          # send an email in all cases (job started, job ended, job aborted)
module load matlab/2018b   
srun matlab -nodisplay -nodesktop -nosplash -r 'simmaster($n, $m, $d, $X_noise, $S_noise, 100); exit' " > $job_file
    			sbatch $job_file
				done 


			done
	
		done

	done
    
done


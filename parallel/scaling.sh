#!/bin/bash

max_procs=16   # Maximum number of processors

procs=1
while [ $procs -le $max_procs ]; do
    # Submit the job with the current number of processors
    sbatch -A bsc32 -q gp_debug --ntasks-per-node=${procs} script.sh
    procs=$(( procs * 2 ))
done

# Exit after submitting all jobs
exit

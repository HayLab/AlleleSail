#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH --time=5:00:00   # walltime
#SBATCH --ntasks=8   # number of processor cores (i.e. tasks)
#SBATCH --nodes=2   # number of nodes
#SBATCH --mem-per-cpu=100MB   # memory per CPU core = mem PER TASK
#SBATCH -J "S gm mc congregate"   # job name
#SBATCH --mail-user=   # email address
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
fn_dir=alleleSail_data/1-11/SvI_germline_mc_mod/
smalldata_file=${fn_dir}last_gen_E_carrier.csv
num_runs=20

echo Count,Sex,fitCost,run,introFreq,totalPop >> $smalldata_file

for freq in $(seq 0 0.025 0.5)
do
    file_name=${fn_dir}introFreq_${freq}
    for fitness in $(seq -0.05 0.01 0.10)
    do
        for run in $(seq 0 $num_runs)
        do
            totalPop=`sed -E -n "s/[[:digit:]]*,([[:digit:]]*),50,Female,.*_${fitness}_${run}$/\1/p" \
                ${file_name}_total.csv`
            
            
            sed -n '/50,E,/p' ${file_name}_allele.csv | \
            sed -n "/_${fitness}_$run$/p" | \
            sed -n '/Female/p' | sed -n 1p | \
            sed -E -n "s/^.*,50,E,([[:digit:]|.]*),(\w*),.*_([[:digit:]|.|-]*)_([[:digit:]]*)$/\1,\2,\3,\4,$freq,$totalPop/p" \
            >> $smalldata_file

            totalPop=`sed -E -n "s/[[:digit:]]*,([[:digit:]]*),50,Male,.*_${fitness}_${run}$/\1/p" \
                ${file_name}_total.csv`
            
            sed -n '/50,E,/p' ${file_name}_allele.csv | \
            sed -n "/_${fitness}_$run$/p" | \
            sed -n '/Male/p' | sed -n 1p | \
            sed -E -n "s/^.*,50,E,([[:digit:]|.]*),(\w*),.*_([[:digit:]|.|-]*)_([[:digit:]]*)$/\1,\2,\3,\4,$freq,$totalPop/p" \
            >> $smalldata_file
        done
    done
done

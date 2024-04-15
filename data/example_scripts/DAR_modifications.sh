#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=12   # number of processor cores (i.e. tasks)
#SBATCH --nodes=12   # number of nodes
#SBATCH --mem-per-cpu=100MB   # memory per CPU core = mem PER TASK
#SBATCH -J "som_v_germ DAR 1"   # job name
#SBATCH --mail-user=  # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --partition=expansion

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
start_time=`date +%s`
echo RUNNING som_germ_MC_1.sh

num_runs=20
file_name=alleleSail_data/3-28/DAR

# start heatmap file
echo generation,sex,fitCost,run,introFreq >> $heatmap_file

# loop through introduction frequencies
for freq in 0.05 0.0 0.1
do
    echo starting fitnesses for freq $freq at `date`
    # loop through fitnesses to generate data
    for fitness in 0 0.05 0.1
    do
        # somatic MC
        srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_Ecost_${fitness}_recessive \
        -fn ${file_name}_somaticMC -r $num_runs -n 50 -i $freq \
        -fc 2 EE $fitness --MC --somatic & # ai is males, genotype, freq, first_gen, gap, repeats

        srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_Ecost_${fitness}_dominant \
        -fn ${file_name}_somaticMC -r $num_runs -n 50 -i $freq \
        -fc 2 E $fitness --dominant_FC --MC --somatic &

        # somatic 

        srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_Ecost_${fitness}_recessive \
        -fn ${file_name}_somatic -r $num_runs -n 50 -i $freq \
        -fc 2 EE $fitness --somatic & # ai is males, genotype, freq, first_gen, gap, repeats

        srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_Ecost_${fitness}_dominant \
        -fn ${file_name}_somatic -r $num_runs -n 50 -i $freq \
        -fc 2 E $fitness --dominant_FC --somatic &

        # germline 

        srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_Ecost_${fitness}_recessive \
        -fn ${file_name}_germline -r $num_runs -n 50 -i $freq \
        -fc 2 EE $fitness & # ai is males, genotype, freq, first_gen, gap, repeats

        srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_Ecost_${fitness}_dominant \
        -fn ${file_name}_germline -r $num_runs -n 50 -i $freq \
        -fc 2 E $fitness --dominant_FC &

        # germline MC

        srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_Ecost_${fitness}_recessive \
        -fn ${file_name}_germlineMC -r $num_runs -n 50 -i $freq \
        -fc 2 EE $fitness --MC & # ai is males, genotype, freq, first_gen, gap, repeats

        srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_Ecost_${fitness}_dominant \
        -fn ${file_name}_germlineMC -r $num_runs -n 50 -i $freq \
        -fc 2 E $fitness --dominant_FC --MC &
    done

    for fitness in 0 0.025 0.05
    do
        # somatic MC
        srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_Ecost_${fitness}_additive \
        -fn ${file_name}_somaticMC -r $num_runs -n 50 -i $freq \
        -fc 2 E $fitness --MC --somatic &

        # somatic
        srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_Ecost_${fitness}_additive \
        -fn ${file_name}_somatic -r $num_runs -n 50 -i $freq \
        -fc 2 E $fitness --somatic &

        # germline MC
        srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_Ecost_${fitness}_additive \
        -fn ${file_name}_germlineMC -r $num_runs -n 50 -i $freq \
        -fc 2 E $fitness --MC &

        # germlines
        srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_Ecost_${fitness}_additive \
        -fn ${file_name}_germline -r $num_runs -n 50 -i $freq \
        -fc 2 E $fitness &


    done
    wait
    echo finished fitnesses for freq $freq at `date`

done

end_time=`date +%s`
echo script time taken in seconds: 
echo $((end_time - start_time))
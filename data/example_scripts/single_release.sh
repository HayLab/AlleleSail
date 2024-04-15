#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH --time=12:00:00   # walltime
#SBATCH --ntasks=12   # number of processor cores (i.e. tasks)
#SBATCH --nodes=12   # number of nodes
#SBATCH --mem-per-cpu=100MB   # memory per CPU core
#SBATCH -J "tpc MC/noMC 3-5"   # job name
#SBATCH --mail-user=   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --partition=expansion

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
start_time=`date +%s`

source /home/mljohnso/zebrafish_project/my-test-venv_3.8/bin/activate
which python

## these ones are for new SM, correct alleles (!!!), full_five
# even tho there are actually 6. lol
echo RUNNING tpc_suppression_5_test.sh

l2=fitCosts

males=1
genotype=0
first_gen=0

num_runs=20
gap=1

filename=alleleSail_data/3-5/full_five_noMC

for repeats in 1
do
    for freq in 0.1 0.2
    do
        srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_gap_${gap}_repeats_${repeats}_XY \
        -fn $filename -sex XY -r $num_runs \
        -n 50 -i 0 -ai $males $genotype $freq $first_gen $gap $repeats &
        srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_gap_${gap}_repeats_${repeats}_ZW \
        -fn $filename -sex ZW -r $num_runs \
        -n 50 -i 0 -ai $males $genotype $freq $first_gen $gap $repeats &
        srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_gap_${gap}_repeats_${repeats}_fsRIDL \
        -fn $filename -a RW -r $num_runs \
        -n 50 -i 0 -ai $males $genotype $freq $first_gen $gap $repeats -fc 0 R 1.0 &
        srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_gap_${gap}_repeats_${repeats}_ZWviableBump \
        -fn $filename -sex ZW_viable -r $num_runs \
        -n 50 -i 0 -ai $males 1 $freq $first_gen $gap $repeats &
        srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_gap_${gap}_repeats_${repeats}_ZWviableNoBump \
        -fn $filename -sex ZW_viable -r $num_runs \
        -n 50 -i 0 -ai $males 0 $freq $first_gen $gap $repeats &
        srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_gap_${gap}_repeats_${repeats}_XYnoBump \
        -fn $filename -sex XY -r $num_runs \
        -n 50 -i 0 -ai $males 1 $freq $first_gen $gap $repeats &
    done
done

wait

end_time=`date +%s`
echo script time taken
echo $((end_time - start_time))
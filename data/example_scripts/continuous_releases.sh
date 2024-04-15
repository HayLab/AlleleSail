#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH --time=6:00:00   # walltime
#SBATCH --ntasks=16   # number of processor cores (i.e. tasks)
#SBATCH --nodes=16   # number of nodes
#SBATCH --mem-per-cpu=100MB   # memory per CPU core
#SBATCH -J "ais full noMC"   # job name
#SBATCH --mail-user=   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --partition=expansion

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
start_time=`date +%s`

## these simulations should be for sex determination, with additional releases
# DEFAULT fc = [] in this 
echo RUNNING ais_full.sh

filename=alleleSail_data/2-23/ais_51releases_noMC
l2=fitCosts

males=1
genotype=0
first_gen=0

num_runs=20
repeats=51
gap=1

for freq in 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5
do
    srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_rd_50_efficiency_1.0_repeats_${repeats}_XY \
    -fn $filename -sex XY -r $num_runs \
    -n 50 -i 0 -ai $males $genotype $freq $first_gen $gap $repeats &
    srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_rd_50_efficiency_1.0_repeats_${repeats}_ZW \
    -fn $filename -sex ZW -r $num_runs \
    -n 50 -i 0 -ai $males $genotype $freq $first_gen $gap $repeats &
    srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_rd_50_efficiency_1.0_repeats_${repeats}_fsRIDL \
    -fn $filename -a RW -r $num_runs \
    -n 50 -i 0 -ai $males $genotype $freq $first_gen $gap $repeats -fc 0 R 1.0 &

    srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_rd_50_efficiency_1.0_repeats_${repeats}_ZWviableBump \
    -fn $filename -sex ZW_viable -r $num_runs \
    -n 50 -i 0 -ai $males 1 $freq $first_gen $gap $repeats &
    srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_rd_50_efficiency_1.0_repeats_${repeats}_ZWviableNoBump \
    -fn $filename -sex ZW_viable -r $num_runs \
    -n 50 -i 0 -ai $males 0 $freq $first_gen $gap $repeats &
    srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_rd_50_efficiency_1.0_repeats_${repeats}_XYnoBump \
    -fn $filename -sex XY -r $num_runs \
    -n 50 -i 0 -ai $males 1 $freq $first_gen $gap $repeats &

    srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_rd_50_efficiency_1.0_repeats_${repeats}_XYsterile \
    -fn $filename -sex XYsterile -r $num_runs \
    -n 50 -i 0 -ai $males $genotype $freq $first_gen $gap $repeats -rd 50 50 50 -e 1.0 -e 1.0 -e 1.0 -e 1.0 -e 1.0 -e 1.0 &

    srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_rd_50_efficiency_1.0_repeats_${repeats}_ZWsterile \
    -fn $filename -sex ZWsterile -r $num_runs \
    -n 50 -i 0 -ai $males $genotype $freq $first_gen $gap $repeats -rd 50 50 50 -e 1.0 -e 1.0 -e 1.0 -e 1.0 -e 1.0 -e 1.0 &
done

wait

end_time=`date +%s`
echo script time taken
echo $((end_time - start_time))

#!/bin/bash
#Submit this script with: sbatch thefilename

#SBATCH --time=12:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=100MB   # memory per CPU core
#SBATCH -J "ais_low_efficiencies"   # job name
#SBATCH --mail-user=    # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
start_time=`date +%s`
module load python3/3.8.5

## these simulations should be for sex determination, with additional releases 
# DEFAULT fc = [] in this 
echo RUNNING ais_high_efficiencies.sh

filename=alleleSail_data/2-23/51releases_highEfficiencies_2
l2=fitCosts

males=1
genotype=0
first_gen=0

num_runs=20
repeats=51
gap=1

for freq in 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5
do
for efficiency in 0.8 0.9 0.95 0.99 1.0
do
    srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_rd_50_efficiency_${efficiency}_repeats_${repeats}_XY \
    -fn $filename --XY -r $num_runs \
    -n 50 -i 0 -ai $males $genotype $freq $first_gen $gap $repeats \
    -e $efficiency -e $efficiency &
    srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_rd_50_efficiency_${efficiency}_repeats_${repeats}_ZW \
    -fn $filename --ZW -r $num_runs \
    -n 50 -i 0 -ai $males $genotype $freq $first_gen $gap $repeats \
    -e $efficiency -e $efficiency &
    srun --ntasks 1 python3 alleleSail_pop_count.py -l introFreq_${freq}_rd_50_efficiency_${efficiency}_repeats_${repeats}_fsRIDL \
    -fn $filename -a RW -r $num_runs \
    -n 50 -i 0 -ai $males $genotype $freq $first_gen $gap $repeats -fc 0 R 1.0 \
    -e $efficiency -e $efficiency &
done
wait
done


end_time=`date +%s`
echo script time taken
echo $((end_time - start_time))

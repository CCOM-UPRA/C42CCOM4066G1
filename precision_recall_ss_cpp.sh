#!/bin/bash 
#SBATCH --job-name=precision_recall_ss
#SBATCH --mem 1
#SBATCH --mail-user=kevin.aviles6@upr.edu
#SBATCH --mail-type=ALL
#SBATCH --error=precision_recall_ss_cpp_error.txt
#SBATCH --output=precision_recall_ss_cpp_output.txt
#SBATCH -A trn035
#SBATCH -N 1
#SBATCH --time=02:00:00

# probe format: ORF1_100_mer_probes.fasta
# 02 refers to the type of compilation this is a more balanced approach
# You can test out different values and compare run times for an additional graph
g++ -std=c++17  -fopenmp -O2 precision_recall.cpp -o precision_recall_ex
# g++ -std=c++17 precision_recall.cpp -o precision_recall_ex
for ((k = 50; k <= 125; k += 25))
        do

        ./precision_recall_ex -g /ccsopen/home/kevinaviles6/genomes/GRCh38.fasta -f /ccsopen/home/kevinaviles6/pr_files_odo/orf12_probes/orf12_${k}mers.fasta -s orf12_sam/ORF12_${k}mers.sam -p precision_recall_results_${k}mers/ -e 5,$((k/2)),5 -t 100,200,25 -c /ccsopen/home/kevinaviles6/pr_files_odo/hum_csv/ -l 13671
        done


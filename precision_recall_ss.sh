#!/bin/bash 
#SBATCH --job-name=precision_recall_ss
#SBATCH --mem 1
#SBATCH --mail-user=kevin.aviles6@upr.edu
#SBATCH --mail-type=ALL
#SBATCH --error=precision_recall_ss_error.txt
#SBATCH --output=precision_recall_ss_output.txt
#SBATCH -A trn035
#SBATCH -N 1
#SBATCH --time=02:00:00

# probe format: ORF1_100_mer_probes.fasta
for ((k = 50; k <= 125; k += 25)) 
	do
	./precision_recall.sh -g /ccsopen/home/kevinaviles6/pr_files_odo/genomes/GRCh38.fasta -f /ccsopen/home/kevinaviles6/pr_files_odo/orf12_probes/orf12_${k}mers.fasta -s orf12_sam/ORF12_${k}mers.sam -p precision_recall_results_${k}mers/ -e 5,$((k/2)),5 -t 100,200,25 -c /ccsopen/home/kevinaviles6/pr_files_odo/hum_csv/ -l 13671
	done



#!/bin/bash
#SBATCH --job-name=75mers
#SBATCH --error=75mer_err.txt
#SBATCH --output=75mer.out
#SBATCH --mem 4000

#https://stackoverflow.com/questions/56962129/how-to-get-original-location-of-script-used-for-slurm-job
if [ -n $SLURM_JOB_ID ] ; then
	# We need to find L1PD files, but SLURM copies the script to a different
	# directory.  Change this directory accordingly when using SLURM.
	DIR="/work/jlopez2/emartinez/L1PD_files/"
else
	# Other L1PD scripts are in the same directory as this script.
	DIR=$(dirname $0)
fi

echo "variable DIR: $DIR"

#cd $PWD

usage() {
	cat << EOF
Usage:

Arguments:
	-g|--genome    Genome file.
	-f|--fasta     Fasta file with all the probes.
	-s|--sam       Directory where SAM files will be generated.
	-c|--meta      Directory with LINE-1 meta data file in .csv format.
	-l|--entries   Entries for the entire genome.
	-p|--print     Directory where results and raw files will be saved in.
	    				[Default: results/]
	-e|--edit      Vary the edit distance for MRFAST.
						[Default: 5,20,5]
	-t|--treshold  Vary threshold for Precision and Recall.
						[Default: 25,500,25]
	-m|--minkmer   Minimum amount of kmers to be used at once.
						[Default: 2]
	-h|--help      Help (shows usage)

		
EOF
}

#############################
# Initialize default values #
#############################

# Because many files are generated, store them in a subdirectory.
printdir="result"
minkmer=2

#####################
# Process arguments #
#####################

while [[ $# -gt 0 ]]; do
    case "$1" in
		-g|--genome)
			genome="$2"
			shift 2
			;;
		-f|--fasta)
			fasta="$2"
			shift 2
			;;
		-s|--sam)
			sam="$2"
			shift 2
			;;
		-c|--meta)
			meta="$2"
			shift 2
			;;
		-l|--entries)
			l1count="$2"
			shift 2
			;;
		-p|--print)
			printdir="$2"
			shift 2
			;;
		-e|--edit)
			edit_range="$2"
			IFS=',' read -r edit_start edit_end edit_step <<< "$edit_range"
			if [[ -z $edit_end || -z $edit_step ]]; then
				echo "Invalid edit range format: $edit_range"
				echo "Use a comma to separate start value, end value, and increment."
				echo "For example: 5,30,5"
				exit 1
			fi
			shift 2
			;;
		-t|--treshold)
			treshold_range="$2"
			IFS=',' read -r treshold_start treshold_end treshold_step <<< "$treshold_range"
			if [[ -z $treshold_end || -z $treshold_step ]]; then
				echo "Invalid treshold range format: $treshold_range"
				echo "Use a comma to separate start value, end value, and increment."
				echo "For example: 100,800,50"
				exit 1
			fi
			shift 2
			;;
		-m|--minkmer)
			minkmer="$2"
			shift 2
			;;
		-h|--help)
			usage
			exit 0
			;;
		*)
			echo "Invalid option: $1"
			exit 1
			;;
	esac
done

##############################################################################
# Initialize default values for edit_start and end and treshold_start and end#
##############################################################################

edit_start="${edit_start:-5}"
edit_end="${edit_end:-20}"
treshold_start="${treshold_start:-25}"
treshold_end="${treshold_end:-500}"
edit_step="${edit_step:-5}"
treshold_step="${treshold_step:-25}"

GENOME=$genome
# The following directory is where the SAM file is generated
SAMfile=$sam
FASTAfile=$fasta
DATA_DIR=$meta

MaxKmers=$(wc -l $FASTAfile | sed 's/ .*//')
MaxKmers=$((MaxKmers / 2))

L1Count=$l1count # <- For the entire genome (FLI+ORF2+FLnI)

pr_dir=$printdir
mkdir -p $pr_dir # -p ignores errors if the directory exists
# Vary the edit distance used by mrFAST
#for ((e = 5; e <= 30; e += 5))
for ((e = $edit_start; e <= $edit_end; e += $edit_step))
do
	if [ ! -f $GENOME.index ] ; then
		mrfast --index $GENOME
	fi
	time mrfast --search $GENOME --seq $FASTAfile -e $e -o $SAMfile
	# Vary the threshold to see how that affects P&R
	for ((t = $treshold_start; t <= $treshold_end ; t += $treshold_step))
	do
		# Create a CSV with tab-separated values
		printf "K-mers min\tPatterns\tTrue +\tFalse +\tFalse -\tPrecision\tRecall\tF1 Score\n" > $pr_dir/results_e${e}_t${t}.csv

		for ((m = $minkmer; m <= $MaxKmers; m += 1))
		do
			echo "edit $e, threshold $t, minKmers $m"
			rawfilename="L1s_raw_e${e}_t${t}_m${m}.csv"
			filteredfilename="L1s_filtered_e${e}_t${t}_m${m}.csv"
			time python3 "${DIR}/L1PD.py" $SAMfile $FASTAfile -t $t -m $m --data_dir $DATA_DIR --csvoutput > $pr_dir/$rawfilename
			time python3 "${DIR}/filter_possible_L1s_CSV.py" $pr_dir/$rawfilename $FASTAfile -t $t --data_dir $DATA_DIR > $pr_dir/$filteredfilename
			# Use wc t get line count, with sed's help to remove trailing file name
			# Line counts correspond to Positives and True Positives.
			PatCount=$(wc -l $pr_dir/$rawfilename | sed 's/ .*//')
			if [ $PatCount -gt 0 ] ; then
				# True positives, false positives, and false negatives
				TPCount=$(grep -v '^DUPLICATE' $pr_dir/$filteredfilename | wc -l | sed 's/ .*//')
				FPCount=$((PatCount - TPCount))
				FNCount=$((L1Count - TPCount))
				# We assume the file with full-length intact L1s has 'FLI-L1'
				# in its name (upper or lower case).
				FLIL1sFound=$(grep -Fci 'FLI-L1' $pr_dir/$filteredfilename)
				# Calculate precision, recall, and F1 Score.
				# Bash has no floating point arithmetic, so bc to the rescue.
				precision=$(echo "scale=5; $TPCount/$PatCount" | bc)
				recall=$(echo "scale=5; $TPCount/$L1Count" | bc)
				F1Score=$(echo "scale=5; 2*$precision*$recall/($precision+$recall)" | bc)
				cat > $pr_dir/results_e${e}_t${t}_m${m}.txt <<- EOF
					FLI-L1s:   $FLIL1sFound 
					L1Count:   $L1Count
					PatCount:  $PatCount
					TPCount:   $TPCount
					Precision: $precision
					Recall:    $recall
					F1 Score:  $F1Score
				EOF
				printf "$m\t$PatCount\t$TPCount\t$FPCount\t$FNCount\t$precision\t$recall\t$F1Score\n" >> $pr_dir/results_e${e}_t${t}.csv
			else
				echo "PatCount = 0" > $pr_dir/results_e${e}_t${t}_m${m}.txt
			fi
		done
	done
done
rm $SAMfile

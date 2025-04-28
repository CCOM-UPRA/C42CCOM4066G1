#!/bin/bash
#SBATCH --job-name=precision_recall
#SBATCH --error=precision_recall_error.txt
#SBATCH --output=precision_recall_output.out
#SBATCH --workdir=.
#SBATCH --mem 1
#SBATCH -A trn035
#SBATCH -N 1
#SBATCH --time=08:00:00
######################################################################
# Precision & Recall
#
# Authors: Juan O. Lopez    (juano.lopez@upr.edu)
#          Emanuel Martinez (emanuel.martinez8@upr.edu)
#          Javier Quinones  (javier.quinones3@upr.edu)
# License: Creative Commons Attribution-ShareAlike 4.0
# http://creativecommons.org/licenses/by-sa/4.0/legalcode
######################################################################

#https://stackoverflow.com/questions/56962129/how-to-get-original-location-of-script-used-for-slurm-job
if [ -n "$SLURM_JOB_ID" ] ; then
	# Need to find dir with script, since SLURM copies scripts to a different directory.
	#command_line=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
	#DIR=$(dirname "$(echo "$command_line" | awk '{print $1}')")
	# User MUST supply workdir for this to work.
	DIR=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/WorkDir=/{print $2}')
else
	# Other L1PD scripts are in the same directory as this script.
	DIR=$(dirname $0)
fi

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
	-m|--MinKmers   Minimum amount of kmers to be used at once.
						[Default: 2]
	-h|--help      Help (shows usage)

	-G|--graph     Generates a graph from the results of precision and recall.

	-x|--x_axis    Specifies the header name for the x axis of the graph.

	-y|--y_axis    Specifies the header name for the y axis of the graph.

	-z|--z_axis    Specifies the header name for the z axis of the graph.
	
	-d|--three_d   Graph will be three dimesnional.
						[Default: false]
	-s|--stackable Graphs will be displayed on the same plain.
						[Default: false] 
	-o|--ouput Name of the output graph
EOF
}

#############################
# Initialize default values #
#############################

# Because many files are generated, store them in a subdirectory.
printdir="result"
MinKmers=2
graph="false"
x_axis="Edit dist"
y_axis="Threshold"
z_axis="F1 score"
stackable="false"
three_d="false"
output="graph.png"
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
			IFS=',' read -r EditStart EditEnd EditStep <<< "$edit_range"
			if [[ -z $EditEnd || -z $EditStep ]]; then
				echo "Invalid edit range format: $edit_range"
				echo "Use a comma to separate start value, end value, and increment."
				echo "For example: 5,30,5"
				exit 1
			fi
			shift 2
			;;
		-t|--treshold)
			treshold_range="$2"
			IFS=',' read -r TresholdStart TresholdEnd TresholdStep <<< "$treshold_range"
			if [[ -z $TresholdEnd || -z $TresholdStep ]]; then
				echo "Invalid treshold range format: $treshold_range"
				echo "Use a comma to separate start value, end value, and increment."
				echo "For example: 100,800,50"
				exit 1
			fi
			shift 2
			;;
		-m|--MinKmers)
			MinKmers="$2"
			shift 2
			;;
		-G|--graph)
                        graph="true"
                        shift 1
                        ;;
		-x|--x_axis) # Since x can have multiple values we must iterate to check the next argument
			shift 1 # move past the -x flagi
			unset x_axis 
                        while [[ $# -gt 0 && "$1" != -* ]]; do # while value isnt another parameter grab the strings
        			x_axis+=("$1")
        			shift 1 # move on to next string
      			done 
                        ;;
		-y|--y_axis)
                        shift 1 # move past the -y flag
			unset y_axis
                        while [[ $# -gt 0 && "$1" != -* ]]; do # while value isnt another parameter grab the strings
        			y_axis+=("$1")
        			shift 1 # move on to next string
      			done 
                        ;;
		-z|--z_axis)
                        shift 1 # move past the -z flag
			unset z_axis
                        while [[ $# -gt 0 && "$1" != -* ]]; do # while value isnt another parameter grab the strings
        			z_axis+=("$1")
        			shift 1 # move on to next string
      			done
			;;
		-d|--three_d)
                        three_d="true"
                        shift 1
                        ;;
		-s|--stackable)
                        stackable="true"
                        shift 1
                        ;;
		-o|--output)
                        output="$2"
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

#############################################################################
# Initialize default values for EditStart and End and TresholdStart and End #
#############################################################################

EditStart="${EditStart:-5}"
EditEnd="${EditEnd:-20}"
TresholdStart="${TresholdStart:-25}"
TresholdEnd="${TresholdEnd:-500}"
EditStep="${EditStep:-5}"
TresholdStep="${TresholdStep:-25}"

##########################################
# Simplify access to parameter variables #
##########################################

GENOME=$genome
FASTAfile=$fasta
SAMfile=$sam
DATA_DIR=$meta
L1Count=$l1count

############################################
# Variable declaration for triple for loop #
############################################

# Extract kmer size from fasta file
k=$(wc -L $FASTAfile | sed 's/ .*//')
# Extract and store total number of probes
MaxKmers=$(wc -l $FASTAfile | sed 's/ .*//') 
MaxKmers=$((MaxKmers / 2))
# As best F1 scores are obtained using more than 50% of the probes, run tests using 50% to 100% of probes
MinKmers=$((MaxKmers / 2))
if (( $MinKmers < 2 )) ; then # Avoid using only 1 probe from 1 ORF in case there are only 2 probes
    MinKmers=2
fi
# Backup of original MaxKmers & MinKmers value to be reset with each edit distance
OriginalMaxKmers=$MaxKmers
OriginalMinKmers=$MinKmers
# Keep track of F1 score to only store values in best_F1_history.csv when an improvement is registered in a given kmer size
maxF1=0 

#########################################################################
# Establish print directory and create header for final results storage #
#########################################################################

pr_dir=$printdir
mkdir -p $pr_dir # -p ignores errors if the directory exists
printf "K-mer size\tEdit dist\tThreshold\tK-mers min\tPrecision\tRecall\tF1 Score\n" > $pr_dir/${k}mers_best_F1_history.csv # Header and file creation for best F1 values csv
printf "K-mer size\tEdit dist\tThreshold\tK-mers min\tPrecision\tRecall\tF1 Score\n" > $pr_dir/${k}mers_F1Scores.csv # Header and file creation for all F1 values csv
#################################################################################################
# Vary edit distance (e), threshold(t) and ammount of probes used (m) to see effect on F1 Score # 
#################################################################################################
cat $GENOME | head -n 10 >> debug_output_file.txt  # Check first few lines of the genome file
cat $FASTAfile | head -n 10 >> debug_output_file.txt # Check first few lines of the fasta file
# Vary the edit distance used by mrFAST
for ((e = $EditStart; e <= $EditEnd; e += $EditStep))
do
	echo $e >> debug_output_file.txt
	# If genome is not indexed, index genome
	if [ ! -f $GENOME.index ] ; then
		mrfast-master/mrfast --index $GENOME
	fi
	# Generate SAM file for current edit distance
	time mrfast-master/mrfast --search $GENOME --seq $FASTAfile -e $e -o $SAMfile

	echo "time mrfast-master/mrfast --search $GENOME --seq $FASTAfile -e $e -o $SAMfile" >> debug_output_file.txt
	# Reset variables with each change in edit distance
	MaxKmers=$OriginalMaxKmers # Reset MaxKmers variable with original value when changing edit distance
	MinKmers=$OriginalMinKmers # Reset MinKmers variable with original value when changing edit distance
	BestM=0 # Keeps track of m value with best F1 value in a given edit distance to reduce run time and redundance
	BestMF1=0 # Keeps track of F1 score for BestM value in a given edit distance

	# Vary the threshold to see how that affects P&R
	for ((t = $TresholdStart; t <= $TresholdEnd ; t += $TresholdStep))
	do
		# Create a CSV with tab-separated values
		printf "K-mers min\tPatterns\tTrue +\tFalse +\tFalse -\tPrecision\tRecall\tF1 Score\n" > $pr_dir/results_e${e}_t${t}.csv

		#######################################################################################################
		#                                       MinKmer Runtime Reducer                                       #
		# As a specific m value provides best F1 Score in the mayority of cases with any threshold in a       #
		# given edit distance, do a test run using whole range of m values (50% to 100% of total probes)      #
		# for each edit distance with first treshold value. After seeing which m value gave highest F1 value  #
		# for current e value, repeat first and next treshold values with one m value for each edit distance. #
		#######################################################################################################
		
		for ((m = $MinKmers; m <= $MaxKmers; m += 1))
		do
			echo "edit $e, threshold $t, MinKmers $m"
			rawfilename="L1s_raw_e${e}_t${t}_m${m}.csv"
			filteredfilename="L1s_filtered_e${e}_t${t}_m${m}.csv"
			# Generate a CSV with possible L1s 
			time python3 "${DIR}/L1PD_files/L1PD.py" $SAMfile $FASTAfile -t $t -m $m --data_dir $DATA_DIR --csvoutput > $pr_dir/$rawfilename
			# Compare possible L1s from L1PD against L1s from L1Base2 to filter L1s with no matches
			time python3 "${DIR}/precision_recall_files/filter_possible_L1s_CSV.py" $pr_dir/$rawfilename $FASTAfile -t $t --data_dir $DATA_DIR > $pr_dir/$filteredfilename
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
				printf "$m\t$PatCount\t$TPCount\t$FPCount\t$FNCount\t$precision\t$recall\t$F1Score\n" >> $pr_dir/results_e${e}_t${t}.csv # Save data for best m score for e and t value
				if (( $MinKmers == $MaxKmers )) ; then # If current run is part of real runs or only 2 probes are availible (no need to run test run)
					printf "$k\t$e\t$t\t$m\t$precision\t$recall\t$F1Score\n" >> $pr_dir/${k}mers_F1Scores.csv # Store all F1 scores with their respectively used parameters
					if (( $(echo "$F1Score > $maxF1" |bc -l) )) ; then
						printf "$k\t$e\t$t\t$m\t$precision\t$recall\t$F1Score\n" >> $pr_dir/${k}mers_best_F1_history.csv # Stores only F1 scores that show improvement in a given kmer value
						maxF1=$F1Score # Stores current best F1 score to be used in if statmement comparison
					fi
					break # As best m value has been selected, skip process of selecting best m value below
				fi
				if (( $(echo "$F1Score > $BestMF1" |bc -l) )) ; then # Check if current m value gave higher F1 score than previous m values in current edit distance during test runs
					BestM=$m # Stores best m value for a given edit distance
					BestMF1=$F1Score # Stores F1 Score for BestM to be compared in if statement
				fi
			else
				echo "PatCount = 0" > $pr_dir/results_e${e}_t${t}_m${m}.txt
			fi
		done
		if (( $MinKmers != $MaxKmers )) ; then # If current run is test run
			t=$(( $t - $TresholdStep )) # Make next run have same t value as current one when going through for loop to start real runs
			# Set all future runs in this edit distance to only run m with BestM instead of range of m values
			MinKmers=$BestM 
			MaxKmers=$BestM
		fi
	done
done
#the SAM FILE can be really big in size, so we delete them. To keep them just delete next line of code.
rm $SAMfile

# Check if we want a graph to be made
if [[ "$graph" == "true" ]]; then
	
	# Case for a 3d stackable graph
	if [[ "$three_d" == "true" && "$stackable" == "true" ]]; then
  		python3 read_csv.py -i "$pr_dir/${k}mers_F1Scores.csv" -x "${x_axis[@]}" -y "${y_axis[@]}" -z "${z_axis[@]}" -d -s -o "$output"
	
	# Case for a 3d non stackable graph
	elif [[ "$three_d" == "true" ]]; then
  		python3 read_csv.py -i "$pr_dir/${k}mers_F1Scores.csv" -x "${x_axis[@]}" -y "${y_axis[@]}" -z "${z_axis[@]}" -d -o "$output"
	
	# Case for a 2d stackable graph
	elif [[ "$stackable" == "true" ]]; then
  		python3 read_csv.py -i "$pr_dir/${k}mers_F1Scores.csv" -x "${x_axis[@]}" -y "${y_axis[@]}" -s -o "$output"
	
	# Case for a 2d non stackable graph
	else
  		python3 read_csv.py -i "$pr_dir/${k}mers_F1Scores.csv" -x "${x_axis[@]}" -y "${y_axis[@]}" -o "$output"

	fi
fi

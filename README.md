# Precision and Recall OpenMP project

>  In this project we will be parallelizing the Precision and recall bash script. This README will inform you of the basic setup to be able to run it on the server.

## Setup Process
1. Rewrite the precision_recall_ss.sh and precision_recall.sh scripts to use your path and your username.
2. Set up the conda virtual environment to download the necesarry dependencies
3. Run the make command on the mrfast-master files to compile the script
4. Run the precision_recall_ss.sh bash script to ensure the script is running properly

> Below you will find instructions on how to set up the project. 

## Rewrite code to work in your directory

> Ensure the paths in precision_recall_ss.sh are properly formatted to your current work directory. You can check your directory by running the `pwd` command.

> Do not change the path for the genome (-g) parameter in the precision_recall_ss.sh file as this is the path that will be used to avoid duplicating large datasets. 

> Ensure the SLURM portions at the top of both the precision_recall_ss.sh and precision_recall.sh script have your email and username for execution purposes.

> Unzip the FLI_ORF2_FLnI-L1.fasta_FLnI-L1.fasta.zip file inside the L1Base2 directory. This data was compressed due to github file size limitations. 

> To unzip the file run the `unzip` command along with the filename
`unzip FLI_ORF2_FLnI-L1.fasta_FLnI-L1.fasta.zip`

## Use the Conda virtual environment to set up the dependencies needed for our scripts
### To install conda run
`wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh`

`bash Miniforge3-Linux-x86_64.sh`

### If miniforge is installed but conda isn't working run the command to initialize it
`source ~/miniforge3/bin/activate`

`conda init bash`


### To create the virtual environment and install the biopython version used in the project run the following command

`conda create -n pr_env python=3.6.10 biopython=1.79 `

### Once created you can activate the environment by running
`conda activate pr_env`

### Create the orf12_sam directory
`mkdir orf12_sam`

## Run the make command to initialize mrFAST

`cd mrfast-master `

`make`

## Run the precision_recall_ss.sh script using SLURM
` sbatch precision_recall_ss.sh `

> Keep in mind this will take some time to complete be patient and check the precision_recall_ss_error.txt for any error messages to trouble shoot


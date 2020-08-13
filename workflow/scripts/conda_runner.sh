#!/bin/bash 

#--------------------------------------------------#
#
#	BASH script for running a snakemake pipeline
#                    
#-------------------------------------------------#

#==============================================
# set diretory
#==============================================


dir=$( pwd )

cur1=$( echo ${dir} | rev | cut -d'/' -f 1 | rev )
cur2=$( echo ${dir} | rev | cut -d'/' -f 2 | rev )
cur3=$(echo ${dir} | rev | cut -d'/' -f 3 | rev )

if  [ "${cur1}" == "snakeAmpSeq" ]
then
	echo "\

        Starting SnakeMake

"
elif [ "${cur2}" == "snakeAmpSeq" ]
then
	cd ../
	echo "\

        Starting SnakeMake

"
elif [ "${cur3}" == "snakeAmpSeq" ]
then
	cd ../../
	echo "\

        Starting SnakeMake

"

else
	echo "\

	this script must be executed from the snakeAmpSeq dir or the
	snakeAmpSeq/workflow/scripts dir. Please make sure you are
	in one of these folders to execute this bash helper script.
	You are currently in: ${dir}

" && exit

fi


#============================================
# set env here we find local install of conda
# in the path and then setup the script to
# find users version of conda
# ${CONDA_PREFIX} is a variable in the BASH env
# If you have conda installed and run into trouble
# check to see that this variable is set in
# your BASH env `env | grep CONDA_PREFIX
# if it is not set, you may need to set it
# on the fly in your shell or in your .bashrc
#============================================

if [ $( echo ${CONDA_PREFIX} | wc -c ) -lt 3 ];
then
	echo "\
	neither conda or nor anaconda found in path
	"&& exit
fi



PATH=${PATH}:${CONDA_PREFIX}/bin
set -x
set -e 
# activate the conda env
source ${CONDA_PREFIX}/etc/profile.d/conda.sh

# determine if the snakemake-mapping has already been installed

if [ ! -d ${CONDA_PREFIX}/envs/snakemake-mapping/ ]
then
    # if not installed go ahead and create it now using the following command
    # that implements environment file which it finds via the relative path.
    # the name is in the snakeAmpSeq_environment.yml
    conda env create --file workflow/envs/snakeAmpSeq_environment.yml
fi
# now we can activate  the file
conda activate ${CONDA_PREFIX}/envs/snakemake-mapping/


#============================================
# Now we can run the snakemake
# pipeline. It takes stems from
#  resources/downloads folder
# and executes snakemake using the stems.
# note that the stems assume that your
# fastq files are still compressed.
# additionally, and mentioned in the
# quick start, they should end in
#`left=<uniq_id>_R1_001.fastq.gz`, and
#`<right=<uniq_id>_R2_002.fastq.gz>`
#============================================


Y=`ls resources/downloads/*gz | rev| cut -d'/' -f -1 | rev | cut -d'_' -f -1 | sort | uniq`

for X in ${Y}
do 
 
	snakemake -p results/4.tsr_report/${X}_tsr_report.tsv
	snakemake -p results/7.tally/${X}.tsv
done


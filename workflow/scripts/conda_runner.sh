#!/bin/bash 

#---------------------------------------------
#
#	BASH script for running snakemake
#	Pipeline (Sure this could eventually
#       be made into a pure python/conda 
#       solution)
#                    
#----------------------------------------------

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


#==============================================
# set env here we find local install of conda
# in the path and then setup the script to
# find users version of conda
#==============================================

conda_path=$( which anaconda  )

if [ $( echo ${conda_path} |wc -c ) -lt 6 ]
then
conda_path=$( which conda )
fi

conda_count=$( echo ${conda_path} | wc -c )
conda_dir=$( echo ${conda_path} | rev |cut -d'/' -f 3- |rev )

if [ ${conda_count} -lt 6 ];
then
	echo "\
	This program requires neither conda or nor anaconda found in pathwgigit
	"&& exit
fi



PATH=${PATH}:${conda_dir}/bin
set -x
set -e 
# activate the conda env
source ${conda_dir}/etc/profile.d/conda.sh
conda activate ${conda_dir}/envs/snakemake-mapping/

# get the stems we need from refs:

Y=`ls resources/downloads/*gz | rev| cut -d'/' -f -1 | rev | cut -d'_' -f -1 | sort | uniq`

for X in ${Y}
do 
 
	snakemake -p results/4.tsr_report/${X}_tsr_report.tsv
	snakemake -p results/7.tally/${X}.tsv
done


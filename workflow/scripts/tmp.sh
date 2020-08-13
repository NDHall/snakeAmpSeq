#!/bin/bash 

dir=$( pwd )
error_msg="\

this script must be executed from the snakeAmpSeq dir or the\n
snakeAmpSeq/workflow/scripts dir. Please make sure you are\n 
in one of these folders to execute this bash helper script.\n
You are currently in: ${dir}\n\n"






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

conda=$( which anaconda |wc -c )

if [ ${conda} -lt 6 ];
then
	echo "\
	This program requires anaconda to be installed
	"&& exit
fi 


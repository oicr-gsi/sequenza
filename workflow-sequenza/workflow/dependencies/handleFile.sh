#!/bin/bash
wdPath=`pwd`
echo $wdPath
extid=$1
inPath=$wdPath/$2
cd $inPath
model=model-fit; mkdir $model

for f in $( ls ); do
	if [[ -d $f ]]; then
		if [[ $f != "model-fit" ]]; then
			cd $inPath
			mv $f $model
		fi
	else
		if [[ ${f: -4} != ".seg" ]] && [[ ${f: -15} != "confints_CP.txt" ]]; then	
                         mv $f $model
		fi 
	fi
done

cd $model
cd $inPath
tar -zcvf ${model}.tar.gz $model
#!/bin/bash
wdPath=`pwd`
echo $wdPath
extid=$1
inPath=$wdPath/$2
cd $inPath
model=$inPath/model-fit; mkdir $model
for f in $( ls ); do
	if [[ -d $f ]]; then
		if [[ $f != "model-fit" ]]; then
			cd $inPath
			mv $f $model
		fi
	else
		if [[ ${f: -4} != ".seg" ]]; then
			mv $f $model
		fi
	fi
done

cd $model
cp ${extid}_genome_vie*.pdf $inPath
cd $inPath
tar -zcvf ${model}.tar.gz $model
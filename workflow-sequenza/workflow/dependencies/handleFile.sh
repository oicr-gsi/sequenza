#!/bin/bash
extid=$1
inPath=$2
cd $inPath
model=$inPath/model-fit; mkdir $model
for f in $( ls ); do
	if [[ -d $f ]]; then
		if [[ $f != "model-fit" ]]; then
			cd $f;
			cp ${extid}_genome_view.pdf ../${extid}_genome_view_${f}.pdf;
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
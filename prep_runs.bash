#!/bin/sh

rm -fr $1_sd* $1_dd*
cat $1* > EX_$1

maxl=200
if [ `grep -c rd0 EX_$1` -gt 0 ]
then
	echo "Treating single diffusion cases..."
	grep -A 1 rd0 EX_$1 | sed '/^--$/d' | sed -r 's/.*(.\/execute).*/\1/g' > $1_sd.tmp

	mkdir $1_sd
	nd=$(((`grep -c "" $1_sd.tmp` / $maxl + 1) / 10 + 1))
	split -dl $maxl -a $nd $1_sd.tmp $1_sd/
	head -3 `head -1 $1_sd.tmp | sed -r 's/[a-z ]*(\/.*)/\1/g'`/execute > headerfile
	cmd=`head -2 EX_$1 | tail -1 | sed -r 's/.\/execute//g'`
	i=0
	for file in $1_sd/*
	do
		cat headerfile $file > tmpfile
		mv tmpfile $file

		echo $cmd$file >> $1_sd_$(($i/500))
		echo "sleep 0.05" >> $1_sd_$(($i/500))
		let i=$i+1
	done
	rm headerfile
	rm $1_sd.tmp
	chmod u+x $1_sd/*
fi

maxl=900
if [ `grep rd[^0] -c EX_$1` -gt 0 ]
then
	echo "Treating double diffusion cases..."
	grep -A 2 rd[^0] EX_$1 | sed '/^--$/d' > $1_dd.tmp
	
	nd=$(((`grep -c "" $1_dd.tmp` / $maxl + 1) / 10 + 1))
	split -dl $maxl -a $nd $1_dd.tmp $1_dd_
	rm $1_dd.tmp
fi


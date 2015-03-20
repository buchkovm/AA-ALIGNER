BASE=$2
GROOMERLIBRARY=$4
GROOMERPATH=$5
BASE=`echo $BASE |sed 's/\.gz//' |sed 's/\.fastq//'`

FORMAT=`unzip -p $3/${BASE}_fastqc.zip "${BASE}_fastqc/fastqc_data.txt" |grep Encoding |cut -f2 |cut -f1 -d' '`
if [ "$FORMAT" = "Sanger" ]; then
	`echo "Base qualities already in $FORMAT format" >$3/${BASE}.groomer_summary.txt`
	cp $1/$2 $3/$2
else
	#BASE=$2
	`echo "Base qualities in $FORMAT format. Converting to Sanger." >$3/${BASE}.groomer_summary.txt`
	FORMAT=`echo $FORMAT |tr '[:upper:]' '[:lower:]'`
	echo $FORMAT
	export PYTHONPATH=$PYTHONPATH:$GROOMERLIBRARY
	if [[ "$2" =~ ".gz" ]];then
		echo "UNZIP"
		`gunzip -c $1/$2 >$3/${BASE}.${FORMAT}.fastq`
		#mv $1/$2 $1/${BASE}.${FORMAT}.fastq.gz
		echo "UNZIP COMPLETE"
		#BASE=`echo $BASE |sed "s/.gz//"`
		`python fastq_groomer.py $3/${BASE}.${FORMAT}.fastq $FORMAT $3/${BASE}.fastq sanger ascii summarize_input >$3/${BASE}.groomer_summary.txt`
		echo "GROOMING COMPLETE"
		gzip $3/${BASE}.fastq
	else
		echo "NOT ZIPPED" 
		cp $1/${BASE}.fastq $3/${BASE}.${FORMAT}.fastq
		`python /nas02/apps/galaxy-prod/galaxy-dist/tools/fastq/fastq_groomer.py $3/${BASE}.${FORMAT}.fastq $FORMAT $3/${BASE}.fastq sanger ascii summarize_input >$3/${BASE}.groomer_summary.txt`	
		gzip $3/${BASE}.fastq
		#gzip $1/${BASE}.${FORMAT}.fastq
	fi
	#BASE=`echo $BASE |sed "s/.fastq//"`
	#echo $BASE
	#gzip $1/${BASE}.fastq
fi;	

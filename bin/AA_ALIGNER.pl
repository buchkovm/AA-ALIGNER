use Getopt::Long;
use Cwd;

my %opts=(all=>1, temp_directory=>"");
my %config_opts=();
GetOptions(\%opts, qw(all create_reference filter_fastq align filter_alignment find_imbalance out_directory:s temp_directory:s out_prefix:s configuration_file:s));


my @files = @ARGV;
if ($opts{filter_fastq} || $opts{align} || $opts{filter_alignment} || $opts{allelic_imbalance} || $opts{create_reference} ||$opts{find_imbalance}){
	$opts{all} = 0;
}
$opts{temp_directory} = $opts{out_directory} if ($opts{temp_directory} eq "");

open(IN, $opts{configuration_file});
while(<IN>){
	chomp;
	if(substr($_,0,1) ne "#"){
		my @cols = split /\t/, $_;
		my @parts = split /=/, $cols[0];
		$config_opts{$parts[0]}=$parts[1];
	}
}
close(IN);

my $JAVAPATH = "java";
$JAVAPATH = $config_opts{JAVA_PATH} . "/" . "java" if ($config_opts{JAVA_PATH} ne "");

my $PERLPATH = "perl";
$PERLPATH = $config_opts{JAVA_PATH} . "/" . "perl" if ($config_opts{JAVA_PATH} ne "");

my $PYTHONPATH = "python";
$PYTHONPATH = $config_opts{JAVA_PATH} . "/" . "python" if ($config_opts{JAVA_PATH} ne "");

my $read_length = $config_opts{READ_LENGTH};
my @alignment_sams=();
my $cur_dir = getcwd;

filterFastq() if ($opts{filter_fastq} || $opts{all});
determineKmer();
createReference() if ($opts{create_reference});
align() if ($opts{align} || $opts{all});
filterAlignment() if ($opts{filter_alignment} || $opts{all});
findImbalance() if ($opts{find_imbalance} || $opts{all});

##Pre-alignment fastqc filtering
sub filterFastq{
	foreach my $file (@files){
		my @parts = split /\//, $file;
		my $dir ="";
		for (my $i=0; $i<$#parts; $i++){
			$dir = $dir . "$parts[$i]\/";
		}	
		$prefix = $parts[$#parts];
		$prefix =~s/.fastq//;
		$prefix =~s/.gz//;
		print STDERR "\nPre-alignment fastq filtering\n\n";
		print "\nPre-alignment fastq filtering\n\n";
		##Generate fastqc report
		print STDERR "Generating Fastqc report ....";
		#my $err = `fastqc --noextract --outdir $opts{out_directory} -t $config_opts{FASTQC_THREADS}  $file`;
		#print STDERR "$err";
		my $fastqc_cmd = "fastqc --noextract --outdir $opts{out_directory} -t $config_opts{FASTQC_THREADS}  $config_opts{FASTQC_OPTIONS} $file";
		`$fastqc_cmd`;
		print "$fastqc_cmd\n";
		##Use fastq report to check base quality format and convert to Sanger/Illumina 1.9 if necessary
		print STDERR "done\n";
		print STDERR "Calculating read length ....";
		my $read_length_cmd = "unzip -p $opts{out_directory}/$prefix" . "_fastqc.zip \"$prefix"."_fastqc/fastqc_data.txt\" |head -n 9 |tail -n 1|cut -f 2";
		my $command = "unzip -p $opts{out_directory}/$prefix" . "_fastqc.zip \"$prefix"."_fastqc/fastqc_data.txt\" |head -n 9 |tail -n 1|cut -f 2";
		my $length = `$command`;
		print "$command\n";
		chomp($length);
		my @length_parts = split /-/, $length;
		$read_length = $length_parts[0];
		print STDERR "$length\n";
	
		print STDERR "Checking base quality score format ....";
		checkAndConvertQualityScores($file);
		##Unzip to temporary file if necessary
		my $base=$parts[$#parts];
		$base =~ s/.gz//;
		my $unzipped_file = $base;
		if (-e "$dir"."$base.gz"){
			print STDERR "Unzipping file ....";
			my $gunzip_cmd ="gunzip -c $dir/$base.gz > $opts{temp_directory}/$base";
			`$gunzip_cmd`;
			print "$gunzip_cmd\n";
			$unzipped_file = "$opts{temp_directory}/$base";
			print STDERR "done\n";
		}
		print STDERR "Trimming fastq bases ....";
		##FASTX Trimmer to trim bases to length in configuration file
		if($config_opts{TRIM_LAST_BASE} eq "SEQUENCE_LENGTH"){
			$config_opts{TRIM_LAST_BASE}=$read_length;
		}
		my $fastx_trim_cmd = "fastx_trimmer -f $config_opts{TRIM_FIRST_BASE} -l $config_opts{TRIM_LAST_BASE} -Q 33 -o $opts{temp_directory}/trimmed_$prefix.fastq $config_opts{FASTX_OPTIONS} -i $unzipped_file";
		$err = `$fastx_trim_cmd`;
		print STDERR "$err";
		print "$fastx_trim_cmd\n";
		print STDERR "done\n";
		
		##FASTX Trimmer to trim put low quality sequences
		print STDERR "Filtering out low quality bases ....";
		my $fast_qual_cmd = "fastq_quality_filter -Q 33 -p $config_opts{FASTQ_QUALITY_FILTER_P} -q $config_opts{FASTQ_QUALITY_FILTER_Q} $config_opts{FASTQ_QUALITY_OPTIONS} -i $opts{temp_directory}/trimmed_$prefix.fastq -o $opts{out_directory}/$base.filtered";
		$err = `$fast_qual_cmd`;
		print "$fast_qual_cmd\n";
		print STDERR "$err";
		print STDERR "done\n";
		print STDERR "Compressing file ....";
		my $gzip ="gzip $opts{out_directory}/$base.filtered";
		`$gzip`;
		print "$gzip\n";
		print STDERR "done\n";
	}
}
sub checkAndConvertQualityScores{
	my $file = shift;
	my $err = 0;	
	
	my @parts = split /\//, $file;
	my $dir ="";
	for (my $i=0; $i<$#parts; $i++){
		$dir = $dir . "$parts[$i]\/";
	}
	my $base = $parts[$#parts];
	$base=~s/\.gz//;
	$base=~s/\.fastq//;
	
	##Extract format from fastqc output
	my $check_format_command = "unzip -p $opts{out_directory}/$base"."_fastqc.zip \"$base"."_fastqc/fastqc_data.txt\" |grep Encoding |cut -f2 |cut -f1 -d' '";
	my $format = `$check_format_command`;
	print "$check_format_command\n";
	chomp($format);
	if ($format eq "Sanger"){
		print STDERR "Format is already in Sanger / Illumina 1.9\n";
	}else{
		print STDERR "Format is $format\n";
		print STDERR "Converting to Sanger/Illumina 1.9 format ....";
		$format=lc($format);
		if (-e "$dir"."$base.gz"){
			##Unzips fastq with wrong format to renamed copy
			my $gunzip_cmd = "gunzip -c $file > $opts{temp_directory}/$base.$format.fastq";
			`$gunzip_cmd`;
			print "$gunzip_cmd\n";
			##Rename original zipped fastq to contain format information 
			my $mv_cmd = "mv $file $dir$base.$format.fastq.gz";
			`$mv_cmd`;
			print "$mv_cmd\n";
			##Add python library
			`export PYTHONPATH:/$PYTHONPATH:$config_opts{GALAXY_LIBRARY_PATH}`;
			##Convert base quality format
			my $groomer_cmd ="$PYTHONPATH $config_opts{GROOMER_PATH}/fastq_groomer.py $opts{temp_directory}/$base.$format.fastq $format $dir/$base.fastq sanger ascii summarize_input >$opts{out_directory}/$base.groomer_summary.txt";
			my $err = `$groomer_cmd`;
			print STDERR "$err\n";
			print "$groomer_cmd\n";
			my $gzip_cmd = "gzip $dir/$base";
			`$gzip_cmd`;
			print "$gzip_cmd\n";
		}elsif(-e $base){
			##Rename original fastq file to contain format information
			my $rename_cmd = "mv $dir$base.fastq $dir$base.$format.fastq";
			`$rename_cmd`;
			print "$rename_cmd\n";
			my $groomer = "$PYTHONPATH $config_opts{GROOMER_PATH}/fastq_groomer.py $dir/$base.$format.fastq $format $dir/$base.fastq sanger ascii summarize_input >$opts{out_directory}/$base.groomer_summary.txt";
			my $err = `$groomer_cmd`;
			print STDERR "$err\n";
			print  "$groomer_cmd\n";
			my $gzip1 = "gzip $dir/$base.fastq";
			`$gzip1_cmd`;
			print "$gzip1_cmd\n";
			my $gzip2 = "gzip $dir/$base.$format.fastq"; 
			`$gzip2_cmd`;
			print "$gzip2_cmd\n";
		}else{
			print STDERR "Having trouble finding $file\n";
			$err = 1;
		}	
		if($err == 0){
			print STDERR "done\n";
		}else{
			print STDERR "error\n";
		}
	}
	
		
}
sub determineKmer{
	print STDERR "Calculating k-mer ....";
	my $mm = $config_opts{MISMATCHES};
	my $k = $read_length/($mm+1);
	if ($k<10){
		print STDERR "$k:K-mer < 10 is not reocommended. Setting k-mer to 10, but this can cause inaccurate sequeunce alignment. Try a lower mismatch threshold or use set FORCE_KMER=1\n";
		$k = 10 if (! $config_opts{FORCE_KMER});
	}elsif ($k>=10 && $k<12){
		$k = 10;
		print STDERR "10\n";
	}elsif ($k>=12 && $k<15){
		$k = 12;
		print STDERR "12\n";
	}elsif ($k>=15){
		$k = 15;
		print STDERR "15\n";	
	}
	$config_opts{KMER}=$k;
	$config_opts{BASESIZE}=$k if ($config_opts{BASESIZE} eq "KMER");	
}
sub createReference{
	print STDERR "\nCreating personalized reference\n\n";
	print  "\nCreating personalized reference\n\n";
	
	my $genome_dir =  "$config_opts{GENOME_DIR}\/$config_opts{GENOME_NAME}";
	mkdir $genome_dir;
	
	createReferenceFa($genome_dir);
	buildGenome($genome_dir);
	createIITFile($genome_dir);
	indexSNPs($genome_dir);

}

sub createReferenceFa{
	my $genome_dir = shift;
	if($config_opts{SPLIT_BY_CHROM}){
	
		#Get prefix for split chr filename
		my @hets_path = split /\//, $config_opts{HET_FILE};
		if ($#hets_path >= 0){
			$hets_path[$#hets_path] =~ s/\.bed$//;
			splitFilesByChrom($config_opts{HET_FILE}, $hets_path[$#hets_path]);
		}
		my @replace_path = split /\//, $config_opts{REPLACE_FILE};
		if ($#replace_path >= 0){
			$replace_path[$#replace_path] =~ s/\.bed$//;
			splitFilesByChrom($config_opts{REPLACE_FILE}, $replace_path[$#het_path]);
		}

		#Create chromosome array
		my @chrom = `ls $config_opts{REF_FASTQ_DIR}/chr*.fa*`;	
		foreach my $file (@chrom){
			$file =~ m/(chr.*).fa/;
			my $chr = $1;
			my $createChrFa = "";
			if($config_opts{REFERENCE_FILE_COMPRESSED}){
				$createChrFa="gunzip -c $config_opts{REF_FASTQ_DIR}/chr$chr.fa.gz |$PERLPATH $config_opts{PATH_TO_AAALIGNER_BIN}/replaceReferenceAlleles.pl -output $genome_dir/chr$chr.fa -hets $opts{temp_directory}/$hets_path[$#hets_path].chr$chr.bed -replace $opts{temp_directory}/$replace_path[$#replace_path].chr$chr.bed >>$genome_dir/$config_opts{GENOME_NAME}.replace.txt 2>>$genome_dir/$config_opts{GENOME_NAME}.replaced.log"; 
			}else{
				$createChrFa="cat $config_opts{REF_FASTQ_DIR}/chr$chr.fa |$PERLPATH $config_opts{PATH_TO_AAALIGNER_BIN}/replaceReferenceAlleles.pl -output $genome_dir/chr$chr.fa -hets $opts{temp_directory}/$het_path[$#het_path].chr$chr.bed -replace $opts{temp_directory}/$replace_path[$#replace_path].chr$chr.bed  >>$genome_dir/$config_opts{GENOME_NAME}.replace.txt 2>>$genome_dir/$config_opts{GENOME_NAME}.replaced.log "; 
			}
			print STDERR "Creating chr$chr.fa ...";
			my $err = `$createChrFa`;
			print STDERR "$err\n";
			print "$createChrFa\n";	
			print STDERR "done\n";
		}
				
	}else{
		my @fastq_file;
		my $createChrFa = "";
		
		if($config_opts{REFERENCE_FILE_COMPRESSED}){
				$createFa="gunzip -c $config_opts{REF_FASTQ_DIR}/$config_opts{REF_FASTQ} |$PERLPATH $config_opts{PATH_TO_AAALIGNER_BIN}/replaceReferenceAlleles.pl -output $genome_dir/$config_opts{GENOME_NAME} -hets $config_opts{HET_FILE} -replace $config_opts{REPLACE_FILE}  >>$genome_dir/$config_opts{GENOME_NAME}.replace.txt 2>>$genome_dir/$config_opts{GENOME_NAME}.replaced.log"; 	
		}else{
				$createFa="cat $config_opts{REF_FASTQ_DIR}/$config_opts{REF_FASTQ} |$PERLPATH $config_opts{PATH_TO_AAALIGNER_BIN}/replaceReferenceAlleles.pl -output $genome_dir/$config_opts{GENOME_NAME} -hets $config_opts{HET_FILE} -replace $config_opts{REPLACE_FILE}  >>$genome_dir/$config_opts{GENOME_NAME}.replace.txt 2>>$genome_dir/$config_opts{GENOME_NAME}.replaced.log"; 	
		}
		print STDERR "Creating reference genome ...";
		my $err = `$createFa`;
		print STDERR "$err\n";
		print "$createFa\n";	
		print STDERR "done\n";
				
	}	
}

sub splitFilesByChrom{
	my $file =shift;
	my $prefix = shift;
	print STDERR "Splitting $file to $prefix.<CHROM>.bed ...";
	my $split_cmd = "$PERLPATH $config_opts{PATH_TO_AAALIGNER_BIN}/splitBedByChrom.pl $file $opts{temp_directory}/$prefix"; 	
	my $err = `$split_cmd`;
	print STDERR "$err\n";
	print "$split_cmd\n";
	print STDERR "done\n";
}

sub createGSNAPFile{
	my @hets_path = split /\//, $config_opts{HET_FILE};
	$hets_path[$#hets_path] =~ s/\.bed//;
	print STDERR "Creating gsnap file ...";	
	my $create_file = "awk '\{OFS=\"\"\; print \">\"\$4\" \"\$1\":\"\$3\" \"\$5\$6\}' $config_opts{HET_FILE} \>$opts{temp_directory}\/$hets_path[$#hets_path].gsnap.txt";
	my $err = `$create_file`;
	print STDERR "$err\n";
	print "$create_file\n";
	print STDERR "done\n";
}

sub buildGenome{
	my $genome_dir = shift;
	print STDERR "Building genome ...";
	my $dir =  "$config_opts{GENOME_DIR}/$config_opts{GENOME_NAME}";
	chdir $dir;
	my $buildGenome = "gmap_build -D $config_opts{GENOME_DIR} -d $config_opts{GENOME_NAME} -k $config_opts{KMER} -b $config_opts{BASESIZE}  -q $config_opts{SAMPLING} $config_opts{GMAP_BUILD_OPTIONS} $genome_dir/*.fa";
	my $err = `$buildGenome`;
	chdir $cur_dir;
	print STDERR "$err\n";
	print "$buildGenome\n";  
	print "done\n"
} 

sub createIITFile{
	my $genome_dir = shift;
	createGSNAPFile();
	print STDERR "Converting gsnap to iit_file\n";
	my @hets_path = split /\//, $config_opts{HET_FILE};
	$hets_path[$#hets_path] =~ s/\.bed//;
	my $iit_cmd = "cat $opts{temp_directory}/$hets_path[$#hets_path].gsnap.txt |iit_store $config_opts{IIT_STORE_OPTS} -o  $genome_dir/$config_opts{SNP_DATABASE_NAME}.iit";
	my $err =`$iit_cmd`;
	print STDERR "$err\n";
	print "$iit_cmd\n";
	print "done\n";	 
}

sub indexSNPs{
	my $genome_dir = shift;
	print STDERR "Indexing SNPs ...";
	my $indexSNPs = "snpindex  -D $config_opts{GENOME_DIR} -d $config_opts{GENOME_NAME} -k $config_opts{KMER} -b $config_opts{BASESIZE}  -q $config_opts{SAMPLING} -v $config_opts{SNP_DATABASE_NAME} $config_opts{SNPINDEX_OPTIONS} $genome_dir/$config_opts{SNP_DATABASE_NAME}.iit";
	my $dir =  "$config_opts{GENOME_DIR}/$config_opts{GENOME_NAME}";
	chdir $dir;
	my $err = `$indexSNPs`;
	chdir $cur_dir;
	print STDERR "$err\n";
	print "$indexSNPs\n";
	print "done\n";	 
}

sub align{
	print STDERR "\nAligning sequences\n\n";
	print  "\nAligning sequences\n\n";
	foreach my $file (@files){
		my @parts = split /\//, $file;
		my $base = $parts[$#parts];
		my $filename = "";
		if($opts{filter_fastq} || $opts{all}){
			$filename = "$opts{out_directory}/$base.filtered.gz";
			$base = "$base.filtered.gz";
		}else{
			$filename = "$file";
		}
		my @parts = split /\./, $filename;
		my $gunzip = "";
		if($parts[$#parts] eq "gz"){
			$gunzip = "--gunzip";
		}
		if($config_opts{INDEL_PENALTY} eq "NOINDELS"){
			$config_opts{INDEL_PENALTY}= $config_opts{MISMATCHES} + 1;
		}
		print STDERR "Aligning reads....";
		print "$align\n";
		my $align = "gsnap --terminal-threshold=10  -A sam --read-group-library $config_opts{READ_GROUP_LIBRARY} --read-group-platform $config_opts{READ_GROUP_PLATFORM} --read-group-id $config_opts{READ_GROUP_ID} -d $config_opts{GENOME_NAME} -D $config_opts{GENOME_DIR} -v $config_opts{SNP_DATABASE_NAME} -k $config_opts{KMER} --basesize $config_opts{BASESIZE} --sampling $config_opts{SAMPLING} -m $config_opts{MISMATCHES} -n $config_opts{NUMBER_READS_OUTPUT} -i $config_opts{INDEL_PENALTY}  --query-unk-mismatch=$config_opts{QUERY_UNK_MISMATCH} --genome-unk-mismatch $config_opts{GENOME_UNK_MISMATCH}  -t $config_opts{GSNAP_THREADS}  --trim-mismatch-score=$config_opts{TRIM_MISMATCH_SCORE} $gunzip $config_opts{GSNAP_OPTIONS} $filename  >$opts{out_directory}/$base.sam"; 
		my  $err = `$align`;
		print STDERR "$err\n";
		print STDERR "done\n";
		push (@alignment_sams, "$opts{out_directory}/$base.sam");
	}
}
	
sub filterAlignment{
	print STDERR "\nPost-alignment filtering\n\n";
	print  "\nPost-alignment filtering\n\n";
	my @samfiles = @files;
	if($opts{align} || $opts{all}){
		@samfiles = @alignment_sams;	
	}
	my @filtered_files;
	
	foreach my $sam (@samfiles){
		my @parts = split /\//, $sam;
		my $base = $parts[$#parts];
		$base =~ s/\.sam//;
		
		##Separates alignment into uniquely, multiple and not mapped reads
		print STDERR "Filtering to keep unique reads...."; 
		my $filter_unique = "$PERLPATH $config_opts{PATH_TO_AAALIGNER_BIN}/printUniqueSam.pl --header -format gsnap -output $opts{temp_directory}/$base $sam";
		print "$filter_unique\n";
		my $err = `$filter_unique`;
		print STDERR "$err done\n";
		
		##Filters out any alignments overalapping regions in specified blacklist (bed format) file
		print STDERR "Filtering blacklist regions....";
		my $blacklist = "$PERLPATH $config_opts{PATH_TO_AAALIGNER_BIN}/filterBlacklistSam.pl -output $opts{temp_directory}/$base.unique.blacklist -filter $config_opts{BLACKLIST_FILE} -input $opts{temp_directory}/$base.unique.sam";
		print "$blacklist\n";
		$err = `$blacklist`;
		print STDERR "$err done\n";
		

		##Converts to bam
		print STDERR "Converting to sorted bam....";
		my $toSortedBam = "samtools view -Sb -h $opts{temp_directory}/$base.unique.blacklist.sam |samtools sort - $opts{temp_directory}/$base.unique.blacklist";	
		print "$toSortedBam\n";
		$err = `$toSortedBam`;
		print STDERR "$err done\n";
	
		##Remove duplicate reads 
		print STDERR "Removing duplicates....";
		my $removeDups = "$JAVAPATH -Xmx4g -jar $config_opts{PATH_TO_PICARD}/MarkDuplicates.jar REMOVE_DUPLICATES=TRUE METRICS_FILE=$opts{out_directory}/$base.unique.blacklist.markDuplicates.metric.txt OUTPUT=$opts{temp_directory}/$base.unique.blacklist.markDuplicates.bam VALIDATION_STRINGENCY=LENIENT INPUT=$opts{temp_directory}/$base.unique.blacklist.bam $config_opts{MARK_DUPLICATES_OPTIONS}";
		print "$removeDups\n";	
		$err = `$removeDups`;
		print STDERR "$err done\n";
		push (@filtered_files,"$opts{temp_directory}/$base.unique.blacklist.markDuplicates.bam");

	}
	##Merge replicate files together
	print STDERR "Merging replicate bam files....";
	my $merge_str = "$JAVAPATH -Xxm4g -jar $config_opts{PATH_TO_PICARD}/MergeSamFiles.jar OUTPUT=$opts{out_directory}/$opts{out_prefix}.bam USE_THREADING=$config_opts{PICARD_USE_THREADING} $config_opts{MERGE_BAM_OPTIONS}";

	foreach my $rep (@filtered_files){
		$merge_str .= " INPUT=$rep";
	}
	print "$merge_str\n";
	$err = `$merge_str`;
	print STDERR "$err done\n";
}
		
sub findImbalance{
	##Initial round of imbalance detection
	print STDERR "\nDetecting Imbalances\n\n";
	print  "\nDetecting Imbalances\n\n";
	if($config_opts{INDEX_BAM}){
		print STDERR "Indexing bam file....";
		my $index = "samtools index $opts{out_directory}/$opts{out_prefix}.bam";
		print "$index\n";
		my $err = `$index`;
		print STDERR "$err done\n";
	}
	
	#Create output directories
	print STDERR "Making directories....";
	`mkdir $opts{out_directory}/chrom_alignments`;
	`mkdir $opts{out_directory}/chrom_alignments/$config_opts{IMBALANCE_DIR_NAME}`;
	
	my $forceMismatches = "";
	$forceMismatches = "-forceMismatches" if ($config_opts{HET_MM_ONLY} == 1);	
	
	##Gets chromosome information from bam header
	my @header = `samtools view -H $opts{out_directory}/$opts{out_prefix}.bam`;
	my @chrom;	
	foreach my $line (@header){
		if($line=~m/@SQ\tSN:(.*)\t/){
			push(@chrom, $1);
		}	
	}
	if($config_opts{INDEL_PENALTY} eq "NOINDELS"){
		$config_opts{INDEL_PENALTY}= $config_opts{MISMATCHES} + 1;
	}
	my $completeImbalances = "";
	$completeImbalances = "-completeImbalance" if ($config_opts{COMPLETE_IMBALANCE} == 1);	
	my @chromosome_files = ();	
	foreach my $chr (@chrom){
		if($config_opts{INDEX_BAM}){
			##Print chromosome bam file
			print STDERR "Creating $chr bam files....";
			my $chrBam = "samtools view -h -b $opts{out_prefix}.bam $chr >$opts{out_directory}/chrom_alignments/$opts{out_prefix}_$chr.bam";
			print "$chrBam\n";
			my $err = `$chrBam`;
			print STDERR "$err done\n";				
		}
		
		my $output =  "$opts{out_directory}/chrom_alignments/$config_opts{IMBALANCE_DIR_NAME}/$opts{out_prefix}_$chr";
		
		##Initial heterozygous site and imbalance detection
		print STDERR "Detect heterozygous site and initial imbalances....";
		my $initialHet = "samtools view $opts{out_directory}/chrom_alignments/$opts{out_prefix}_$chr.bam | $PERLPATH $config_opts{PATH_TO_AAALIGNER_BIN}/detectAllelicImbalance.pl -output $output.initial -minimumCoverage $config_opts{MINIMUM_READ_DEPTH}  -mismatches $config_opts{MISMATCHES} -format gsnap -alignedHets $config_opts{HET_FILE} $forceMismatches $completeImbalances";
		print "$initialHet\n";
		$err =`$initialHet`;
		print STDERR "$err done\n";

		##Align reads with alternate allele
		print STDERR "Aligning alternate reads....";
		my $alignAltReads = "gsnap --terminal-threshold=10  -A sam --read-group-library=$config_opts{READ_GROUP_LIBRARY} --read-group-platform=$config_opts{READ_GROUP_PLATFORM} --read-group-id=$config_opts{READ_GROUP_ID} -d $config_opts{GENOME_NAME} -D $config_opts{GENOME_DIR} -v $config_opts{SNP_DATABASE_NAME} -k $config_opts{KMER} --basesize=$config_opts{BASESIZE} --sampling=$config_opts{SAMPLING} -m $config_opts{MISMATCHES} -n $config_opts{NUMBER_READS_OUTPUT} -i $config_opts{INDEL_PENALTY}  --query-unk-mismatch=$config_opts{QUERY_UNK_MISMATCH} --genome-unk-mismatch=$config_opts{GENOME_UNK_MISMATCH}  -t $config_opts{GSNAP_ALTREAD_THREADS}  --trim-mismatch-score=$config_opts{TRIM_MISMATCH_SCORE} $output.initial.aider.$config_opts{MINIMUM_READ_DEPTH}.altReads.fa  >$output.initial.aider.$config_opts{MINIMUM_READ_DEPTH}.altReads.sam"; 
		print "$alignAltReads\n";
		$err = `$alignAltReads`;
		print STDERR "$err done\n";
		
		##Calculate read mappability
		print STDERR "Calculating read mappability....";
		my $determineMappability = "$PERLPATH  $config_opts{PATH_TO_AAALIGNER_BIN}/determineMappability.pl $output.initial.aider.$config_opts{MINIMUM_READ_DEPTH}.altReads";
		print "$determineMappability\n";
		$err = `$determineMappability`; 
		print STDERR "$err done\n";

		##Final heterozygous site and imbalance detection
		print STDERR "Detect heterozygous site and final imbalances....";
		my $finalHet = "samtools view $opts{out_directory}/chrom_alignments/$opts{out_prefix}_$chr.bam | $PERLPATH  $config_opts{PATH_TO_AAALIGNER_BIN}/detectAllelicImbalance.pl -output $output.final -minimumCoverage $config_opts{MINIMUM_READ_DEPTH}  -mismatches $config_opts{MISMATCHES} -format gsnap -alignedHets $config_opts{HET_FILE} $forceMismatches -mappability  $output.initial.aider.$config_opts{MINIMUM_READ_DEPTH}.altReads.mappability.txt $completeImbalances";
		print "$finalHet\n";
		$err =`$finalHet`;
		print STDERR "$err done\n";
		push (@chromosome_files, "$output.final.aider.$config_opts{MINIMUM_READ_DEPTH}.heterozygotes.txt");
	}
	my $mergeHets = "$config_opts{PATH_TO_AAALIGNER_BIN}/mergeAIDerOutput.pl -pvalue $config_opts{IMBALANCE_PVALUE} -filter $config_opts{HET_FILE} -output $opts{out_directory}/$opts{out_prefix}.aider.$config_opts{MINIMUM_READ_DEPTH} @chromosome_files";
	
}

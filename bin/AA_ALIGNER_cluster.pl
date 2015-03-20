use Getopt::Long;
use Cwd;

my %config_opts=();
GetOptions(\%opts, qw(all create_reference filter_fastq align filter_alignment find_imbalance out_directory:s temp_directory:s out_prefix:s configuration_file:s cluster_configuration=s help));

my @files = @ARGV;
my $JAVA_PATH = "java";
my $PERL_PATH = "perl";
my $PYTHON_PATH = "python";
my $RUN = 1; ##For debugging purposes. Set to 0 to see commmands but not to run them.
my $cur_dir = getcwd;

if ($opts{help} || $#files<0 ||(!$opts{filter_fastq} && !$opts{align} && !$opts{filter_alignment} && !$opts{allelic_imbalance} && !$opts{create_reference} &&! $opts{find_imbalance} && !$opts{all})){
	my $usage = "\nAA-ALIGNER v1.0\n\nperl AA_ALIGNER_cluster.pl [options]... [fastq files(s)]\n\nGENERAL OPTIONS\n-configuration_file <configuration file>\t\tconfiguration file (required)\n-cluster_configuration <cluster configuration file>\tconfiguration file for cluster settings\n-out_directory\t<output directory>\t\t\tdirectory for permanent files (default=./)\n-temp_directory <temporary directory>\t\t\tdirectory for intermediate files (default=./)\n-out_prefix <output prefix>\t\t\t\tprefix for output bam file and imbalance files\n\nRUN OPTIONS:\n-create_reference\tcreate a sample-specific reference genome (default=OFF)\n-filter_fastq\t\tfilters fastq files before alignment to convert base qualities to Sanger format, trim sequences (if necessary), and removes low quality sequence reads (default=OFF)\n-align\t\t\taligns sequences to specified genomes (default=OFF)\n-filter_alignment\tfilters alignments to remove reads that map multiple sites,fall in blacklist regions, and are possible PCR duplicates (default = OFF)\n-find_imbalance\t\tidentifies sites allelic imbalance in sequence alignment (default=OFF)\n-all\t\t\truns entire alignment pipeline filter_fastq, align, filter_alignments, and find_imbalance (default=ON)\n";
	print "$usage\n";
}else{
	#if (($opts{filter_fastq} || $opts{align} || $opts{filter_alignment} || $opts{allelic_imbalance} || $opts{create_reference} ||$opts{find_imbalance}){
	#	$opts{all} = 0;
	#}
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

	$JAVAPATH = $config_opts{JAVA_PATH}  if ($config_opts{JAVA_PATH} ne "");

	$PERLPATH = $config_opts{PERL_PATH} if ($config_opts{PERL_PATH} ne "");

	$PYTHONPATH = $config_opts{PYTHON_PATH} if ($config_opts{PYTHON_PATH} ne "");

	my %cluster_opt= ();
	open(IN, $opts{cluster_configuration});
	while(<IN>){
		chomp;
		if(substr($_,0,1) ne "#"){
			my @cols = split /\t/, $_;
			my @parts = split /=/, $cols[0];
			$cluster_opts{$parts[0]}=$parts[1];
		}
	}
	close(IN);
	
	determineKmer();
	my $previous_job = "";
	$previous_job = createReference() if ($opts{create_reference});

	my @filter_alignment_jobs = ();
	foreach my $file (@files){
		print "$file\n\n";
		print STDERR "$file\n\n";
		my @filter_jobs = ();
		push (@filter_jobs, $previous_job) if ($previous_job ne '');

		push(@filter_jobs, filterFastq($file)) if ($opts{filter_fastq} || $opts{all});
		if($opts{align} || $opts{all}){
			push (@align_jobs, align($file, @filter_jobs)) if ($opts{align} || $opts{all});
		}else{
			@align_jobs = @filter_jobs;
		}
		if($opts{filter_alignment} || $opts{all}){
			push(@filter_alignment_jobs,filterAlignment($file, @align_jobs)) if ($opts{filter_alignment} || $opts{all});
		}else{
			push (@filter_alignment_jobs, @align_jobs);
		}
	}	

	my $previous_job = "";
	$previous_job = mergeFilteredAlignment(@filter_alignment_jobs) if ($opts{filter_alignment} || $opts{all});
	findImbalance($previous_job) if ($opts{find_imbalance} || $opts{all});
}

##Pre-alignment fastqc filtering
sub filterFastq{
	my $file = shift;	
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
	my $submit_command = createSubmitCommand("FASTQC", "$prefix"."_fastqc","");
	my $fastqc_cmd = "$submit_command \"fastqc --noextract --outdir $opts{temp_directory} -t $config_opts{FASTQC_THREADS}  $config_opts{FASTQC_OPTIONS} $file\"";
	my $log = `$fastqc_cmd` if ($RUN ==1);
	print "$fastqc_cmd\n";
	print STDERR "$log done\n";
	my $previous_job = $prefix ."_fastqc";	

	##Check base quality score format
	my @parts = split /\//, $file;
	my $dir ="";
	for (my $i=0; $i<$#parts; $i++){
		$dir = $dir . "$parts[$i]\/";
	}
	my $base = $parts[$#parts];
	$base=~s/\.gz//;
	$base=~s/\.fastq//;
		
	print STDERR "Checking base quality score format ....";
	$submit_command = createSubmitCommand("", $prefix . "_check_format", $previous_job); 
	my $quality_command =  "$submit_command \"sh $config_opts{PATH_TO_AAALIGNER_BIN}/qualityFormat2sanger.sh $dir $parts[$#parts] $opts{temp_directory} $config_opts{GALAXY_LIBRARY_PATH} $config_opts{GROOMER_PATH}\"";
	print "$quality_command\n";
	$err = `$quality_command` if ($RUN==1);
	my $previous_job = $prefix ."_check_format";
	print STDERR "$err...done";
		
	##Unzip to temporary file if necessary
	my $base=$parts[$#parts];
	$base =~ s/.gz//;
	my $unzipped_file = $base;
	if (-e "$dir"."$base.gz"){
		print STDERR "Unzipping file ....";
		$submit_command = createSubmitCommand("",$prefix."_unzip", $previous_job);
		my $gunzip_cmd ="$submit_command \"gunzip -c $opts{temp_directory}/$base.gz > $opts{temp_directory}/$base\"";
		$previous_job = $prefix. "_unzip";
		`$gunzip_cmd` if ($RUN==1);
		print "$gunzip_cmd\n";
		$unzipped_file = "$opts{temp_directory}/$base";
		print STDERR "done\n";
	}

	print STDERR "Trimming fastq bases ....";
	##FASTX Trimmer to trim bases to length in configuration file
	if($config_opts{TRIM_LAST_BASE} eq "SEQUENCE_LENGTH"){
		$config_opts{TRIM_LAST_BASE}=$config_opts{READ_LENGTH};
	}
	$submit_command = createSubmitCommand("FASTX", $prefix."_fastx_trim", $previous_job);
	my $fastx_trim_cmd = "$submit_command \"fastx_trimmer -f $config_opts{TRIM_FIRST_BASE} -l $config_opts{TRIM_LAST_BASE} -Q 33 -o $opts{temp_directory}/trimmed_$prefix.fastq $config_opts{FASTX_OPTIONS} -i $unzipped_file\"";
	$err = `$fastx_trim_cmd` if ($RUN == 1);
	print STDERR "$err";
	print "$fastx_trim_cmd\n";
	print STDERR "done\n";
		
	##FASTQ_QUALITY_FILTER Trimmer to trim out low quality sequences
	print STDERR "Filtering out low quality bases ....";
	$submit_command = createSubmitCommand("FASTQ_QUALIY",$prefix."_fastq_qual", $prefix."_fastx_trim");
	my $fast_qual_cmd = "$submit_command \"fastq_quality_filter -Q 33 -p $config_opts{FASTQ_QUALITY_FILTER_P} -q $config_opts{FASTQ_QUALITY_FILTER_Q} $config_opts{FASTQ_QUALITY_OPTIONS} -i $opts{temp_directory}/trimmed_$prefix.fastq -o $opts{out_directory}/$base.filtered\"";
	$err = `$fast_qual_cmd` if ($RUN==1);
	print "$fast_qual_cmd\n";
	print STDERR "$err";
	print STDERR "done\n";
	print STDERR "Compressing file ....";

	$submit_command = createSubmitCommand("", $prefix."_gzip", $prefix."_fastq_qual");
	my $gzip ="$submit_command \"gzip $opts{out_directory}/$base.filtered\"";
	$err = `$gzip`;
	print "$gzip\n";
	print STDERR "$err done\n";
	push (@filtered_fastqs, "$opts{out_directory}/$base.filtered.gz");
	return $prefix. "_gzip";
}
sub determineKmer{
	print STDERR "Calculating k-mer ....";
	my $mm = $config_opts{MISMATCHES};
	my $k = $config_opts{READ_LENGTH}/($mm+1);
	if ($k<10){
		print STDERR "$k:K-mer < 10 is not reocommended. Setting k-mer to 10, but this can cause inaccurate sequeunce alignment. Try a lower mismatch threshold or use set FORCE_KMER=1\n";
		$k = 10 if (! $confiq_opts{FORCE_KMER});
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
	my @previous_jobs = ();	
	
	my $create_reference_job = createReferenceFa($genome_dir);
	push (@previous_jobs, buildGenome($genome_dir, $create_reference_job));
	push (@previous_jobs, createIITFile($genome_dir));
	$previous_job = indexSNPs($genome_dir,@previous_jobs);

}

sub createReferenceFa{
	my $genome_dir = shift;
	if($config_opts{SPLIT_BY_CHROM}){
	
		#Get prefix for split chr filename
		my @hets_path = split /\//, $config_opts{HET_FILE};
		my @previous_jobs = ();
		if ($#hets_path >= 0){
			$hets_path[$#hets_path] =~ s/\.bed$//;
			push (@previous_jobs,  splitFilesByChrom($config_opts{HET_FILE}, $hets_path[$#hets_path]));
		}
		my @replace_path = split /\//, $config_opts{REPLACE_FILE};
		if ($#replace_path >= 0){
			$replace_path[$#replace_path] =~ s/\.bed$//;
			push (@previous_jobs, splitFilesByChrom($config_opts{REPLACE_FILE}, $replace_path[$#het_path]));
		}

		#Create chromosome array
		my @chrom = `ls $config_opts{REF_FASTQ_DIR}/chr*.fa*`;	
		
		foreach my $file (@chrom){
			$file =~ m/chr(.*).fa/;
			my $chr = $1;
			my $createChrFa = "";
			if($config_opts{REFERENCE_FILE_COMPRESSED}){
				$createChrFa="gunzip -c $config_opts{REF_FASTQ_DIR}/chr$chr.fa.gz |$PERLPATH $config_opts{PATH_TO_AAALIGNER_BIN}/replaceReferenceAlleles.pl -output $genome_dir/chr$chr.fa -hets $opts{temp_directory}/$hets_path[$#hets_path].chr$chr.bed -replace $opts{temp_directory}/$replace_path[$#replace_path].chr$chr.bed  >>$genome_dir/$config_opts{GENOME_NAME}.replace.txt 2>>$genome_dir/$config_opts{GENOME_NAME}.replaced.log"; 
			}else{
				$createChrFa="cat $config_opts{REF_FASTQ_DIR}/chr$chr.fa |$PERLPATH $config_opts{PATH_TO_AAALIGNER_BIN}/replaceReferenceAlleles.pl -output $genome_dir/chr$chr.fa -hets $opts{temp_directory}/$het_path[$#het_path].chr$chr.bed -replace $opts{temp_directory}/$replace_path[$#replace_path].chr$chr.bed   >>$genome_dir/$config_opts{GENOME_NAME}.replace.txt 2>>$genome_dir/$config_opts{GENOME_NAME}.replaced.log "; 
			}
			print STDERR "Creating chr$chr.fa ...";
			$submit_command = createSubmitCommand("REFERENCE_CREATION", "create_genome_chr$chr", @previous_jobs);
			my $err = `$submit_command \"$createChrFa\"` if ($RUN == 1);
			print STDERR "$err\n";
			print "$submit_command \"$createChrFa\"\n";	
			print STDERR "done\n";
			@previous_jobs= ("create_genome_chr$chr");
		}
		$previous_job = $previous_jobs[0];
				
	}else{
		my @fastq_file;
		my $createChrFa = "";
		
		if($config_opts{REFERENCE_FILE_COMPRESSED}){
				$createFa="gunzip -c $config_opts{REF_FASTQ_DIR}/$config_opts{REF_FASTQ} |$PERLPATH $config_opts{PATH_TO_AAALIGNER_BIN}/replaceReferenceAlleles.pl -output $genome_dir/$config_opts{GENOME_NAME} -hets $config_opts{HET_FILE} -replace $config_opts{REPLACE_FILE} -append  >>$genome_dir/$config_opts{GENOME_NAME}.replace.txt 2>>$genome_dir/$config_opts{GENOME_NAME}.replaced.log"; 	
		}else{
				$createFa="cat $config_opts{REF_FASTQ_DIR}/$config_opts{REF_FASTQ} |$PERLPATH $config_opts{PATH_TO_AAALIGNER_BIN}/replaceReferenceAlleles.pl -output $genome_dir/$config_opts{GENOME_NAME} -hets $config_opts{HET_FILE} -replace $config_opts{REPLACE_FILE} -append >>$genome_dir/$config_opts{GENOME_NAME}.replace.txt 2>>$genome_dir/$config_opts{GENOME_NAME}.replaced.log"; 	
		}
		print STDERR "Creating reference ...";
		$submit_command = createSubmitCommand("REFERENCE_CREATION", "create_genome", @previous_jobs);
		my $err = `$submit_command "$createFa"` if ($RUN == 1);
		print STDERR "$err\n";
		print "$submit_command \"$createFa\"\n";	
		print STDERR "done\n";
		$previous_job= "create_genome";
	}
	return $previous_job;
}

sub splitFilesByChrom{
	my $file =shift;
	my $prefix = shift;
	my $previous_job = "";
	print STDERR "Splitting $file to $prefix.<CHROM>.bed ...";
	$submit_command = createSubmitCommand("", $prefix."_split", "");
	my $split_cmd = "$PERLPATH $config_opts{PATH_TO_AAALIGNER_BIN}/splitBedByChrom.pl $file $opts{temp_directory}/$prefix"; 	
	my $err = `$submit_command \"$split_cmd\"` if ($RUN == 1);
	print STDERR "$err\n";
	print "$submit_command \"$split_cmd\"\n";
	print STDERR "done\n";
	return $prefix . "_split";
}

sub createGSNAPFile{
	my @hets_path = split /\//, $config_opts{HET_FILE};
	$hets_path[$#hets_path] =~ s/\.bed//;
	print STDERR "Creating gsnap file ...";	
	$submit_command = createSubmitCommand("", $hets_path[$#hets_path]."_gsnapFile","");
	my $create_file = "$submit_command \"$PERLPATH $config_opts{PATH_TO_AAALIGNER_BIN}/bedToGSNAP.pl $config_opts{HET_FILE} \>$opts{temp_directory}\/$hets_path[$#hets_path].gsnap.txt\"";
	my $err = `$create_file` if ($RUN == 1);
	print STDERR "$err\n";
	print "$create_file\n";
	print STDERR "done\n";
	#return $hets_path[$#hets_path];
	return $hets_path[$#hets_path] . "_gsnapFile";
}

sub buildGenome{
	my $genome_dir = shift;
	my $previous_job = shift;
	print STDERR "Building genome ...";
	$submit_command = createSubmitCommand("GMAP_BUILD", $config_opts{GENOME_NAME}."_build",$previous_job);
	my $buildGenome = "$submit_command \"gmap_build -T . -D $config_opts{GENOME_DIR} -d $config_opts{GENOME_NAME} -k $config_opts{KMER} -b $config_opts{BASESIZE}  -q $config_opts{SAMPLING} $config_opts{GMAP_BUILD_OPTIONS} $genome_dir/*.fa\"";
	my $dir = "$config_opts{GENOME_DIR}/$config_opts{GENOME_NAME}";
	#chdir $dir;
	my $err = `$buildGenome` if ($RUN == 1);
	#chdir $cur_dir;
	print STDERR "$err\n";
	print "$buildGenome\n";  
	print STDERR "done\n";
	return $config_opts{GENOME_NAME} ."_build";
} 

sub createIITFile{
	my $genome_dir = shift;
	my $previous_job=createGSNAPFile();
	print STDERR "Converting gsnap to iit_file\n";
	my @hets_path = split /\//, $config_opts{HET_FILE};
	$hets_path[$#hets_path] =~ s/\.bed//;
	$submit_command = createSubmitCommand("", $hets_path[$#hets_path]."_toIIT", $previous_job);
	my $iit_cmd = "$submit_command \"cat $opts{temp_directory}/$hets_path[$#hets_path].gsnap.txt |iit_store $config_opts{IIT_STORE_OPTS} -o  $genome_dir/$config_opts{SNP_DATABASE_NAME}.iit\"";
	my $err =`$iit_cmd` if ($RUN == 1);
	print STDERR "$err\n";
	print "$iit_cmd\n";
	print STDERR "done\n";	 
	return $hets_path[$#hets_path] . "_toIIT";
}

sub indexSNPs{
	my $genome_dir = shift;
	my @previous_jobs = @_;
	print STDERR "Indexing SNPs ...";
	my @hets_path = split /\//, $config_opts{HET_FILE};
	$hets_path[$#hets_path] =~ s/\.bed//;
	$submit_command = createSubmitCommand("SNPINDEX", $config_opts{SNP_DATABASE_NAME}."_snp_index", @previous_jobs);
	my $indexSNPs = "$submit_command \"snpindex  -D $config_opts{GENOME_DIR} -d $config_opts{GENOME_NAME} -k $config_opts{KMER} -b $config_opts{BASESIZE}  -q $config_opts{SAMPLING} -v $config_opts{SNP_DATABASE_NAME} $config_opts{SNPINDEX_OPTIONS} $genome_dir/$config_opts{SNP_DATABASE_NAME}.iit\"";
	my $dir = "$config_opts{GENOME_DIR}/$config_opts{GENOME_NAME}";
	#chdir $dir;
	my $err = `$indexSNPs` if ($RUN == 1);
	#chdir $cur_dir;
	print STDERR "$err\n";
	print "$indexSNPs\n";
	print STDERR "done\n";	 
	return $config_opts{SNP_DATABASE_NAME} . "_snp_index";
}

sub align{
	my $file = shift;
	my @previous_jobs = @_;
	print STDERR "\nAligning sequences\n";
	print  "\nAligning sequences\n";
	
	my @parts = split /\//, $file;
	my $base = $parts[$#parts];
	$base =~ s/.gz//;
	my $filename = "";
	if($opts{filter_fastq} || $opts{all}){
		$filename = "$opts{out_directory}/$base.filtered.gz";
		#$base = "$base.filtered.gz";
	}else{
		$filename = "$file";
	}
	my @file_parts = split /\./, $filename;
	my $gunzip = "";
	if($file_parts[$#file_parts] eq "gz"){
		$gunzip = "--gunzip";
	}
	if($config_opts{INDEL_PENALTY} eq "NOINDELS"){
		$config_opts{INDEL_PENALTY}= $config_opts{MISMATCHES} + 1;
	}
	print STDERR "Aligning reads....";
	$submit_command = createSubmitCommand("GSNAP_ALIGN", $base."_align", @previous_jobs);
	
	my $align = "$submit_command \"gsnap --terminal-threshold=10  -A sam --read-group-library $config_opts{READ_GROUP_LIBRARY} --read-group-platform=$config_opts{READ_GROUP_PLATFORM} --read-group-id=$config_opts{READ_GROUP_ID} -d $config_opts{GENOME_NAME} -D $config_opts{GENOME_DIR} -v $config_opts{SNP_DATABASE_NAME} -k $config_opts{KMER} --basesize=$config_opts{BASESIZE} --sampling=$config_opts{SAMPLING} -m $config_opts{MISMATCHES} -n $config_opts{NUMBER_READS_OUTPUT} -i $config_opts{INDEL_PENALTY}  --query-unk-mismatch=$config_opts{QUERY_UNK_MISMATCH} --genome-unk-mismatch=$config_opts{GENOME_UNK_MISMATCH}  -t $config_opts{GSNAP_THREADS}  --trim-mismatch-score=$config_opts{TRIM_MISMATCH_SCORE} $gunzip $config_opts{GSNAP_OPTIONS} $filename  >$opts{out_directory}/$base.sam\""; 
	my  $err = `$align` if ($RUN == 1);
	print "$align\n";
	print STDERR "$err\n";
	print STDERR "done\n";
	#push (@alignment_sams, "$opts{temp_directory}/$base.sam");
	my $previous_job = $base . "_align";
	return $previous_job;

}
	
sub filterAlignment{
	my $file = shift;
	my @previous_jobs = shift;
	print STDERR "\nPost-alignment filtering\n\n";
	print  "\nPost-alignment filtering\n\n";
	$file =~ s/.gz//;
	my @parts = split /\//, $file;
	my $base = $parts[$#parts];
	my $sam = "$opts{out_directory}/$base.sam";
	
	print STDERR "Filtering to keep unique reads...."; 
	$submit_command = createSubmitCommand("", $base."_filterUnique", @previous_jobs);
	my $filter_unique = "$submit_command \"$PERLPATH $config_opts{PATH_TO_AAALIGNER_BIN}/printUniqueSam.pl --header -format gsnap -output $opts{temp_directory}/$base $sam\"";
	print "$filter_unique\n";
	my $err = `$filter_unique` if ($RUN == 1);
	print STDERR "$err done\n";
		
	##Filters out any alignments overalapping regions in specified blacklist (bed format) file
	print STDERR "Filtering blacklist regions....";
	$submit_command = createSubmitCommand("FILTER_DUPLICATES", $base . "_filterBlacklist" ,$base."_filterUnique");
	my $blacklist = "$submit_command \"$PERLPATH $config_opts{PATH_TO_AAALIGNER_BIN}/filterBlacklistSam.pl -output $opts{temp_directory}/$base.unique.blacklist -filter $config_opts{BLACKLIST_FILE} -input $opts{temp_directory}/$base.unique.sam\"";

	print "$blacklist\n";
	$err = `$blacklist` if ($RUN == 1);
	print STDERR "$err done\n";
	
	##Converts to bam
	print STDERR "Converting to sorted bam....";
	$submit_command = createSubmitCommand("", $base . "_blacklistToBam" ,$base."_filterBlacklist");
	my $toSortedBam = "$submit_command \"samtools view -Sb -h $opts{temp_directory}/$base.unique.blacklist.sam |samtools sort - $opts{temp_directory}/$base.unique.blacklist\"";	
	print "$toSortedBam\n";
	$err = `$toSortedBam` if ($RUN == 1);
	print STDERR "$err done\n";
	
	##Remove duplicate reads 
	print STDERR "Removing duplicates....";
	$submit_command = createSubmitCommand("MARKDUPLICATES", $base . "_markDuplicates", $base . "_blacklistToBam");
	my $removeDups = "$submit_command \"$JAVAPATH -Xmx4g -jar $config_opts{PATH_TO_PICARD}/MarkDuplicates.jar REMOVE_DUPLICATES=TRUE METRICS_FILE=$opts{out_directory}/$base.unique.blacklist.markDuplicates.metric.txt OUTPUT=$opts{temp_directory}/$base.unique.blacklist.markDuplicates.bam VALIDATION_STRINGENCY=LENIENT INPUT=$opts{temp_directory}/$base.unique.blacklist.bam $config_opts{MARK_DUPLICATES_OPTIONS} \"";
	print "$removeDups\n";	
	$err = `$removeDups` if ($RUN == 1);
	print STDERR "$err done\n";
	return $base . "_markDuplicates";
}

sub mergeFilteredAlignment{ 
	my @previous_jobs = @_;
	##Merge replicate files together
	print STDERR "\nMerging replicate bam files....";
	print  "\nMerging replicate bam files\n\n";
	$submit_command = createSubmitCommand("", $opts{out_prefix}. "_mergeReplicateBams" ,@previous_jobs);
	my $merge_str = "$submit_command \"$JAVAPATH -Xmx4g -jar $config_opts{PATH_TO_PICARD}/MergeSamFiles.jar OUTPUT=$opts{out_directory}/$opts{out_prefix}.bam USE_THREADING=$config_opts{PICARD_USE_THREADING} $config_opts{MERGE_BAM_OPTIONS}\"";

	foreach my $file (@files){
		$file =~ s/.gz//;
		my @parts = split /\//, $file;
		my $base = $parts[$#parts];
		my $rep = "$opts{temp_directory}/$base.unique.blacklist.markDuplicates.bam";
		$merge_str .= " INPUT=$rep";
	}
	print "$merge_str\n";
	$err = `$merge_str` if ($RUN == 1);
	print STDERR "$err done\n";
	return $opts{out_prefix} . "_mergeReplicateBams";
}
		
sub findImbalance{
	my $previous_job = shift;
	##Initial round of imbalance detection	
	print STDERR "\nDetecting Imbalances\n\n";
	print  "\nDetecting Imbalances\n\n";
	my $submit_command = "";
	if($config_opts{INDEX_BAM}){
		print STDERR "Indexing bam file....";
		$submit_command = createSubmitCommand("", $opts{out_prefix} . "_indexBam" ,$previous_job);
		my $index = "$submit_command samtools index $opts{out_directory}/$opts{out_prefix}.bam";
		print "$index\n";
		my $err = `$index` if ($RUN == 1);
		print STDERR "$err done\n";
		$previous_job = $opts{out_prefix} . "_indexBam";
	}
	
	#Create output directories
	print STDERR "Making directories....";
	`mkdir $opts{out_directory}/chrom_alignments`;
	`mkdir $opts{out_directory}/chrom_alignments/$config_opts{IMBALANCE_DIR_NAME}`;
	
	my $forceMismatches = "";
	$forceMismatches = "-forceMismatches" if ($config_opts{HET_MM_ONLY} == 1);	
	
	##Gets chromosome information from bam header
	#my @header = `samtools view -H $opts{out_directory}/$opts{out_prefix}.bam`;
	my @header = `cut -f1 $config_opts{GENOME_DIR}/$config_opts{GENOME_NAME}/$config_opts{GENOME_NAME}.chromosome`;

	my @chrom;	
	foreach my $line (@header){
		chomp $line;	
		#if($line=~m/@SQ\tSN:(.*)\t/){
			#push(@chrom, $1);
		#}	
		push(@chrom, $line);
	}
	if($config_opts{INDEL_PENALTY} eq "NOINDELS"){
		$config_opts{INDEL_PENALTY}= $config_opts{MISMATCHES} + 1;
	}
	
	my @chromosome_files = ();
	foreach my $chr (@chrom){
		my $previous_chrom_job = $previous_job;
		if($config_opts{INDEX_BAM}){
			##Print chromosome bam file
			print STDERR "Creating $chr bam files....";
			$submit_command = createSubmitCommand("", $opts{out_prefix} . "_index_$chr" ,$previous_job);
			my $chrBam = "$submit_command \"samtools view -h -b $opts{out_prefix}.bam $chr >$opts{out_directory}/chrom_alignments/$opts{out_prefix}_$chr.bam\"";
			print "$chrBam\n";
			my $err = `$chrBam` if ($RUN == 1);
			print STDERR "$err done\n";				
			$previous_chrom_job = $opts{out_prefix} . "_index_$chr";
		}
		
		my $output =  "$opts{out_directory}/chrom_alignments/$config_opts{IMBALANCE_DIR_NAME}/$opts{out_prefix}_$chr";
		
		##Initial heterozygous site and imbalance detection
		print STDERR "Detect heterozygous site and initial imbalances....";
		$submit_command = createSubmitCommand("AIDER_SUBMIT", $opts{out_prefix} . "_initial_imbalance_$chr" ,$previous_chrom_job);
		my $initialHet = "$submit_command \"samtools view $opts{out_directory}/chrom_alignments/$opts{out_prefix}_$chr.bam | $PERLPATH $config_opts{PATH_TO_AAALIGNER_BIN}/AIDer.pl -output $output.initial -minimumCoverage $config_opts{MINIMUM_READ_DEPTH}  -mismatches $config_opts{MISMATCHES} -format gsnap -alignedHets $config_opts{HET_FILE} $forceMismatches\"";
		print "$initialHet\n";
		$err =`$initialHet` if ($RUN == 1);
		print STDERR "$err done\n";

		##Align reads with alternate allele
		print STDERR "Aligning alternate reads....";
		$submit_command = createSubmitCommand("AIDER_ALIGN", $opts{out_prefix} . "_align_$chr", $opts{out_prefix} . "_initial_imbalance_$chr");
		my $alignAltReads = "$submit_command \"gsnap --terminal-threshold=10  -A sam --read-group-library $config_opts{READ_GROUP_LIBRARY} --read-group-platform $config_opts{READ_GROUP_PLATFORM} --read-group-id $config_opts{READ_GROUP_ID} -d $config_opts{GENOME_NAME} -D $config_opts{GENOME_DIR} -v $config_opts{SNP_DATABASE_NAME} -k $config_opts{KMER} --basesize $config_opts{BASESIZE} --sampling $config_opts{SAMPLING} -m $config_opts{MISMATCHES} -n $config_opts{NUMBER_READS_OUTPUT} -i $config_opts{INDEL_PENALTY}  --query-unk-mismatch=$config_opts{QUERY_UNK_MISMATCH} --genome-unk-mismatch=$config_opts{GENOME_UNK_MISMATCH}  -t $config_opts{GSNAP_ALTREAD_THREADS}  --trim-mismatch-score=$config_opts{TRIM_MISMATCH_SCORE} $output.initial.aider.$config_opts{MINIMUM_READ_DEPTH}.altReads.fa  >$output.initial.aider.$config_opts{MINIMUM_READ_DEPTH}.altReads.sam\""; 
		print "$alignAltReads\n";
		$err = `$alignAltReads` if ($RUN == 1);
		print STDERR "$err done\n";
		
		##Calculate read mappability
		print STDERR "Calculating read mappability....";
		$submit_command = createSubmitCommand("AIDER_MAPPABILITY", $opts{out_prefix} . "_mappability_$chr", $opts{out_prefix} . "_align_$chr");
		my $determineMappability = "$submit_command \"$PERLPATH  $config_opts{PATH_TO_AAALIGNER_BIN}/determineMappability.pl $output.initial.aider.$config_opts{MINIMUM_READ_DEPTH}.altReads\"";
		print "$determineMappability\n";
		$err = `$determineMappability` if ($RUN == 1); 
		print STDERR "$err done\n";

		##Final heterozygous site and imbalance detection
		print STDERR "Detect heterozygous site and final imbalances....";
		$submit_command = createSubmitCommand("AIDER_SUBMIT", $opts{out_prefix} . "_final_imbalance_$chr", $opts{out_prefix} . "_mappability_$chr");
		my $finalHet = "$submit_command \"samtools view $opts{out_directory}/chrom_alignments/$opts{out_prefix}_$chr.bam | $PERLPATH  $config_opts{PATH_TO_AAALIGNER_BIN}/AIDer.pl -output $output.final -minimumCoverage $config_opts{MINIMUM_READ_DEPTH}  -mismatches $config_opts{MISMATCHES} -format gsnap -alignedHets $config_opts{HET_FILE} $forceMismatches -mappability  $output.initial.aider.$config_opts{MINIMUM_READ_DEPTH}.altReads.mappability.txt\"";
		print "$finalHet\n";
		$err =`$finalHet` if ($RUN == 1);
		print STDERR "$err done\n";
		push (@imbalance_jobs, $opts{out_prefix} . "_final_imbalance_$chr");
		push (@chromosome_files, "$output.final.aider.$config_opts{MINIMUM_READ_DEPTH}.heterozygotes.txt");
	}
	print STDERR "Merge heterozygous sites....";
	$submit_command = createSubmitCommand("", $opts{out_prefix} . "_mergeHeterozygousSites", @imbalance_jobs);	
	my $mergeHets = "$submit_command \"$PERLPATH $config_opts{PATH_TO_AAALIGNER_BIN}/mergeAIDerOutput.pl -pvalue $config_opts{IMBALANCE_PVALUE} -filter $config_opts{HET_FILE} -output $opts{out_directory}/$opts{out_prefix}.aider.$config_opts{MINIMUM_READ_DEPTH} @chromosome_files\"";
	print "$mergeHets\n";
	$err = `$mergeHets` if ($RUN == 1);
	print STDERR "$err done\n";	
}
sub createSubmitCommand{
	my $program = shift;
	my $job_name = shift;
	my @previous_jobs = @_;
	my $option_name = $program . "_SUBMIT_OPTIONS";
	my $other_opts = $cluster_opts{$option_name};
	my $submit = "";
	my $queue = "-q $cluster_opts{QUEUE}";
	$queue = "" if ($cluster_opts{QUEUE} eq "");
	
	if($cluster_opts{PLATFORM} eq "lsf"){
		
		my $dependency = "";
		
		if ($previous_jobs[0] ne ""){
			my $state = "done";
			$state = "ended" if ($program eq "AIDER_MAPPABILITY");
			$dependency = "-w \"$state(". $previous_jobs[0] .")";
			for (my $i = 1; $i<=$#previous_jobs; $i++){
				$dependency .= " && $state($previous_jobs[$i])"; 	
			}
			$dependency .= "\"";
		}
		$submit = "bsub -J $job_name $queue $dependency -o $job_name.out $other_opts"; 
	}
	return $submit;	
}

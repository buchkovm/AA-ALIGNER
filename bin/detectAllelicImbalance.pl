##################################################################
# Allelic Imbalance Detecter (AIDer) v1.0
# Detect allelic imbalance in an alignment bam file
#
# Martin Buchkovich
# January 2015
####################################################################
my %opts=(mismatches=>100);
use Math::CDF;
use Getopt::Long;
GetOptions(\%opts,qw(alignedHets=s forceMismatches  mismatches=s format=s output=s minimumCoverage=s mappability=s));

my $min_qual = 30;
my $mismatch_threshold = $opts{mismatches};
my $format = $opts{format} || "gsnap";
my $min_cov = $opts{minimumCoverage} || 5;
my $output = $opts{output} || "scan";
my $forceMismatches = $opts{forceMismatches} || 0;

open (HETOUT, ">$output.aider.$min_cov.heterozygotes.txt");

if(defined ($opts{mappability})){
	open(IN, $opts{mappability});
	while(<IN>){
		chomp;
		my @cols = split /\t/, $_;
		my $name = $cols[3];
		$name =~ s/_alt//;
		#my @parts = split /_/, $cols[3];
		#my $name = "$parts[0]";
		#for(my $i=1; $i<$#parts; $i++){
		#	$name .= "_$parts[$i]";
		#}
		$map{$name}=$cols[5];	
	}
}

if($forceMismatches){
	open(IN, $opts{alignedHets});
	while(<IN>){
		chomp;
		my @cols = split /\t/, $_;
		$align_hets{$cols[0]}{$cols[2]}="$cols[4]$cols[5]";
	}
}
if(!defined($opts{mappability})){
	open (READOUT, ">$output.aider.$min_cov.reads.txt");
	open (ALTREADOUT, ">$output.aider.$min_cov.altReads.txt");
	open (ALTFAOUT, ">$output.aider.$min_cov.altReads.fa") ; 
}

my $binStart = 0;
my $binEnd = 0;
my $curChr = "";
my @reads_mm = ();
my @reads = ();
my %bin_snps=();
my $bin_hets=();

my $debug_pos = 15922; ##FOR DEBUGGING ONLY

while (<>){
	my @cols = split /\t/, $_;
	$chrom = $cols[2];
	##Find the edit distance from SAM
	$_=~m/NM:i:(\d+)/;
	my $edit_distance = $1;
	my $mismatch = $edit_distance;
	my $nh = 0;
	
	if($pos > ($debug_pos-49)){
		my $tmp =  "";
	}
	##Gets the number of mismatched (XW) and unique hits (NH) from GSNAP 
	if($format eq "gsnap"){
		$_=~m/XW:i:(\d+)/;
		$mismatch = $1;
		$_ =~ m/NH\:\w+\:(\d+)/;
		$nh = $1;
	}elsif($format eq "star"){
		$_ =~ m/NH\:\w+\:(\d+)/;
		$nh = $1;
	##Gets the number of mismatche unique hits (X0) from GSNAP 
	}elsif($format eq "bwa"){
		$_ =~ m/X0\:\w+\:(\d+)/;
		$nh = $1;
	}		
	
	##Gets the MD string 
	$_=~m/MD:Z:([ATCG0-9]*)/;
	my $md = $1; 

	##If mappability file is included looks up the mappability of the read
	my $mappable = 1;
	if(defined $map{$cols[0]}){
		$mappable = $map{$cols[0]};
	}
	if($mappable == 2){
		$tmp=0;
	}	
	##The binEnds when there are no more overlapping sequences or a new chromosome
	if($curChr ne $cols[2] || $cols[3] > $binEnd){
		##Checks snps in bin to see if they meet qualifications to be a het
		findBinHets(\%bin_snps, \%bin_hets);	
		
		##If only allowing mismatch at heterozygous sites
		if($forceMismatches){
			foreach my $read (@reads_mm){
				my $keepRead = forceMismatchAtHets($read, \%bin_snps);
				push(@reads, $read) if ($keepRead == 1);	
			}
			%bin_hets = ();
			findBinHets(\%bin_snps, \%bin_hets);	
		}else{
			push(@reads, @reads_mm);	
		}
		
		createMappabilityReads(\@reads,\%bin_hets) if (!defined($opts{mappability}));
		##Reset bin and chromosome variables
		printBinHets(\%bin_hets);
		$binStart = $cols[3];
		$curChr = $cols[2];
		%bin_snps=();
		%bin_hets=();
		@reads = ();
		@reads_mm = ();
	}

	if($nh==1 && $mismatch<=$mismatch_threshold && $mappable == 1){	##Aligns to unique location, is below mismatch threshold, and alternate read (if included) also aligns uniquely
		##Splits reads based on whether there is a mismatch
		if ($mismatch == 0){
			push (@reads, $_);
		}elsif($mismatch <= $opts{mismatches}){
			push (@reads_mm, $_);
		}
		
		findBinSNPs(\%bin_snps, $_, $edit_distance, $md);	
		
	}
	##Counts the allele of each position in the sequence

	##Extend the binEnd to the last position of the sequence
	my $readEnd = $cols[3] + length($cols[9])-1;
	$binEnd = $readEnd;
	
}
##Checks snps in bin to see if they meet qualifications to be a het
findBinHets(\%bin_snps, \%bin_hets);	
		
##If only allowing mismatch at heterozygous sites
if($forceMismatches){
	foreach my $read (@reads_mm){
		my $keepRead = forceMismatchAtHets($read, \%bin_snps);
		push(@reads, $read) if ($keepRead == 1);	
		findBinHets(\%bin_snps, \%bin_hets);	
	}
}else{
	push(@reads, @reads_mm);	
}
	
createMappabilityReads(\@reads,\%bin_hets) if (!defined($opts{mappability}));
##Reset bin and chromosome variables
printBinHets(\%bin_hets);

sub findBinSNPs{
	my $bin_snp_ref = shift;
	my $read = shift;
	my $edit_distance= shift;
	my $mdstring = shift;

	my @cols = split /\t/, $read;
	
	my @seq = split //, $cols[9];
	my @quals = split //, $cols[10];
	
	my $start = $cols[3];
	if($start > $debug_pos){
		$temp = "";
	}
	##Process MD:Z string
	my %alt_alleles = ();
	my @md_parts = split /([ATCG])/, $mdstring;	
	my $md_pos = 0;
	foreach my $x (@md_parts){
		if($x eq "A" || $x eq "T" || $x eq "C" || $x eq "G" || $x eq "N" ){ ##Different  from reference 
			$alt_alleles{$md_pos}=$x;
			$md_pos++;
		}else{
			$md_pos+=$x;
		}
	}

	my @cigar = split /([MIHDNSHP=X])/, $cols[5];
	my $cur_cnt=0;
	my $cur_pos=$start;
	my $cur_seq_ind = 0;	
	for (my $cigar_ind = 0; $cigar_ind<=$#cigar; $cigar_ind++){
	#if($cigar[1] eq "H"){ ##Hard clippling at start of read
	#	$start = $start + $cigar[0];
	#}
		my $part = $cigar[$cigar_ind];
		if($part eq "M" || $part eq "=" || $part eq "X"){
			my $end = $cur_pos+$cur_cnt;
			for (my $i=$cur_pos; $i<($end); $i++){
				if($i == $debug_pos){
					$tmp==0;
				}
				my $qual = ord($quals[$cur_seq_ind])-33;
				if(defined ($alt_alleles{$cur_seq_ind})){
					my $base=$seq[$cur_seq_ind];
					#if($cigar_ind < ($#cigar-2) && $cigar[$cigar_ind+2] eq "I"){
					#	for (my $k = 2; $k < $cigar[$cigar_ind+1]; $k++){
					#		$base .= $seq[$cur_seq_ind];
					#		$cur_seq_ind++;
					#	} 
					#	$cigar_ind = $cigar_ind + 2;
					#}
					addAlt($bin_snp_ref, $cur_pos, $base, $alt_alleles{$cur_seq_ind}) if ($qual >= $min_qual);	
					
				}else{
					addRef($bin_snp_ref, $cur_pos, $seq[$cur_seq_ind]) if ($qual >= $min_qual);
				}	
				$cur_seq_ind++;
				$cur_pos++;
			} 			
		}elsif($part eq "I"){	
			
		}elsif($part eq "D"){	
			
		}elsif($part eq "N"){	
			$cur_pos+=$cur_cnt;	
		}elsif($part eq "S"){	
			$cur_seq_ind += $cur_cnt;	
		}elsif($part eq "H"){
			##Do nothing
		}else{
			$cur_cnt = $part;
		}	
	}
}

sub addRef{
	my $bin_snp_ref = shift;
	my $pos = shift;
	my $base = shift;
	$bin_snp_ref->{$pos}{ref}=$base;
	$bin_snp_ref->{$pos}{ref_cnt}++;
}

sub addAlt{
	my $bin_snp_ref = shift;
	my $pos = shift;
	my $base = shift;
	my $ref = shift;
	$bin_snp_ref->{$pos}{ref}=$ref;
	$bin_snp_ref->{$pos}{alt}{$base}++;
}

sub findBinHets{
	my $bin_snp_ref = shift;
	my %bin_snp = %{$bin_snp_ref};
	my $bin_hets_ref=shift;
	foreach my $pos (sort {$a<=>$b} keys %bin_snp){
		my %alt = %{$bin_snp{$pos}{alt}};
		my $alt_cnt =0;
		my %alt_min;
		if($pos == $debug_pos){
			my $tmp =  "";
		}
		foreach my $alt_allele (sort {$alt{$b}<=>$alt{$a}} keys %alt){
			if($bin_snp{$pos}{alt}{$alt_allele} >= $min_cov){
				$alt_min{$alt_allele}=$bin_snp{$pos}{alt}{$alt_allele};
				$alt_cnt++;
				if($bin_snp{$pos}{ref_cnt} >= $min_cov){
					$bin_hets_ref->{$pos}{alt}{$alt_allele} = $bin_snp{$pos}{alt}{$alt_allele};
					$bin_hets_ref->{$pos}{ref} = $bin_snp{$pos}{ref};
					$bin_hets_ref->{$pos}{ref_cnt} = $bin_snp{$pos}{ref_cnt};
				}elsif($alt_cnt>1){
					#delete ($bin_homAlts{$pos}) if (exists $bin_homAlts{$pos});
					foreach my $aname (keys %alt_min){
						$bin_hets_ref->{$pos}{alt}{$aname}=$alt_min{$aname};
					} 
					$bin_hets_ref->{$pos}{ref} = $bin_snp{$pos}{ref};
					$bin_hets_ref->{$pos}{ref_cnt} = $bin_snp{$pos}{ref_cnt};
					
				}	
			}		
		}			
	}
}

sub forceMismatchAtHets{
	my $read = shift;
	my $bin_snp_ref = shift;

	my @cols = split /\t/, $read;
	
	my @seq = split //, $cols[9];
	my @quals = split //, $cols[10];
	
	my $start = $cols[3];
	my @cigar = split /(H)/, $cols[5];
	if($cigar[1] eq "H"){ ##Hard clippling at start of read
		$start = $start + $cigar[0];
	}
	
	$read=~m/MD:Z:([ATCG0-9]*)/;
	my $md = $1; 

	my @md_parts = split /([ATCG])/, $md;	
	my $pos = $start;
	my $cnt = 0;
	my $keepRead = 0;

	my %alt_hash = ();
	foreach my $x (@md_parts){
		if($x eq "A" || $x eq "T" || $x eq "C" || $x eq "G"){ ##Different  from reference 
			if(!exists $align_hets{$cols[2]}{$pos}){ ##SNP not known during alignment
				if(exists $bin_hets{$pos}){
					#This signal to keep read will be flagged whenever any of the mismatches are at a predicted heterozygous site
					$keepRead=1;	
				}	
			}
			$pos++;	
			$cnt++;	
		}else{ ##Matches reference
			$pos+=$x;
		}
	}
	if($keepRead==0){
		my $pos = $start;
		my $cnt = 0;
		foreach my $x (@md_parts){
			if($x eq "A" || $x eq "T" || $x eq "C" || $x eq "G"){ ##Different  from reference 
				my $qual = ord($quals[$cnt])-33;
				if ($qual >= $min_qual && $seq[$cnt] ne "N"){
					removeAlt($bin_snp_ref, $pos, $seq[$cnt]);
				}
				$pos++;	
				$cnt++;	
			}else{ ##Matches reference
				for(my $i = 1; $i<=$x; $i++){
					my $qual = ord($quals[$cnt])-33;
					if ($qual >= $min_qual && $seq[$cnt] ne "N"){
						removeRef($bin_snp_ref, $pos);
					}
					$pos++;
					$cnt++;
				}
			}
		}
	}
	return $keepRead;
		
}

sub removeRef{
	my $bin_snp_ref = shift;
	my $pos = shift;
	$bin_snp_ref->{$pos}{ref_cnt}--;

}

sub removeAlt{
	my $bin_snp_ref = shift;
	my $pos = shift;
	my $base = shift;
	$bin_snp_ref->{$pos}{alt}{$base}--;
}

sub createMappabilityReads{
	my $reads_array_ref = shift;
	my @reads_array = @$reads_array_ref;
	my $bin_hets_ref =shift;
	 
	foreach my $read (@reads_array){
		my @cols = split /\t/, $read;
		my @seq = split //, $cols[9];
		my $start = $cols[3];
		my $mdstring = '';
		my $alt_seq = "";
		my $mdcount = 0;
		my $insnp = 0;
		
		my @cigar = split /([MIHDNSHP=X])/, $cols[5];
		my $cur_cnt=0;
		my $cur_pos=$start;
		my $cur_seq_ind = 0;	
		for (my $cigar_ind = 0; $cigar_ind<=$#cigar; $cigar_ind++){
		#if($cigar[1] eq "H"){ ##Hard clippling at start of read
		#	$start = $start + $cigar[0];
		#}
			my $part = $cigar[$cigar_ind];
			if($part eq "M" || $part eq "=" || $part eq "X"){
				my $end = $cur_pos+$cur_cnt;
				for (my $i=$cur_pos; $i<$end; $i++){
					if(exists $bin_hets_ref->{$i}){
						my $base = $seq[$cur_seq_ind];
						my $ref = "";
						my %alt_hash = %{$bin_hets_ref->{$cur_pos}{alt}};
						my @alt_bases = sort {$alt_hash{$b}<=>$alt_hash{$a}} keys %alt_hash;
					
						my $ref_cnt = $bin_hets_ref->{$cur_pos}{ref_cnt} || 0; 
						if ($ref_cnt >= $min_cov){
							$ref = $bin_hets_ref->{$cur_pos}{ref};
						}else{
							$ref = $alt_bases[1];
						}
						if ($base eq $ref){
							$alt_seq .= $alt_bases[0];
						}else{
							$alt_seq .= $ref;
						}
				 		$mdstring .= "$mdcount$base";
						$mdcount = "";
						$insnp = 1;
					}else{
						$alt_seq .= $seq[$cur_seq_ind];
						$mdcount++;
					}
					$cur_seq_ind++;
					$cur_pos++; 			
				}	
			}elsif($part eq "I"){	
				##Do nothing for now	
			}elsif($part eq "D"){	
				##Do nothing for now	
			}elsif($part eq "N"){	
				$cur_pos+=$cur_cnt;	
			}elsif($part eq "S"){	
				$cur_seq_ind += $cur_cnt;	
			}elsif($part eq "H"){
				##Do nothing
			}else{
				$cur_cnt = $part;
			}	
		}
		$mdstring .= "$mdcount" if ($mdcount >0);
		
		if($insnp == 1){
			my $bedStart = $start - 1;	
			my $end = $bedStart + length($cols[9]);	
		 	print READOUT "$cols[2]\t$bedStart\t$end\t$cols[0]\t$cols[9]\n";
                        print ALTREADOUT "$cols[0]" . "_alt\t$bedStart\t$end\t$cols[0]\t$mdstring\t$alt_seq\n";
                        print ALTFAOUT ">$cols[0]"."_alt\n$alt_seq\n";	
		}
	}	
}


sub printBinHets{
	my $het_hash_ref = shift;
	my %het_hash = %{$het_hash_ref};	
	foreach my $pos (sort {$a<=>$b} keys %het_hash){
		my $start = $pos -1;
		my $alleles = "$het_hash{$pos}{ref}:$het_hash{$pos}{ref_cnt}";
		my %alt_hash = %{$het_hash{$pos}{alt}};
		my @keys = sort {$alt_hash{$b}<=>$alt_hash{$a}} keys %alt_hash;
		foreach my $alt_allele (@keys){
			$alleles .= " $alt_allele:$alt_hash{$alt_allele}";
		}		
		my $pvalue = 1;
		if ($het_hash{$pos}{ref_cnt} >= $min_cov){
			$pvalue = findPvalue($alt_hash{$keys[0]},$het_hash{$pos}{ref_cnt});
	        }elsif($#keys >0){
			$pvalue = findPvalue($alt_hash{$keys[0]}, $alt_hash{$keys[1]});
		}
		print HETOUT "$chrom\t$start\t$pos\t$chrom:$pos\t$alleles\t$het_hash{$pos}{ref}\t$pvalue\n";
	}
}

sub findPvalue(){
	my $a1 = shift;
	my $a2 = shift;
	my $total = $a1 + $a2;
	my $p = 1;
	
	if($a1 < $a2){
		$p=2 * Math::CDF::pbinom($a1, $total, .5);
	}elsif($a1 == $a2){
		$p=1;
	}else{
		$p=2 * Math::CDF::pbinom($a2, $total, .5);
	}
	$p = 1 if ($p>1);
	return $p;
}

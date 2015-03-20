#!/usr/perl/bin
use Getopt::Long;
my %opts=(replace=>"", hets=>"");
GetOptions(\%opts, qw(replace=s hets=s output:s chr=s append));

my $file1 = shift;
my $output = $opts{output} || "chrom.fa";
my %snps;
my $sel_chrom = $opts{chr} || "";
warn "Het file excluded, tri-allelic heterozygous variants may not be accuratele represented.\n" if ($opts{hets} eq "");
warn "HomAlt file excluded.\n" if ($opts{replace} eq "");
open(IN1, $opts{replace});
while(<IN1>){
	chomp;
	my @cols= split /\t/, $_;
	if($sel_chrom eq "" || $sel_chrom eq $cols[0]){
	  $snps{$cols[0]}{$cols[2]}{name} = $cols[3];
	  $snps{$cols[0]}{$cols[2]}{chr} = $cols[0];
	  $snps{$cols[0]}{$cols[2]}{allele1} = $cols[4];
	  $snps{$cols[0]}{$cols[2]}{allele2} = $cols[5];
	}
}
close(IN1);

open(IN1, $opts{hets});
while(<IN1>){
	chomp;
	my @cols= split /\t/, $_;
	if($sel_chrom eq "" || $sel_chrom eq $cols[0]){
		$snps{$cols[0]}{$cols[2]}{name} = $cols[3];
		$snps{$cols[0]}{$cols[2]}{chr} = $cols[0];
		$snps{$cols[0]}{$cols[2]}{allele1} = $cols[4];
		$snps{$cols[0]}{$cols[2]}{allele2} = $cols[5];
	}
}
close(IN1);

my $matching;
my $notmatching;
my $total;
if($opts{append} == 0){
	open(OUT, ">$output");
}else{
	open(OUT, ">>$output");
	
}
my $line =1;
my $cnt=1;
my $in_chrom;
my $cur_chrom;
while(<>){
	chomp;
	if(substr($_,0,1) eq ">"){
		my $chr = $_;
		$chr=~s/>//;
		if($sel_chrom eq "" || $sel_chrom eq $chr){
			print OUT "$_\n";
			$in_chrom = 1;
			$cur_chrom = $chr;
			$cnt=1;
		}else{
			$in_chrom = 0;	
		}
	}
	elsif($in_chrom){
		foreach my $base (split //, $_){
			$base = uc($base);
			my $base_comp = $base;
			$base_comp=~tr/ATCG/TAGC/;
			if(exists $snps{$cur_chrom}{$cnt}){
				$total++;
				if($base eq $snps{$cur_chrom}{$cnt}{allele2}){
					print OUT $snps{$cur_chrom}{$cnt}{allele1};
					$notmatching++;
					print "$snps{$cur_chrom}{$cnt}{chr}\t$cnt\t$snps{$cur_chrom}{$cnt}{name}\t$base\t$snps{$cur_chrom}{$cnt}{allele1}\t$snps{$cur_chrom}{$cnt}{allele2}\t1\n"
				}elsif($base eq $snps{$cur_chrom}{$cnt}{allele1}){
					$matching++;
					print "$snps{$cur_chrom}{$cnt}{chr}\t$cnt\t$snps{$cur_chrom}{$cnt}{name}\t$base\t$snps{$cur_chrom}{$cnt}{allele1}\t$snps{$cur_chrom}{$cnt}{allele2}\t0\n";
					print OUT $snps{$cur_chrom}{$cnt}{allele1};
				#}elsif($base_comp eq $snps{$cur_chrom}{$cnt}{allele1}){
				#	$matching_comp++;
				#	print "$snps{$cur_chrom}{$cnt}{chr}\t$cnt\t$snps{$cur_chrom}{$cnt}{name}\t$base\t$snps{$cur_chrom}{$cnt}{allele1}$snps{$cur_chrom}{$cnt}{allele2}\t2\n";
				#	warn "$snps{$cur_chrom}{$cnt}{allele1} $snps{$cur_chrom}{$cnt}{allele2} does not match reference $base at $snps{$cur_chrom}{$cnt}{chr}:$cnt reverse complementing\n";					
				#	$snps{$cur_chrom}{$cnt}{allele1}=~tr/ATCG/TAGC/;
				#	print OUT $snps{$cur_chrom}{$cnt}{allele1};
				#}elsif($base_comp eq $snps{$cur_chrom}{$cnt}{allele2}){
				#	$notmatching_comp++;
				#	print "$snps{$cur_chrom}{$cnt}{chr}\t$cnt\t$snps{$cur_chrom}{$cnt}{name}\t$base\t$snps{$cur_chrom}{$cnt}{allele1}$snps{$cur_chrom}{$cnt}{allele2}\t4\n";
				#	warn "$snps{$cur_chrom}{$cnt}{allele1} $snps{$cur_chrom}{$cnt}{allele2} does not match reference $base at $snps{$cur_chrom}{$cnt}{chr}:$cnt reverse complementing\n";					
				#	$snps{$cur_chrom}{$cnt}{allele1}=~tr/ATCG/TAGC/;
				#	print OUT $snps{$cur_chrom}{$cnt}{allele1};
				
				}else{
					warn "$snps{$cur_chrom}{$cnt}{allele1} $snps{$cur_chrom}{$cnt}{allele2} does not match reference at $base $snps{$cur_chrom}{$cnt}{chr}:$cnt and may be triallelic\n";					
					print "$snps{$cur_chrom}{$cnt}{chr}\t$cnt\t$snps{$cur_chrom}{$cnt}{name}\t$base\t$snps{$cur_chrom}{$cnt}{allele1}\t$snps{$cur_chrom}{$cnt}{allele2}\t1\n";
					if($snps{$cur_chrom}{$cnt}{allele1} ne ""){
						print OUT "$snps{$cur_chrom}{$cnt}{allele1}";
					}else{
						print OUT $base;
					}
				}
			}else{
				print OUT $base;
			}
			$cnt++;
		}
		print OUT "\n";
			
	}
	$line++;	
}
close(OUT);

warn "Matching: $matching Matching_comp:$matching_comp\nNot Matching: $notmatching Not Matching_comp:$notmatching_comp\nTotal: $total\n";

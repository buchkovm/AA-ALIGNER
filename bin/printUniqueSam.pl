#!/usr/bin/env perl
use Getopt::Long;
my %opts=(header=>0, format=>"gsnap", output=>"");
GetOptions(\%opts, qw(header format:s output=s));
my $unmapped;
my $unique;
my $ambiguous;
my $total_reads;
my $hits_ind;
open(OUT, ">$opts{output}.unmapped.sam");
open(OUT1, ">$opts{output}.unique.sam");
open(OUT2, ">$opts{output}.multiple.sam");
open(OUT3, ">$opts{output}.sam.log");

if (lc($opts{format}) eq "bwa"){
	$hits_ind = 13;
}else{
	$hits_ind = 12;
}

while(<>){
	if (substr($_,0,1) eq "@"){
		if($opts{header}){	
			print OUT "$_";
			print OUT1 "$_";
			print OUT2 "$_";
		}
	}else{
		$_=~m/HI\:\w+\:(\d+)/;
		my $hi = $1;
		$hi = 1 if ($opts{format} eq "bwa");
		my @cols = split /\t/, $_;
		my @hits = split /:/, $cols[$hits_ind];
		if(lc($opts{format}) eq "bwa"){
			$_ =~ m/X0\:\w+\:(\d+)/;
			$hits[2] = $1;
		}elsif(lc($opts{format}) eq "gsnap"){
			$_ =~ m/NH\:\w+\:(\d+)/;
			$hits[2] = $1;
		}
		$hits[2] = 0 if ($cols[1] & 0x0004);
		if($hits[2]==0 || $cols[1] & 0x0004){
			$unmapped++;
			$total_reads++;
			print OUT "$_";
		}elsif($hits[2]==1){
			$unique++;
			$total_reads++;
			print OUT1 "$_";
		}elsif($hits[2]>1){
			print OUT2 "$_";
			if($hi ==1){
				$ambiguous++;
				$total_reads++;
			}
		}
	}
}	
print OUT3 "Total Reads: $total_reads\n";
print OUT3 "Unmapped: $unmapped (" . sprintf("%.2f", ($unmapped*100/$total_reads)) . ")\n";
print OUT3 "Uniquely Mapped: $unique (" . sprintf("%.2f", ($unique*100/$total_reads)) . ")\n";
print OUT3 "Ambiguously Mapped: $ambiguous (" . sprintf("%.2f", ($ambiguous*100/$total_reads)) . ")\n";


#!/bin/env perl
use Getopt::Long;
my %opts={filter=>""};
GetOptions(\%opts, qw(output=s filter=s input=s)); 
die "No filter list specified" if ($opts{filter} eq '');

$filter_file = $opts{filter};
$removed_sequence_file = "$opts{output}.removed.txt";
open(LOG, ">$opts{output}.log") || die("Could not write to file!");

print LOG "Starting...\n";

# Generate filter hash table
%filter_hash = ();
open(FLT, $filter_file) || die("Could not open filter file!");
while($record = <FLT>) {
	chomp($record);
	@temp = split('\t', $record);
    for($i = $temp[1]; $i <= $temp[2]; $i++) {
        $filter_hash{$temp[0]}{$i} = 1;
    }
}
close(FLT);

print LOG "Hash table created...\n";

open(OUTF, ">$removed_sequence_file") || die("Could not write to file!");
open(OUT, ">$opts{output}.sam") || die("Could not write to file!");

%other_hash = ();
if($opts{input}=~m/.*\.bam$/){
    open (IN, "samtools view $opts{input} |");
}else{
    open(IN, $opts{input});
} 
while($record = <IN>) {
    chomp($record);
    $overlap = 0;
    my $size = 0;
    if(substr($record, 0, 1) ne "@"){
   	 @temp = split('\t', $record);
         my $stop = $temp[3];
         my @cigars = split /([MIX=SD])/, $temp[5]; 
         my $cur_ind=0;
         for(my $i=0; $i<=$#cigars; $i=$i+2){
        	 if($cigars[$i+1] eq "D"){
              	 	  $stop += $cigars[$i];
       		  }elsif($cigars[$i+1] eq "I"){
         		#Do nothing
		  }else{
             	 	$stop = $stop + $cigars[$i];                         
		  }
	}    

    	$size = length($temp[9]);
   	for($i = $temp[3]; $i <= $stop; $i++) {
         if(exists $filter_hash{$temp[2]}{$i}) {
            $overlap = 1;
       	 }
    	}
    	if($overlap == 0) { #keep it if there is no overlap
       		 print OUT "$record\n";
   	} else {
       		 print OUTF "$record\n";
    	}
    }else{
	print OUT "$record\n";	
    }
}
close(OUTF);
close(OUT);
close(LOG);

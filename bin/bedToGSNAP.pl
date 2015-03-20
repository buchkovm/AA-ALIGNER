my $file = shift;

open (IN, $file);
while(<IN>){
	chomp;
	my @cols = split /\t/, $_;
	print ">$cols[3] $cols[0]:$cols[2] $cols[4]$cols[5]\n";
}
close(IN);

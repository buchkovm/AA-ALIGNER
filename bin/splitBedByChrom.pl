my $file = shift;
my $output_prefix = shift;

my %snps;
open (FILE, $file);
while(<FILE>){
	chomp;
	my @cols =split /\t/, $_;
	$snps{$cols[0]}{$cols[2]}=$_;
}
close(FILE);

foreach my $chrom (keys %snps){
	open (OUT, ">$output_prefix.$chrom.bed");
	my %chr_hash = %{$snps{$chrom}};
	foreach my $pos (sort {$a<=>$b} keys %chr_hash){
		print OUT "$chr_hash{$pos}\n";
	}
	close OUT;
}

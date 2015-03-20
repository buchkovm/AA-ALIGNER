use Getopt::Long;
my %opts = (pvalue=>.01, output=>"Imbalances");
GetOptions(\%opts,qw(output=s filter=s pvalue=s));

my @files = @ARGV;
my %known_hets;
open(IN, $opts{filter});
while(<IN>){
	chomp;
	my @cols = split /\t/, $_;
	$filter{$cols[0]}{$cols[2]}=$cols[3];
}
close (IN);

my $pvalue = $opts{pvalue};
$pvalue =~s/0*\.//;
open(HET, ">$opts{output}.heterozygotes.txt");
open(IMB, ">$opts{output}.heterozygotes.imbalanced.p$pvalue.txt");

if ($opts{filter} ne ""){
	open(KNOWNHET, ">$opts{output}.heterozygotes.included.txt");
	open(KNOWNIMB, ">$opts{output}.heterozygotes.imbalanced.p$pvalue.included.txt");
	open(UNKHET, ">$opts{output}.heterozygotes.predicted.txt");
	open(UNKIMB, ">$opts{output}.heterozygotes.imbalanced.p$pvalue.predicted.txt");
}

foreach my $file (@files){
	open(IN, $file);
	while(<IN>){	
		chomp;
		my @cols = split /\t/, $_;
		my $known = 0;
		$known = 1 if(defined $filter{$cols[0]}{$cols[2]});
		
		print HET "$_\n";
		if ($opts{filter} ne ""){
			if ($known == 1){
				print KNOWNHET "$_\n";
			}else{
				print UNKHET "$_\n";	
			}
		}
		
		if ($cols[6] <= $opts{pvalue}){
			print IMB "$_\n";
			if ($opts{filter} ne ""){
				if ($known == 1){
					print KNOWNIMB "$_\n";
				}else{
					print UNKIMB "$_\n";	
				}
			}
		}
	}
}

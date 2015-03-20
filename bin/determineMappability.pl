my $file = shift;
my $format = shift || "gsnap";

open (IN, "$file.sam");
open (UNI, ">$file.mappability.txt");

while(<IN>){
	chomp;
	my @cols = split /\t/, $_;
	if(substr($_,0,1) ne "@"){
		my $chr = $cols[2];
		my $start = $cols[3];
		my $bedStart = $cols[3]-1;
		my $hits = 0;
		if($format eq "gsnap"){
			#$_=~m/NH:i:(\d+)/;
			$_=~m/NH\:\w+\:(\d+)/;	
			$hits = $1 || 0;
		}else{
			$_=~m/X0\:\w+\:(\d+)/;	
			$hits = $1 || 0;
		}
		my $readlength = length($cols[9]);
		my $bedEnd = $bedStart + $readlength;
		if($chr eq $cols[2] && $start == $cols[3]){
			if($hits ==1){
				print UNI "$chr\t$bedStart\t$bedEnd\t$cols[0]\t$cols[9]\t1\n";
			}else{
				print UNI "$chr\t$bedStart\t$bedEnd\t$cols[0]\t$cols[9]\t$hits\n";
			}
		}else{
			if($hits ==1){
				print UNI "$chr\t$bedStart\t$bedEnd\t$cols[0]\t$cols[9]\t-1\n";
			}else{
				print UNI "$chr\t$bedStart\t$bedEnd\t$cols[0]\t$cols[9]\t-$hits\n";
			}
		
		} 
	}
}
close(IN);

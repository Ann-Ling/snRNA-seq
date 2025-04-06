#!/usr/bin/perl -w
use strict;
open IN,"$ARGV[0]";
my$g=<IN>;
chomp$g;
my@title=split(/\t/,$g);
my%umi;
my%total;
while(<IN>){
	chomp;
	my@A=split(/\t/);
	for(my$i=1;$i<=$#A;$i++){
		$umi{$A[0]}{$title[$i]}=$A[$i];
		$total{$title[$i]}+=$A[$i];
	}
}
close IN;
open IN,"$ARGV[1]";
my%cluster;
my%species;
<IN>;
while(<IN>){
	chomp;
	my@A=split(/\t/);
	push @{$cluster{$A[2]}},$A[0];
	$species{$A[2]}=$A[1];
}
close IN;
open OUT1, ">merged.pheno.cell.txt";
print OUT1 "Sample_ID\tStudy_ID\tCelltype\n";
my%combined_umi;
my%combined_total;
my@new;
foreach my$key(keys %cluster){
	shuffle(\@{$cluster{$key}});
	for(my$i=0;$i<=$#{$cluster{$key}};$i+=10){
		my$id=$key.".".$i;
		print OUT1 "$id\t$species{$key}\t$key\n";
		push @new,$id;
		for(my$j=$i;$j<=$i+9;$j++){
			if($j>$#{$cluster{$key}}){next}
			foreach my$ke(keys %umi){
				$combined_umi{$ke}{$id}+=$umi{$ke}{$cluster{$key}[$j]};
			}
			$combined_total{$id}+=$total{$cluster{$key}[$j]};
		}
	}
}
close OUT1;
open OUT3,">merged.umi.txt";
open OUT2, ">merged.exp.txt";
foreach my$ele(@new){
	print OUT2 "\t$ele";
	print OUT3 "\t$ele";
}
print OUT2 "\n";
print OUT3 "\n";
foreach my$key(keys %combined_umi){
	print OUT2 "$key";
	print OUT3 "$key";
	foreach my$ele(@new){
		my$exp=$combined_umi{$key}{$ele}/$combined_total{$ele}*1000;
		print OUT2 "\t$exp";
		print OUT3 "\t$combined_umi{$key}{$ele}";
	}
	print OUT2 "\n";
	print OUT3 "\n";
}
close OUT2;
close OUT3;
print "done\n";
sub shuffle {
	my $array = shift;
	my $i;
	for ($i = @$array; --$i; ) {
		my $j = int rand ($i+1);
		next if $i == $j;
		@$array[$i,$j] = @$array[$j,$i];
	}
}

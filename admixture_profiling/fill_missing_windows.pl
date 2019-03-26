#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qq(GetOptions);
use Data::Dumper;

#adds missing windows to admixture output
#usage: perl fill_missing_windows.pl --windows 20kb_windows.bed --admixture introgression.tsv --isolates 30 --out f_introgression.tsv

my $windows;
my $admixture;
my $out;
my $isolates;

GetOptions(
	'windows=s' => \$windows,
	'admixture=s' => \$admixture,
	'out=s' => \$out,
	'isolates=i' => \$isolates,
) or die;

my %index;

open (IN1, "$admixture") or die;

while (my $line = <IN1>) {
	chomp $line;
	my @cols = split(/\t/, $line);
	my $cor = "$cols[0]\t$cols[1]\t$cols[2]";
	$index{$cor} = $line; 
}

close IN1;

open (IN2, "$windows") or die;
open (OUT, ">$out") or die;

while (my $window = <IN2>) {
	chomp $window;
	my @win_cols = split(/\t/, $window);
	my $target = "$win_cols[0]\t$win_cols[1]\t$win_cols[2]";
	if (exists $index{$target}) {
		print OUT "$index{$target}\n";
	}
	else {
		my @fill = ("NA") x $isolates;
		print OUT "$target\t";
		print OUT join("\t", @fill), "\n";
	}
}

close IN2;
close OUT;

exit;



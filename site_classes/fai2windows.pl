#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#script to take faidx for a genome, and produce a BED file for a given window size
#usage: perl fai2windows.pl --fai genome.fa.fai --window 100000 --out 100kb_windows.bed

my $fai;
my $window;
my $out;

GetOptions(
        'fai=s' => \$fai,
        'window=s' => \$window,
        'out=s' => \$out,  
) or die "missing input\n";

open (IN, "$fai") or die;
open (OUT, ">$out") or die;

while (my $line = <IN>) {
	chomp $line;
	my @cols = split(/\t/, $line);
	my $length = $cols[1];
	my $sum = 0;
	my $count = 0;
	until ( ($sum + $window) > $length) {
		my $end = $sum + $window;
		$count++;
		print OUT "$cols[0]\t$sum\t$end\t$cols[0].$count\t.\t+\n";
		$sum = $sum + $window;
	}
	$count++;
	print OUT "$cols[0]\t$sum\t$length\t$cols[0].$count\t.\t+\n";
}

close IN;
close OUT;

exit;

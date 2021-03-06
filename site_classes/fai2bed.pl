#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#script to take faidx for a genome, and produce a BED6 file(single line in fai = single line in BED)
#usage: perl fai2bed.pl --fai genome.fa.fai --out chromosomes.bed 

my $fai;
my $out;

GetOptions(
	'fai=s' => \$fai,
        'out=s' => \$out,
) or die "missing input\n";

open (IN, "$fai") or die;
open (OUT, ">$out") or die;

while (my $line = <IN>) {
	chomp $line;
	my @cols = split(/\t/, $line);
	my $length = $cols[1];
	print OUT "$cols[0]\t0\t$length\t$cols[0]\t.\t+\n"
}

close IN;
close OUT;

exit;

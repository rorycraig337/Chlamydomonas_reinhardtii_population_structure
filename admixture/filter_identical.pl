#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qq(GetOptions);

#takes text file of consensus bases for two clades, and outputs lines where both clade have a consensus base that differs
#usage: perl filter_identical.pl --in consensus.txt --out consensus_polymorphic.txt

my $in;
my $out;

GetOptions(
	'in=s' => \$in,
	'out=s' => \$out,
) or die;

open (IN, "$in") or die;
open (OUT, ">$out") or die;

while (my $line = <IN>) {
	chomp $line;
	my @cols = split(/\t/, $line);
	if ( ($cols[1] ne $cols[2]) and ($cols[1] ne "N") and ($cols[2] ne "N") ) {
		print OUT "$line\n";
	}
}

close IN;
close OUT;

exit;

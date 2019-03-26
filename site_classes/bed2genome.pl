#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#converts bed file to genome file
#usage bed2genome.pl

my $bed;
my $genome;

GetOptions(
	'bed=s' => \$bed,
	'genome=s' => \$genome,
) or die "missing input\n";

open (IN, "$bed") or die;
open (OUT, ">$genome") or die;

while (my $line = <IN>) {
	chomp $line;
	my @cols = split(/\t/, $line);
	print OUT "$cols[0]\t$cols[2]\n";
}

close IN;
close OUT;

exit;

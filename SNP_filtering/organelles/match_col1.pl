#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#script filters tsv file based on a exact pattern match in the first column, preserving comment lines
#usage: perl match_col1.pl --tsv in.tsv --match word --out

my $tsv;
my $match;
my $out;

GetOptions(
	'tsv=s' => \$tsv,
	'match=s' => \$match,
	'out=s' => \$out,
) or die "missing input\n";

open (IN, "$tsv") or die;
open (OUT, ">$out") or die;

while (my $line = <IN>) {
	chomp $line;
	unless ($line =~ /^#/) {
		my @cols = split(/\t/, $line);
		if ($cols[0] eq "$match") {
			print OUT "$line\n";
		}
	}
	else {
		print OUT "$line\n";
	}
}

close IN;
close OUT;

exit;

#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#takes list of isolates and output of IBD analysis, and extracts only pairs in the list
#usage: perl extract_pairs.pl --list pairs.txt --IBD IBD.txt --out out.txt

my $pairs;
my $IBD;
my $out;

GetOptions(
        'pairs=s' => \$pairs,
        'IBD=s' => \$IBD,
        'out=s' => \$out,
) or die "missing input\n";

my %index;

open (LIST, "$pairs") or die;

while (my $iso = <LIST>) {
	chomp $iso;
	$index{$iso} = 1; #store isolate IDs
}

close LIST;

open (IN, "$IBD") or die;
open (OUT, ">$out") or die;

while (my $line = <IN>) {
	chomp $line;
	my @cols = split(/\t/, $line);
	if ( (exists $index{$cols[0]}) and (exists $index{$cols[1]}) ) {
		print OUT "$line\n";
	}
}

close IN;
close OUT;

exit;

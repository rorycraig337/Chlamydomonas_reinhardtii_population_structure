#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#takes two lists of isolates and output of IBD analysis, and extracts only pairs between the list
#usage: perl extract_between_pairs.pl --list pairs.txt --IBD IBD.txt --out out.txt

my $pairs1;
my $pairs2;
my $IBD;
my $out;

GetOptions(
        'pairs1=s' => \$pairs1,
        'pairs2=s' => \$pairs2,
        'IBD=s' => \$IBD,
        'out=s' => \$out,
) or die "missing input\n";

my %index1;

open (LIST1, "$pairs1") or die;

while (my $iso1 = <LIST1>) {
	chomp $iso1;
	$index1{$iso1} = 1; #store isolate IDs
}

close LIST1;

my %index2;

open (LIST2, "$pairs2") or die;

while (my $iso2 = <LIST2>) {
	chomp $iso2;
	$index2{$iso2} = 1; #store isolate IDs
}

close LIST2;

open (IN, "$IBD") or die;
open (OUT, ">$out") or die;

while (my $line = <IN>) {
	chomp $line;
	my @cols = split(/\t/, $line);
	if ( ( (exists $index1{$cols[0]}) and (exists $index2{$cols[1]}) ) or ( (exists $index2{$cols[0]}) and (exists $index1{$cols[1]}) ) ) {
		print OUT "$line\n";
	}
}

close IN;
close OUT;

exit;

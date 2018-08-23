#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#script to take output of hmmIBD, and parse IBD tracts of a minimum length ("min" flag)
#additional lengths file is output
#usage: perl extract_IBD.pl --IBD in.IBD.txt --out out.IBD.txt --min 100000 --lengths IBD.lengths.txt

my $IBD;
my $out;
my $min;
my $lengths;

GetOptions(
        'IBD=s' => \$IBD,
        'out=s' => \$out,
	'min=i' => \$min,
	'lengths=s' => \$lengths,
) or die "missing input\n";

open (IN, "$IBD") or die;
open (OUT1, ">$out") or die;
open (OUT2, ">$lengths") or die;

my $header = <IN>;

while (my $line = <IN>) {
	chomp $line;
	my @cols = split(/\t/, $line);
	my $length = $cols[4] - $cols[3];
	print OUT2 "$length\n";
	if ( ($cols[5] == 0) and ($length >= $min) ) {
		print OUT1 "$line\n";
	}
}

close IN;
close OUT1;
close OUT2;

exit;

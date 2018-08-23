#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qq(GetOptions);

#takes a list of isolate pairs and extracts total length of IBD between pairs from IBD tracts file
#usage: perl cumulate_IBD.pl --pairs pairs.txt --IBD tracts.txt --out cumulative_IBD.txt

my $pairs;
my $IBD;
my $out;

GetOptions(
	'pairs=s' => \$pairs,
	'IBD=s' => \$IBD,
	'out=s' => \$out,
) or die "missing input\n";

my %isolate_pairs;

open (IN1, "$pairs") or die;

while (my $line1 = <IN1>) {
	chomp $line1;
	my @cols1 = split(" ", $line1);
	my $pair = "$cols1[0]\t$cols1[1]";
	$isolate_pairs{$pair} = 1;
}

close IN1;

open (IN2, "$IBD") or die;

while (my $line2 = <IN2>) {
	chomp $line2;
	my @cols2 = split(" ", $line2);
	my $match = "$cols2[0]\t$cols2[1]";
	my $length = $cols2[4] - $cols2[3] + 1;
	$isolate_pairs{$match} = $isolate_pairs{$match} + $length;
}

close IN2;

open (OUT, ">$out") or die;

foreach my $key (keys %isolate_pairs) {
	print OUT "$key\t$isolate_pairs{$key}\n";
} 

close OUT;

exit;

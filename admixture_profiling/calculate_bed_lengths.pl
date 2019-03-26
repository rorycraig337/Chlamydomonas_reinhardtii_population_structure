#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#script takes list of isolates and file suffix for bed files, and calcualtes the total number and length of elements in each bed file
#usage: calculate_bed_lengths.pl --isolates all.txt --suffix introgressed.bed --out lengths.txt

my $isolates;
my $suffix;
my $out;

GetOptions(
	'isolates=s' => \$isolates,
	'suffix=s' => \$suffix,
	'out=s' => \$out,
) or die "missing input\n";

open (IN1, "$isolates") or die;
open (OUT, ">$out") or die;

while (my $isolate = <IN1>) {
	chomp $isolate;
	open (IN2, "$isolate.$suffix") or die;
	my $length = 0;
	my $count = 0;
	while (my $line = <IN2>) {
		chomp $line;
		$count++;
		my @cols = split(/\t/, $line);
		$length += $cols[2] - $cols[1];
	}
	close IN2;
	print OUT "$isolate\t$count\t$length\n";
}

close IN1;
close OUT;

exit;

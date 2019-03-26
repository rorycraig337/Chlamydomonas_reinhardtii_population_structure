#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);
use Data::Dumper;

#combines output of IBD summary for different tract lengths, outputting one file for copies and one for tract number
#usage: perl combine_IBD_summary.pl --tract1 1.txt --tract2 2.txt --out_copies copies.tsv --out_tracts tracts.tsv

my $tract1;
my $tract2;
my $out_copies;
my $out_tracts;

GetOptions(
	'tract1=s' => \$tract1,
	'tract2=s' => \$tract2,
	'out_copies=s' => \$out_copies,
	'out_tracts=s' => \$out_tracts,
) or die "missing input\n";

my %index1;
my %index2;

open (IN1, "$tract1") or die;

while (my $line1 = <IN1>) {
	chomp $line1;
	my @cols1 = split(/\t/, $line1);
	$index1{$cols1[3]} = $cols1[6];
	$index2{$cols1[3]} = $cols1[7];
}

close IN1;

open (IN2, "$tract2") or die;;
open (OUT1, ">$out_tracts") or die;
open (OUT2, ">$out_copies") or die;

while (my $line2 = <IN2>) {
	chomp $line2;
	my @cols2 = split(/\t/, $line2);
	print OUT1 "$cols2[3]\t$index1{$cols2[3]}\t$cols2[6]\n";
	print OUT2 "$cols2[3]\t$index2{$cols2[3]}\t$cols2[7]\n";
}

close IN2;
close OUT1;
close OUT2;

exit;

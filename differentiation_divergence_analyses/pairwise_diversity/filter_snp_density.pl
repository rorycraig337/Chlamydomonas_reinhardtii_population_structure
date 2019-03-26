#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qq(GetOptions);

#takes file with SNP density per window, and filters windows in a pairwise comparison based min_SNPs flag
#usage: perl filter_snp_density.pl --density density.tsv --min_SNP 10 --in in.tsv --out out.tsv 

my $density;
my $in;
my $min_SNPs;
my $out;

GetOptions(
	'density=s' => \$density,
	'min_SNPs=i' => \$min_SNPs,
	'in=s' => \$in,
	'out=s' => \$out,
) or die;

my %low_SNP;

open (IN1, "$density") or die;

while (my $line1 = <IN1>) {
	chomp $line1;
	my @cols1 = split(/\t/, $line1);
	my $win = "$cols1[0]\t$cols1[1]\t$cols1[2]";
	if ($cols1[3] < $min_SNPs) {
		$low_SNP{$win} = 1;
	}
}

close IN1;

open (IN2, "$in") or die;
open (OUT, ">$out") or die;

while (my $con = <IN2>) {
	chomp $con;
	my @cols = split(" ", $con);
	my $win_match = "$cols[0]\t$cols[1]\t$cols[2]";
	unless (exists $low_SNP{$win_match}) {
		print OUT "$con\n";
	}
}

close IN2;
close OUT;

exit;

#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qq(GetOptions);

#script to take VCF and table of SNPs build from VCF with positions based on SNPs (1, 2, 3, 4 etc.) and add in correct site information
#usage: perl add_genomic_positions.pl --vcf in.vcf --snps in.tsv --out out.tsv

my $vcf;
my $snps;
my $out;

GetOptions(
	'vcf=s' => \$vcf,
	'snps=s' => \$snps,
	'out=s' => \$out,
) or die;

open (IN1, "$vcf") or die;

my %positions;
my $count = 0;

while (my $line1 = <IN1>) {
	chomp $line1;
	unless ($line1 =~ /^#/) {
		$count++;
		my @cols1 = split(/\t/, $line1);
		$positions{$count} = "$cols1[0]\t$cols1[1]";
	}
}

close IN1;

open (IN2, "$snps") or die;
open (OUT, ">$out") or die;

while (my $line2 = <IN2>) {
	chomp $line2;
	my @cols2 = split(/\t/, $line2);
	my $pos = $positions{$cols2[0]};
	print OUT "$pos\t$cols2[1]\t$cols2[2]\n";
}

close IN2;
close OUT;

exit;

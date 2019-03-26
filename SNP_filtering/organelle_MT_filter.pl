#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qq(GetOptions);

#filters cpDNA, mtDNA and mating type locus: --mating should supply single line bed file with mating type locus
#usage: perl organelle_MT_filter.pl --vcf in.vcf --mating MT.bed --out out.vcf

my $vcf;
my $mating;
my $out;

GetOptions(
	'vcf=s' => \$vcf,
	'mating=s' => \$mating,
	'out=s' => \$out,
) or die "missing input\n";

open (MT, "$mating") or die;

my $mt_locus = <MT>;
my @mt_cols = split(/\t/, $mt_locus);

open (IN, "$vcf") or die;
open (OUT, ">$out") or die;

while (my $line = <IN>) {
	chomp $line;
	my @cols = split(/\t/, $line);
	unless ( ($line =~ /^cpDNA/) or ($line =~ /^mtDNA/) or ($line =~ /^mtMinus/) or ( ($cols[0] eq $mt_cols[0]) and ($cols[1] > $mt_cols[1]) and ($cols[1] <= $mt_cols[2]) ) ) {
		print OUT "$line\n";
	}
}

close IN;
close OUT;

exit;

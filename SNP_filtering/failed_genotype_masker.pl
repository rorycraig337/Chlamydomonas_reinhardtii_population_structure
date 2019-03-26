#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qq(GetOptions);

#script masks genotypes for individuals that fail a GQ and DP filter, inserting word "null" in place of genotype info
#usage: perl failed_genotype_masker.pl --vcf in.vcf --GQ 20 --DP_min 3 --DP_max 60 --out out.vcf

my $vcf;
my $GQ;
my $DP_min;
my $DP_max;
my $out;

GetOptions(
    'vcf=s' => \$vcf,
    'GQ=i' => \$GQ,
    'DP_min=i' => \$DP_min,
    'DP_max=i' => \$DP_max,
    'out=s' => \$out,
) or die "missing input\n";

open (IN, "$vcf") or die;
open (OUT, ">$out") or die;

while (my $line = <IN>) {
	chomp $line;
	unless ($line =~ /^#/) {
		my @cols = split(/\t/, $line);
		my $count = 0;
		foreach my $col (@cols) {
			if ($count > 8) {
				my @calls = split(/:/, $col);
				unless ( ($calls[3] >= $GQ) and ($calls[2] >= $DP_min) and ($calls[2] <= $DP_max) ) {
					$col = "null";
				}
			}
			$count++;
		}
		print OUT join("\t", @cols), "\n";
	}
	else {
		print OUT "$line\n";
	}
}

close IN;
close OUT;

exit;

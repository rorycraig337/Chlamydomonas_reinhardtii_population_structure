#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qq(GetOptions);

#takes a VCF of SNPs and window file and outputs number of SNPs per window
#usage: perl tally_window_snps.pl --windows in.bed --vcf in.vcf --out out.tsv

my $windows;
my $vcf;
my $out;

GetOptions(
	'windows=s' => \$windows,
	'vcf=s' => \$vcf,
	'out=s' => \$out,
) or die;

my %window_hash;

open (IN1, "$windows") or die;

while (my $window = <IN1>) {
	chomp $window;
	my @win_cols = split(/\t/, $window);
	my $chr = "$win_cols[0]";
	my $bait = "$win_cols[1]\t$win_cols[2]";
	$window_hash{$chr}{$bait} = 0;
}

close IN1;

open (IN2, "$vcf") or die;

while (my $line = <IN2>) {
	chomp $line;
	unless ($line =~ /^#/) {
		my @cols = split(/\t/, $line);
		my $chrom = $cols[0];
		my $site = $cols[1];
		foreach my $win (keys %{ $window_hash{$chrom} }) {
			my @win_info = split(/\t/, $win);
			if ( ($site > $win_info[0]) and ($site <= $win_info[1]) ) {
				$window_hash{$chrom}{$win} = $window_hash{$chrom}{$win} + 1;
			}
		}
	unless ($site % 10000) {
		print "$site\n";
	}
	}
}

close IN2;

open (OUT, ">$out") or die;

foreach my $chromosome (sort keys %window_hash) {
	foreach my $window (keys %{ $window_hash{$chromosome} }) {
		print OUT "$chromosome\t$window\t$window_hash{$chromosome}{$window}\n";
	}
}

close OUT;

exit;

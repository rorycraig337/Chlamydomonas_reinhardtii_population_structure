#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);
use Data::Dumper;

#checks number of SNPs occuring in each IBD tract
#usage: perl IBD_check.pl --IBD ibd.tsv --snps hmmIBD.txt --out IBD_snps.txt

my $IBD;
my $snps;
my $out;

GetOptions(
	'IBD=s' => \$IBD,
	'snps=s' => \$snps,
	'out=s' => \$out,
) or die "missing input\n"; 

my %index;
my %IDs;

open (IN1, "$snps") or die;

print "ready to read...\n";

my $header = <IN1>;
chomp $header;
my @hcols = split (/\t/, $header);
my $hcol_count = 0;
foreach my $hcol (@hcols) {
	$hcol_count++;
	if ($hcol_count > 2) {
		$IDs{$hcol_count} = $hcol;
	}
}

my $chr = 1;

while (my $line = <IN1>) {
	chomp $line;
	my @cols = split (/\t/, $line);
	if ($cols[0] ne "$chr") {
		print "$chr\n";
	}
	$chr = $cols[0]; 
	my $col_count = 0;
	foreach my $col (@cols) {
		$col_count++;
		if ($col_count > 2) {
			my $iso = $IDs{$col_count};
			$index{$cols[0]}{$cols[1]}{$iso} = $col;
		}
	}
}

close IN1;

print "nice!\nnow the big check...\n";

open (IN2, "$IBD") or die;
open (OUT, ">$out") or die;

while (my $tracts = <IN2>) {
	chomp $tracts;
	my @tcols = split (/\t/, $tracts);
	my $snp_count = 0;
	for (my $i=$tcols[3]; $i <= $tcols[4]; $i++) {
		if (exists $index{$tcols[2]}{$i}) {
			if ( ( ($index{$tcols[2]}{$i}{$tcols[0]} == 0) or ($index{$tcols[2]}{$i}{$tcols[0]} == 1) ) and ( ($index{$tcols[2]}{$i}{$tcols[1]} == 0) or ($index{$tcols[2]}{$i}{$tcols[1]} == 1) ) ) {
				if ($index{$tcols[2]}{$i}{$tcols[0]} != $index{$tcols[2]}{$i}{$tcols[1]}) {
					$snp_count++;
				}
			}
		}
	}
	print OUT "$tracts\t$snp_count\n";
}

close IN2;
close OUT;

exit;

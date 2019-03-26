#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);
use Data::Dumper;

#script to pre-filter VCF prior to popgen statistics, removing sites that fail for minumum good calls for one or both clades
#neccessary as otherwise popgen will be done on sites that pass for one clade, but fail for the other
#assumes filtering on quality/depth has already been performed i.e. takes presence of 1/0 genotype as a good call

my $vcf;
my $pop_list;
my $minA;
my $minB;
my $out;

GetOptions(
	'vcf=s' => \$vcf,
	'pop_list=s' => \$pop_list,
	'minA=i' => \$minA,
	'minB=i' => \$minB,
	'out=s' => \$out,
) or die "missing input\n";

open (POP, "$pop_list") or die "pop list file name invalid\n";

my %cladeA;
my %cladeB;

my $pop1_inds = <POP>;
my $pop2_inds = <POP>;
chomp $pop1_inds;
chomp $pop2_inds;
my @pop1_list = split(/,/, $pop1_inds);
my @pop2_list = split(/,/, $pop2_inds);
foreach my $ind1 (@pop1_list) {
	$cladeA{$ind1} = 1;
}
foreach my $ind2 (@pop2_list) {
	$cladeB{$ind2} = 1;
}

close POP;	

open (IN, "$vcf") or die;
open (OUT, ">$out") or die;

my %clade_index;

my $A_total;
my $B_total;
my $sites;

while (my $line = <IN>) {
	chomp $line;
	if ($line =~ /^#CHROM/) {
		my @head_cols = split(/\t/, $line);
		my $head_col_count = 0;
		foreach my $head_col (@head_cols) {
			if ($head_col_count > 8) {
				if (exists $cladeA{$head_col}) {
					$clade_index{A}{$head_col_count} = 1;
				}
				elsif (exists $cladeB{$head_col}) {
					$clade_index{B}{$head_col_count} = 1;
				}
			}
			$head_col_count++;
		}
		print OUT "$line\n";
	}
	elsif ($line !~ /^#/) {
		my @vcf_cols = split(/\t/, $line);
		my $vcf_count = 0;
		my $A_count = 0;
		my $B_count = 0;
		foreach my $vcf_col (@vcf_cols) {
			if ( ($vcf_count > 8) and (exists $clade_index{A}{$vcf_count}) ) {
				my @calls = split(/:/, $vcf_col);
				if ( ($calls[0] eq "0") or ($calls[0] eq "1") ) {
					$A_count++;
				}
			}
			if ( ($vcf_count > 8) and (exists $clade_index{B}{$vcf_count}) ) {
				my @calls = split(/:/, $vcf_col);
				if ( ($calls[0] eq "0") or ($calls[0] eq "1") ) {
					$B_count++;
				}
			}
			$vcf_count++;
		}
		if ( ($A_count >= $minA) and ($B_count >= $minB) ) {
			$A_total += $A_count;
			$B_total += $B_count;
			$sites++;
			print OUT "$line\n";
		}
	}
	else {
		print OUT "$line\n";
	}

}

my $A_mean = $A_total / $sites;
my $B_mean = $B_total / $sites;

print "$sites passed, mean A isolates = $A_mean, mean B isolates = $B_mean\n";

close IN;
close OUT;

exit;



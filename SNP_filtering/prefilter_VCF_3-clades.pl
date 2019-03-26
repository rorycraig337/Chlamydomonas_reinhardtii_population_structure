#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);
use Data::Dumper;

#script to pre-filter VCF prior to popgen statistics, removing sites that fail for minumum good calls for one or both clades
#neccessary as otherwise popgen will be done on sites that pass for one clade, but fail for the other
#assumes filtering on quality/depth has already been performed i.e. takes presence of 1/0 genotype as a good call
#pop_list takes three line file, with comma separated list of isolates in each clade on each line
#usage: perl prefilter_VCF.pl --vcf all_sites.vcf --pop_list clades.txt --minA 10 --minB 4 --minC 2 --out snps.min_isolates.vcf

my $vcf;
my $pop_list;
my $minA;
my $minB;
my $minC;
my $out;

GetOptions(
	'vcf=s' => \$vcf,
	'pop_list=s' => \$pop_list,
	'minA=i' => \$minA,
	'minB=i' => \$minB,
	'minC=i' => \$minC,
	'out=s' => \$out,
) or die "missing input\n";

#collect clade info

open (POP, "$pop_list") or die "pop list file name invalid\n";

my %cladeA;
my %cladeB;
my %cladeC;

my $pop1_inds = <POP>;
my $pop2_inds = <POP>;
my $pop3_inds = <POP>;
chomp $pop1_inds;
chomp $pop2_inds;
chomp $pop3_inds;
my @pop1_list = split(/,/, $pop1_inds);
my @pop2_list = split(/,/, $pop2_inds);
my @pop3_list = split(/,/, $pop3_inds);
foreach my $ind1 (@pop1_list) {
	$cladeA{$ind1} = 1;
}
foreach my $ind2 (@pop2_list) {
	$cladeB{$ind2} = 1;
}
foreach my $ind3 (@pop3_list) {
        $cladeC{$ind3} = 1;
}
close POP;	

#loop through VCF, and filter sites with too few isolates for either clade

open (IN, "$vcf") or die;
open (OUT, ">$out") or die;

my %clade_index;

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
				if (exists $cladeB{$head_col}) {
					$clade_index{B}{$head_col_count} = 1;
				}
                                elsif (exists $cladeC{$head_col}) {
                                        $clade_index{C}{$head_col_count} = 1;
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
		my $C_count = 0;
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
                        if ( ($vcf_count > 8) and (exists $clade_index{C}{$vcf_count}) ) {
                                my @calls = split(/:/, $vcf_col);
                                if ( ($calls[0] eq "0") or ($calls[0] eq "1") ) {
                                        $C_count++;
                                }
                        }
			$vcf_count++;
		}
		if ( ($A_count >= $minA) and ($B_count >= $minB) and ($C_count >= $minC) ) {
			print OUT "$line\n";
		}
	}
	else {
		print OUT "$line\n";
	}

}

close IN;
close OUT;

exit;


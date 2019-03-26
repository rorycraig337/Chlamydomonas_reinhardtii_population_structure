#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qq(GetOptions);
use Data::Dumper;

#script masks genotypes for individuals that fail a GQ and DP filter, inserting word "null" in place of genotype info
#GQ is only filtered for SNPs, invariant sites are only filtered on depth
#DP_max is a two column tsv with isolate in column one, and mean converage in column two
#DP_max_int is used to change the max depth, which will be set as mean coverage + (DP_max_int * sqrt(mean coverage))
#usage: perl failed_genotype_masker_maxDP.pl --vcf in.vcf --GQ 20 --DP_min 3 --DP_max cov.tsv --DP_max_int 3 --out out.vcf

my $vcf;
my $GQ;
my $DP_min;
my $DP_max;
my $DP_max_int;
my $out;

GetOptions(
    'vcf=s' => \$vcf,
    'GQ=i' => \$GQ,
    'DP_min=i' => \$DP_min,
    'DP_max=s' => \$DP_max,
    'DP_max_int=i' => \$DP_max_int,
    'out=s' => \$out,
) or die "missing input\n";

my %max;

open (TSV, "$DP_max") or die;

while (my $tsv_line = <TSV>) {
	chomp $tsv_line;
	my @tsv_cols = split(/\t/, $tsv_line);
	my $cutoff = $tsv_cols[1] + ($DP_max_int * sqrt($tsv_cols[1]));
	$max{$tsv_cols[0]} = $cutoff; #store hash of isolate names and per isolate max DP
}

close TSV;

open (IN, "$vcf") or die;
open (OUT, ">$out") or die;

my %isolates;
my $ind_index_check = 0;

while (my $line = <IN>) {
	chomp $line;
	if ($line =~ /^#CHROM/) {
		my @head_cols = split(/\t/, $line);
		my $head_col_count = 0;
		foreach my $head_col (@head_cols) {
			if ($head_col_count > 8) {
				if (exists $max{$head_col}) {
					$isolates{$head_col_count} = $max{$head_col}; #store hash of per isolate column number and per isolate max DP
				}
				else {
					die "isolate present in VCF is missing in DP_max file\n";
				}
			}
			$head_col_count++;
		}
		$ind_index_check = 1;
	}
	unless ($line =~ /^#/) {
		if ($ind_index_check == 0) {
			die "VCF file does not contain header with valid information\n";
		}
		my @cols = split(/\t/, $line);
		my $SNP_flag = 0;
		my $count = 0;
		if ($cols[8] !~ "DP") { #some rare sites have 0 depth for all isolates, can be completely excluded
			next;
		}
		if ($cols[4] ne ".") { #site is a SNP
			$SNP_flag = 1;
		}
		foreach my $col (@cols) {
			if ($count > 8) {
				my @calls = split(/:/, $col);
				my $isolate_DP_max = $isolates{$count}; #for each isolate retrieve the associated max DP
				if ($SNP_flag == 1) {
					if ( ($calls[0] eq ".") or ($calls[0] eq "null") ) {
						$col = "null"; #site has no call, mask 
					}
					else { #site has a genotype
						unless ( ($calls[2] ne ".") and ($calls[3] >= $GQ) and ($calls[2] >= $DP_min) and ($calls[2] <= $isolate_DP_max) ) {
							$col = "null"; #mask to 'null' if genotype fails any filters 
						}
						else { #passed all filters, now check if >=90% of reads support the called allele 
							my $ad_count = 0;
							my @ad_vals = split(/,/, $calls[1]);
							my $ad_total = 0;
							my $ad_snp = 0;
							foreach my $ad (@ad_vals) {
								$ad_total += $ad;
								if ( ( ($calls[0] eq "0") and ($ad_count == 0) ) or ( ($calls[0] eq "1") and ($ad_count == 1) ) ) { #relevant AD for the allele
									$ad_snp = $ad;
								}
								$ad_count++;
							}
							unless ($ad_count > 2) { #ignore these sites, were dealt with earlier in pipeline
								if ($ad_snp == 0) {
									$col = "null"; #mask to 'null', mix of reads at this site
								}
								else {
									my $ad_percent = ($ad_snp / $ad_total) * 100;
									if ($ad_percent < 90) {
										$col = "null"; #mask to 'null', mix of reads at this site
									}
								}
							}
						}
					}
				}
				else { #site is invariant, only filter on depth here
					if ( ($calls[0] eq ".") or ($calls[0] eq "null") ) {
						$col = "null"; #site has no call, mask 
					}
					else {
						unless ( ($calls[2] ne ".") and ($calls[2] >= $DP_min) and ($calls[2] <= $isolate_DP_max) ) {
							$col = "null"; #mask to 'null' if genotype fails any filters 
						}
					}
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



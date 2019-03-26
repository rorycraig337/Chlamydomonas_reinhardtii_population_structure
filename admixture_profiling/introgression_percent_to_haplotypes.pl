#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#script to take file of windows with percentage of introgressed SNPs per window per isolate, and join windows into introgressed haplotypes
#haplotype must start and end with window with >= 'flank' percent introgressed SNPs, and can contain windows with >= 'inner' percent introgressed SNPs
#can also allow haplotypes to pass over NA blocks, up to a maximum of 'NA' (e.g. so centromere is exlcuded)
#'break' sets the number of < 'inner' windows needed before haplotype is broken (for sliding windows a single window may be too conservative)
#assumes input file is sorted by chromosome and start coordinate of window
#usage: perl introgression_percent_to_haplotypes.pl --admixture snps_percentages.tsv --isolates names.txt --flank 60 --inner 40 --NA 12 --break 5 --suffix introgressed.bed

my $admixture;
my $isolates;
my $flank;
my $inner;
my $NA;
my $break;
my $suffix;

GetOptions(
	'admixture=s' => \$admixture,
	'isolates=s' => \$isolates,
	'flank=i' => \$flank,
	'inner=i' => \$inner,
	'NA=i' => \$NA,
	'break=i' => \$break,
	'suffix=s' => \$suffix,
) or die "missing input\n";

my %isolate_index;

open (IN1, "$isolates") or die;

my $counter = 0;

while (my $isolate = <IN1>) { #loop through isolate list and store in hash
	chomp $isolate;
	$counter++;
	$isolate_index{$counter} = $isolate;
}

close IN1;

for (my $i=1; $i <= $counter; $i++) { #foreach isolate
	open (IN2, "$admixture") or die;
	my $id = $isolate_index{$i}; #get isolate name
	open (OUT, ">$id.$suffix") or die;
	my $chrom = 0;
	my $hap = 0;
	my $start = 0;
	my $end = 0;
	my $NA_count = 0;
	my $low_count = 0;
	while (my $line = <IN2>) { #loop through admixture file
		chomp $line;
		my @cols = split(/\t/, $line);
		my $window = "$cols[0]\t$cols[1]\t$cols[2]";
		if ($cols[0] ne "$chrom") { #check if a new chromosome, if so break haplotype
			if ($hap == 1) {
				print OUT "$chrom\t$start\t$end\n";
			}
			$hap = 0;
			$start = 0;
			$end = 0;
			$NA = 0;
			$low_count = 0;
		}
		$chrom = $cols[0];
		my $col_count = 0;
		foreach my $col (@cols) {
			$col_count++;
			if ($col_count == ($i + 3)) { #correct column for isolate
				if ($col eq "NA") { 
					if ($hap == 1) { #currently in haplotype, check number of preceeding NAs
						$NA_count++;
						if ($NA_count > $NA) { #if max NA reached, break haplotype
							print OUT "$chrom\t$start\t$end\n";
							$hap = 0;
							$start = 0;
							$end = 0;
							$NA_count = 0;	
						}
					}
				}
				elsif ( ($hap == 1) and ($col > $inner) ) { #currently in haplotype, and window passes
					if ($col >= $flank) { #could be end of haplotype, store coordinate
						$end = $cols[2];
					}
					$NA_count = 0;
					$low_count = 0;
				}
				elsif ( ($hap == 0) and ($col > $flank) ) { #start of new haplotype
					$start = $cols[1];
					$end = $cols[2]; #could be single window haplotype so store end too
					$hap = 1;	
				}
				elsif ( ($hap == 1) and ($col <= $inner) ) { #potentially end of haplotype
					$low_count++;
					if ($low_count > $break) { #end of haplotype
						print OUT "$chrom\t$start\t$end\n";
						$hap = 0;
						$start = 0;
						$end = 0;
						$NA_count = 0;
					}	
				}
			}
		}
	}
	close IN2;
	close OUT;
}

exit;

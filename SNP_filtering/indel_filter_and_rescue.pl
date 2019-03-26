#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);
use Data::Dumper;

#script creates per-isolate mask bed files for x bp either side of an indel, set with flag 'flank'
#requires all sites VCF, outputs per isolate bed file of coordinates to be masked and a new VCF
#mode must be set to discard or rescue, discard removes any site with an indel
#rescue will change indel sites that are otherwise invariant/biallelic for SNPs to that state
#usage: perl indel_filter_and_rescue.pl --vcf in.vcf --out indel_processed.vcf --flank 5 --suffix indel_mask.bed --mode rescue

my $vcf;
my $flank;
my $suffix;
my $out;
my $mode;

GetOptions(
	'vcf=s' => \$vcf,
	'flank=s' => \$flank,
	'out=s' => \$out,
	'suffix=s' => \$suffix,
	'mode=s' => \$mode,
) or die "missing input\n";

open (IN, "$vcf") or die;
open (OUT, ">$out") or die;

my %isolate_index;
my %indel_index;
my $ind_index_check = 0;

while (my $line = <IN>) {
	chomp $line;
	my @cols = split(/\t/, $line);
	if ($line =~ /^#CHROM/) {
		print OUT "$line\n";
		my $head_col_count = 0;
		foreach my $head_col (@cols) {
			if ($head_col_count > 8) {
				$isolate_index{$head_col_count} = "$head_col"; #create hash linking column number to isolate ID
			}
			$head_col_count++;
		}
		$ind_index_check = 1;
	}
	elsif ($line !~ /^#/) {
		if ( (length $cols[3] == 1) and (length $cols[4] > 1) and ($cols[4] !~ /,/) ) { #biallelic insertion
			my $start = $cols[1] - $flank; #bed is 0-based, but only need to minus flank as coordinate base is same between ref and indel
			my $end = $cols[1] + $flank;
			my $col_count = 0;
			foreach my $col (@cols) {
				if ($col_count > 8) {
					my @calls = split(/:/, $col);
					if ( ($calls[0] ne ".") and ($calls[0] == 1) ) { #indel is present
						my $ID = $isolate_index{$col_count};
						my $mask = "$cols[0]\t$start\t$end";
						push(@{$indel_index{$ID}}, $mask);
					}
				}
				$col_count++;
			}
			if ($mode eq "rescue") {
				$cols[4] = ".";
				print OUT join("\t", @cols), "\n";
			}
			elsif ($mode ne "discard") {
				die "mode not set correctly, discard or rescue?\n";
			}
		}
		elsif ( (length $cols[3] > 1) and (length $cols[4] == 1) and ($cols[4] !~ /,/) ) { #biallelic deletion
			if ($cols[4] ne ".") {
				my $start;
				if ($cols[4] eq '*') {
					$start = $cols[1] - $flank - 1; #minus 1 as base is deleted itself 
				}
				else {
					$start = $cols[1] - $flank;
				} 
				my $end = $cols[1] + (length $cols[3]) + $flank - 1;
				my $col_count = 0;
				foreach my $col (@cols) {
					if ($col_count > 8) {
						my @calls = split(/:/, $col);
						if ( ($calls[0] ne ".") and ($calls[0] == 1) ) { #indel is present
							my $ID = $isolate_index{$col_count};
							my $mask = "$cols[0]\t$start\t$end";
							push(@{$indel_index{$ID}}, $mask);
						}
					}
					$col_count++;
				}
				if ($mode eq "rescue") {
					my $base = substr($cols[3], 0, 1);
					$cols[3] = $base;
					if ( ($cols[4] eq "$base") or ($cols[4] eq '*') ) {
						$cols[4] = ".";
						print OUT join("\t", @cols), "\n";
					}
					else {
						print "WARNING: multi-allelic? $line\n";
					}
				}
			}
			else { #was an indel in the raw version of the VCF, not in this reduced version
				if ($mode eq "rescue") {
					my $base = substr($cols[3], 0, 1);
					$cols[3] = $base;
					print OUT join("\t", @cols), "\n";
				}
				my $col_count = 0;
				foreach my $col (@cols) {
					if ($col_count > 8) {
						my @calls = split(/:/, $col);
						if ( ($calls[0] ne ".") and ($calls[0] != 0) ) { #still variation at this site, there shouldn't be
							die "error: $line\n";
						}
					}
					$col_count++;
				}
			}
		}
		elsif ( (length $cols[3] > 1) and (length $cols[4] > 1) and ($cols[4] !~ /,/) ) {
			die "does this happen? $line\n";
		}
		elsif ( (length $cols[3] == 1) and ($cols[4] eq '*') ) { #biallelic single bp deletion
			my $start = $cols[1] - $flank - 1;
			my $end = $cols[1] + $flank;
			my $col_count = 0;
			foreach my $col (@cols) {
				if ($col_count > 8) {
					my @calls = split(/:/, $col);
					if ( ($calls[0] ne ".") and ($calls[0] == 1) ) { #indel is present
						my $ID = $isolate_index{$col_count};
						my $mask = "$cols[0]\t$start\t$end";
						push(@{$indel_index{$ID}}, $mask);
					}
				}
				$col_count++;
			}
			if ($mode eq "rescue") {
				$cols[4] = ".";
				print OUT join("\t", @cols), "\n";
			}
		}
		elsif ( (length $cols[3] == 1) and ($cols[4] =~ /,/) ) { #multi-allelic, may be indel
			my @alleles;
			push @alleles, $cols[3];
			my @allele_cols = split(/,/, $cols[4]);
			foreach my $allele1 (@allele_cols) {
				push @alleles, $allele1;
			}
			my $allele_count = 0;
			my $SNP_count = 0;
			my $SNP;
			my $SNP_mark;
			foreach my $allele2 (@alleles) {
				if ($allele_count > 0) {
					if (length $allele2 > 1) { #allele is an insertion
						my $start = $cols[1] - $flank;
						my $end = $cols[1] + $flank;
						my $col_count = 0;
						foreach my $col (@cols) {
							if ($col_count > 8) {
								my @calls = split(/:/, $col);
								if ( ($calls[0] ne ".") and ($calls[0] == $allele_count) ) { #indel is present
									my $ID = $isolate_index{$col_count};
									my $mask = "$cols[0]\t$start\t$end";
									push(@{$indel_index{$ID}}, $mask);
								}
							}
							$col_count++;
						}
					}
					elsif ($allele2 eq '*') { #allele is a single bp deletion
						my $start = $cols[1] - $flank - 1;
						my $end = $cols[1] + $flank;
						my $col_count = 0;
						foreach my $col (@cols) {
							if ($col_count > 8) {
								my @calls = split(/:/, $col);
								if ( ($calls[0] ne ".") and ($calls[0] == $allele_count) ) { #indel is present
									my $ID = $isolate_index{$col_count};
									my $mask = "$cols[0]\t$start\t$end";
									push(@{$indel_index{$ID}}, $mask);
								}
							}
							$col_count++;
						}
					}
					elsif ($allele2 ne ".") { #allele is single base
						$SNP_count++;
						$SNP = $allele2;
						$SNP_mark = $allele_count;
					}
				}
				$allele_count++;
			}
			if ($mode eq "rescue") {
				if ($SNP_count == 0) { #restore site as invariant
					$cols[4] = ".";
					print OUT join("\t", @cols), "\n";	
				}
				elsif ($SNP_count == 1) { #restore site as SNP
					$cols[4] = $SNP;
					my $col_count = 0;
					foreach my $col (@cols) {
						if ($col_count > 8) {
							my @calls = split(/:/, $col);
							if ( ($calls[0] ne ".") and ( ($calls[0] == $SNP_mark) or ($calls[0] eq "0") ) ) { #SNP needs to now take alt value of 1, invariant must be checked for AD
								my $ad_count = 0;
								my @ad_vals = split(/,/, $calls[1]);
								my $ad_total = 0;
								my $ad_snp = 0;
								foreach my $ad (@ad_vals) {
									$ad_total += $ad;
									if ( ( ($calls[0] eq "0") and ($ad_count == 0) ) or ( ($calls[0] ne "0") and ($ad_count == $SNP_mark) ) ) { #relevant AD for the allele
										$ad_snp = $ad;
									}
									$ad_count++;
								}
								if ($ad_snp == 0) {
									print OUT "\tnull"; #mask this site for this isolate, no informative reads
								}
								else {
									my $ad_percent = ($ad_snp / $ad_total) * 100;
									if ($ad_percent >= 90) { #allele is at frequency >= 90%, retain
										if ($calls[0] eq "0") { #site is invariant
											print OUT "\t";
											print OUT join(":", @calls);											
										}
										else { #site is a SNP, change allele to 1
											$calls[0] = 1;
											print OUT "\t";
											print OUT join(":", @calls);
										}
									}
									else {
										print OUT "\tnull"; #mask this site for this isolate
									}
								}	
							}
							else { #site is indel, will be masked later
								print OUT "\t$col";
							}
						}
						else {
							if ($col_count == 0) {
								print OUT "$col";
							}
							else {
								print OUT "\t$col";
							}
						}
						$col_count++;
					}
					print OUT "\n";
				}
			}
		}
		elsif ( (length $cols[3] > 1) and ($cols[4] =~ /,/) ) { #multi-allelic, is an indel
			my $print_count = 0;
			my @alleles;
			push @alleles, $cols[3];
			my $ref_length = length $cols[3];
			my @allele_cols = split(/,/, $cols[4]);
			foreach my $allele1 (@allele_cols) {
				push @alleles, $allele1;
			}
			my $allele_count = 0;
			my $skip = 0;
			foreach my $allele2 (@alleles) { #preliminary loop to check for site to skip completely
				if ( ($allele_count > 0) and (length $allele2 == $ref_length) ) { #complex substitution? Ignore
					$skip = 1;
				}
				$allele_count++;
			}
			$allele_count = 0;
			foreach my $allele2 (@alleles) {
				if ($allele_count > 0) {
					if (length $allele2 == 1) { #alt is single base
						if ( ($allele2 ne ".") and ($allele2 eq '*') ) { #base itself is deleted
							my $start = $cols[1] - $flank - 1;
							my $end = $cols[1] + (length $cols[3]) + $flank - 1;
							my $col_count = 0;
							foreach my $col (@cols) {
								if ($col_count > 8) {
									my @calls = split(/:/, $col);
									if ( ($calls[0] ne ".") and ($calls[0] == $allele_count) ) { #indel is present
										my $ID = $isolate_index{$col_count};
										my $mask = "$cols[0]\t$start\t$end";
										push(@{$indel_index{$ID}}, $mask);
									}
								}
								$col_count++;
							}
						}
						elsif ( ($allele2 ne ".") and ($allele2 ne '*') ) { #base is not deleted
							my $start = $cols[1] - $flank;
							my $end = $cols[1] + (length $cols[3]) + $flank - 1;
							my $col_count = 0;
							foreach my $col (@cols) {
								if ($col_count > 8) {
									my @calls = split(/:/, $col);
									if ( ($calls[0] ne ".") and ($calls[0] == $allele_count) ) { #indel is present
										my $ID = $isolate_index{$col_count};
										my $mask = "$cols[0]\t$start\t$end";
										push(@{$indel_index{$ID}}, $mask);
									}
								}
								$col_count++;
							}
						}
						if ( ($mode eq "rescue") and ($print_count == 0) and ($skip == 0) ) { #restore as invariant, can't be a SNP 
							my $base = substr($cols[3], 0, 1);
							$cols[3] = $base;
							$cols[4] = ".";
							print OUT join("\t", @cols), "\n";	
							$print_count++;
						}
					}
					else { #both ref and alt are multiple bp, complex variation, try to mask around location of indel 
						my $alt_length = length $allele2;
						if ($alt_length > $ref_length) { #insertion, could be at start or end of both so mask from each to be conservative
							my $start = $cols[1] - $flank - 1;
							my $end = $cols[1] + $ref_length + $flank - 1;
							my $col_count = 0;
							foreach my $col (@cols) {
								if ($col_count > 8) {
									my @calls = split(/:/, $col);
									if ( ($calls[0] ne ".") and ($calls[0] == $allele_count) ) { #indel is present
										my $ID = $isolate_index{$col_count};
										my $mask = "$cols[0]\t$start\t$end";
										push(@{$indel_index{$ID}}, $mask);
									}
								}
								$col_count++;
							}
							if ( ($mode eq "rescue") and ($print_count == 0) and ($skip == 0) ) { #restore as invariant, can't be a SNP 
								my $base = substr($cols[3], 0, 1);
								$cols[3] = $base;
								$cols[4] = ".";
								print OUT join("\t", @cols), "\n";
								$print_count++;	
							}
						}
						elsif ($alt_length < $ref_length) { #deletion, could be at start or end of both so mask from each to be conservative
							my $start = $cols[1] - $flank - 1;
							my $end = $cols[1] + $ref_length + $flank - 1;
							my $col_count = 0;
							foreach my $col (@cols) {
								if ($col_count > 8) {
									my @calls = split(/:/, $col);
									if ( ($calls[0] ne ".") and ($calls[0] == $allele_count) ) { #indel is present
										my $ID = $isolate_index{$col_count};
										my $mask = "$cols[0]\t$start\t$end";
										push(@{$indel_index{$ID}}, $mask);
									}
								}
								$col_count++;
							}	
							if ( ($mode eq "rescue") and ($print_count == 0) and ($skip == 0) ) { #restore as invariant, can't be a SNP 
								my $base = substr($cols[3], 0, 1);
								$cols[3] = $base;
								$cols[4] = ".";
								print OUT join("\t", @cols), "\n";
								$print_count++;	
							}						
						}
					}
				}
				$allele_count++;
			}	
		}
		elsif ( (length $cols[3] == 1) and (length $cols[4] == 1) ) { #normal site, print
			print OUT "$line\n";
		}	
		else {
			die "what's this? $line\n";
		}
	}
	else {
		print OUT "$line\n"; #print out comment lines
	}
}

close IN;
close OUT;

foreach my $isolate (keys %indel_index) {
	open (BED, ">$isolate.$suffix") or die;
    foreach my $mask_range (@{$indel_index{$isolate}}) {
        print BED "$mask_range\n";
    }
    close BED;
}

exit;

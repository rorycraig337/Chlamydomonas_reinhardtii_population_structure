#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qq(GetOptions);
use Data::Dumper;

#takes text file of consensus markers for two clades, a BED file of windows and a VCF as input
#calculates proportion of other clade marker SNPs in each window for each isolate
#usage: perl tally_marker_snps.pl --cladeA cladeA.txt --cladeB cladeB.txt --bed windows.bed --vcf snps.vcf --in consensus.txt --out out.tsv

my $cladeA;
my $cladeB;
my $in;
my $bed;
my $vcf;
my $out;

GetOptions(
	'cladeA=s' => \$cladeA,
	'cladeB=s' => \$cladeB,
	'in=s' => \$in,
	'bed=s' => \$bed,
	'vcf=s' => \$vcf,
	'out=s' => \$out,
) or die;

#get clade information

my %pop;

my $clade1_count = 0;
my $clade2_count = 0;

open (IN1, "$cladeA") or die;
open (IN2, "$cladeB") or die;

while (my $isolate1 = <IN1>) {
	chomp $isolate1;
	$pop{A}{$isolate1} = 1;
	$clade1_count++;
}

while (my $isolate2 = <IN2>) {
	chomp $isolate2;
	$pop{B}{$isolate2} = 1;
	$clade2_count++;
}

my $total_count = $clade1_count + $clade2_count;

close IN1;
close IN2;

#get marker positions

my %markers;

open (IN3, "$in") or die;

while (my $con = <IN3>) {
	chomp $con;
	my @con_cols = split(/\t/, $con);
	my $marker_site = "$con_cols[0]\t$con_cols[1]";
	my $marker_bases = "$con_cols[2]\t$con_cols[3]";
	$markers{$marker_site} = $marker_bases;
}

close IN3;

#get windows

my %windows;

open (IN4, "$bed") or die;

while (my $win = <IN4>) {
	chomp $win;
	my @win_cols = split(/\t/, $win);
	my $chr = $win_cols[0];
	my $win_bases = "$win_cols[1]\t$win_cols[2]";
	push @{ $windows{$chr} }, $win_bases;
}

close IN4;

#loop through VCF, foreach isolate matching bases to consensus for appropriate clade

my %pop_index;
my %isolate_index;

open (IN5, "$vcf") or die;

my $ind_index_check = 0;

while (my $vcf_line = <IN5>) {
	chomp $vcf_line;
	if ($vcf_line =~ /^#CHROM/) {
		my @head_cols = split(/\t/, $vcf_line);
		my $head_col_count = 0;
		#store column information foreach isolate
		foreach my $head_col (@head_cols) {
			if ($head_col_count > 8) {
				if (exists $pop{A}{$head_col}) {
					$pop_index{A}{$head_col_count} = $head_col;
				}
				elsif (exists $pop{B}{$head_col}) {
					$pop_index{B}{$head_col_count} = $head_col;
				}
			}
			$head_col_count++;
		}
		$ind_index_check = 1;
	}
	unless ($vcf_line =~ /^#/) {
		if ($ind_index_check == 0) {
			die "VCF file does not contain header with valid information\n";
		}
		my @vcf_cols = split(/\t/, $vcf_line);
		my $vcf_chr = $vcf_cols[0];
		my $vcf_base = $vcf_cols[1];
		my $vcf_site = "$vcf_cols[0]\t$vcf_cols[1]";
		my $ref = $vcf_cols[3];
		my $alt = $vcf_cols[4];
		#VCF site is a marker
		if (exists $markers{$vcf_site}) {
			my $con_site = $markers{$vcf_site};
			my @con_bases = split(/\t/, $con_site);
			my $A_con = $con_bases[0];
			my $B_con = $con_bases[1];
			my @window_match;
			#find appropriate window for site
			foreach my $chr_key (keys %windows) {
				if ($chr_key eq "$vcf_chr") {
					foreach my $win_array (@{$windows{$chr_key} } ) {
						my @win_cor = split(/\t/, $win_array);
						if ( ($vcf_base > $win_cor[0]) and ($vcf_base <= $win_cor[1]) ) {
							push @window_match, $win_array;
						}
				    }
				}
			}
			my $vcf_count = 0;
			foreach my $vcf_col (@vcf_cols) {
				#foreach isolate compare base to clade consensus
				if ( ($vcf_count > 8) and (exists $pop_index{A}{$vcf_count}) ) {
					my $isolate = $pop_index{A}{$vcf_count};
					my $isolate_base;
					unless ($vcf_col eq "null") {
						my @calls = split(/:/, $vcf_col);
						if ($calls[0] eq "0") {
							$isolate_base = $ref;
							if ($ref eq "$A_con") {
								foreach my $target_window (@window_match) {
									push @{ $isolate_index{$vcf_chr}{$target_window}{$isolate} }, "0"; 
								}
							}
							else {
								foreach my $target_window (@window_match) {
									push @{ $isolate_index{$vcf_chr}{$target_window}{$isolate} }, "1"; 
								}
							}
						}
						elsif ($calls[0] eq "1") {
							$isolate_base = $alt;
							if ($alt eq "$A_con") {
								foreach my $target_window (@window_match) {
									push @{ $isolate_index{$vcf_chr}{$target_window}{$isolate} }, "0"; 
								}
							}
							else {
								foreach my $target_window (@window_match) {
									push @{ $isolate_index{$vcf_chr}{$target_window}{$isolate} }, "1"; 
								}
							}
						}
					}
					else {
						foreach my $target_window (@window_match) {
							push @{ $isolate_index{$vcf_chr}{$target_window}{$isolate} }, "N";
						}
					}
				}	
				if ( ($vcf_count > 8) and (exists $pop_index{B}{$vcf_count}) ) {
					my $isolate = $pop_index{B}{$vcf_count};
					my $isolate_base;
					unless ($vcf_col eq "null") {
						my @calls = split(/:/, $vcf_col);
						if ($calls[0] eq "0") {
							$isolate_base = $ref;
							if ($ref eq "$B_con") {
								foreach my $target_window (@window_match) {
									push @{ $isolate_index{$vcf_chr}{$target_window}{$isolate} }, "0"; 
								}
							}
							else {
								foreach my $target_window (@window_match) {
									push @{ $isolate_index{$vcf_chr}{$target_window}{$isolate} }, "1"; 
								}
							}
						}
						elsif ($calls[0] eq "1") {
							$isolate_base = $alt;
							if ($alt eq "$B_con") {
								foreach my $target_window (@window_match) {
									push @{ $isolate_index{$vcf_chr}{$target_window}{$isolate} }, "0"; 
								}
							}
							else {
								foreach my $target_window (@window_match) {
									push @{ $isolate_index{$vcf_chr}{$target_window}{$isolate} }, "1"; 
								}
							}
						}
					}
 					else {
 						foreach my $target_window (@window_match) {
							push @{ $isolate_index{$vcf_chr}{$target_window}{$isolate} }, "N";
						}
					}
				}
				$vcf_count++;
			}
		}
	}
}

close IN5;

#print to file percentage of opposite clade SNPs foreach window

open (OUT, ">$out") or die;

foreach my $key (sort keys %isolate_index) {
	foreach my $sub_key (sort keys %{ $isolate_index{$key} } ) {
		print OUT "$key\t$sub_key";
		foreach my $sub2_key (sort keys %{ $isolate_index{$key}{$sub_key} } ) {
			my $length = scalar @{$isolate_index{$key}{$sub_key}{$sub2_key}};
			my $array_string = join("", @{$isolate_index{$key}{$sub_key}{$sub2_key}});
			my $same = ($array_string =~ tr/0//);
			my $other = ($array_string =~ tr/1//);
			my $snps = $same + $other;
			unless ($snps == 0) {
				my $percentage = ($other / $snps ) * 100;
				print OUT "\t$percentage";
			}
			else {
				print OUT "\tNA";
			}
		}
		print OUT "\n";
	}
}

close OUT;

exit;




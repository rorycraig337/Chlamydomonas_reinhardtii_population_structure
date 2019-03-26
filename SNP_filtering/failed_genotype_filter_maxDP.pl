#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qq(GetOptions);

#script removes sites where any individuals fail GQ and DP filters
#GQ is only filtered for SNPs, invariant sites are only filtered on depth
#DP_max is a two column tsv with isolate in column one, and mean converage in column two
#DP_max_int is used to change the max depth, which will be set as mean coverage + (DP_max_int * sqrt(mean coverage))
#usage: perl failed_genotype_masker_maxDP.pl --vcf in.vcf --GQ 20 --DP_min 3 --DP_max cov.tsv --DP_max_int 4 --out out.vcf

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
my $strain_number = 0;

while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /^#CHROM/) {
        my @head_cols = split(/\t/, $line);
        my $head_col_count = 0;
        my $total_col = scalar @head_cols;
        $strain_number = $total_col - 9;
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
        if ($cols[4] ne ".") { #site is a SNP
            $SNP_flag = 1;
        }
        my $alt_count = 0;
        my $pass_count = 0;
        my $depth_count = 0;
        my $ref_count = 0;
        my $allele_count = 0;
        foreach my $col (@cols) {
            no warnings 'numeric';
            if ($count > 8) {
                my @calls = split(/:/, $col);
                my $isolate_DP_max = $isolates{$count}; #for each isolate retrieve the associated max DP    
                if ($calls[0] eq "1") {
                    $alt_count++;
                    if ($SNP_flag == 1) {
                        if ($calls[3] >= $GQ) {
                            $pass_count++;
                        }
                    }
                    else {
                        print "WARNING: alt allele found but not a SNP? : $line\n";
                    }
                    if ( ($calls[2] >= $DP_min) and ($calls[2] <= $isolate_DP_max) ) {
                        $depth_count++;
                    }
                }
                elsif ($calls[0] eq "0") {
                    $ref_count++;
                    if ($SNP_flag == 1) { #is a SNP, check GQ
                        if ($calls[3] >= $GQ) {
                            $pass_count++;
                        }
                    }
                    else { #invariant
                        $pass_count++;
                    }
                    if ( ($calls[2] >= $DP_min) and ($calls[2] <= $isolate_DP_max) ) {
                        $depth_count++;
                    }
                }
                unless ($calls[0] eq "null") { #now check AD 
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
                	unless ( ($ad_count > 2) or ($ad_snp == 0) ) { #ignore previously multi-allelic, or sites with no informative reads
                    		my $ad_percent = ($ad_snp / $ad_total) * 100;
				if ($ad_percent >= 90) {
                        		$allele_count++;
                    		}
                	}
		}
            }
            $count++;
        }
        my $strain_count = $ref_count + $alt_count;
        if ( ($pass_count == $strain_number) and ($depth_count == $strain_number) and ($allele_count == $strain_number) and ($strain_count == $strain_number) ) {
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


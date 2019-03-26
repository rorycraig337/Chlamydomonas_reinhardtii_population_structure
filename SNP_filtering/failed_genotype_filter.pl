#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qq(GetOptions);

#script removes sites where any individuals fail GQ and DP filters
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

my $strain_number = 0;

while (my $line = <IN>) {
        chomp $line;
        if ($line =~ /^#CHROM/) {
                my @split = split(/\t/, $line);
                my $col_count = scalar @split;
                $strain_number = $col_count - 9;
        }
        unless ($line =~ /^#/) {
                my @cols = split(/\t/, $line);
                my $count = 0;
                my $alt_count = 0;
                my $pass_count = 0;
                my $depth_count = 0;
                my $ref_count = 0;
                foreach my $col (@cols) {
		    no warnings 'numeric';
                    if ($count > 8) {
                        my @calls = split(/:/, $col);
                        if ($calls[0] eq "1") {
                            $alt_count++;
                            if ($calls[3] >= $GQ) {
                                $pass_count++;
                            }
                            if ( ($calls[2] >= $DP_min) and ($calls[2] <= $DP_max) ) {
                                $depth_count++;
                            }
                        }
                        elsif ($calls[0] eq "0") {
                            $ref_count++;
                            if ($calls[3] >= $GQ) {
                                $pass_count++;
                            }
                            if ( ($calls[2] >= $DP_min) and ($calls[2] <= $DP_max) ) {
                                $depth_count++;
                            }
                        }
                    }
                    $count++;
                }
                my $strain_count = $ref_count + $alt_count;
                if ( ($pass_count == $strain_number) and ($depth_count == $strain_number) and ($strain_count == $strain_number) ) {
                    print OUT1 "$line\n";
                }
        }
    else {
                print OUT1 "$line\n";
        }
}

close IN;
close OUT1;

exit;



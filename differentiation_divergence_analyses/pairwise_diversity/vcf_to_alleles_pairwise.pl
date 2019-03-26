#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#script to convert VCF format to Rob's big table (.tsv) format
#optionally takes subset of samples, and outputs allele frequencies only for these samples
#usage: perl vcf_to_alleles.pl --vcf in.vcf --mode all/subset [--samples sample1,sample2] --prefix prefix

my $vcf;
my $mode;
my @samples = "null";
my $prefix;

GetOptions(
        'vcf=s' => \$vcf,
        'mode=s' => \$mode,
        'samples=s' => \@samples,
        'prefix=s' => \$prefix,
) or die "missing input\n";

@samples = split(/,/,join(',',@samples));

open (IN, "$vcf") or die;
open (OUT, ">$prefix.tsv") or die;	

my @strain_index;
my $sample_count;

while (my $line = <IN>) {
	chomp $line;
	if ($line =~ /^#CHROM/) {
		my @headers = split(/\t/, $line);
		my $header_count = 0;
		foreach my $header (@headers) {
			if ($header_count > 8) {
				if ($mode eq "all") {
					push @strain_index, $header_count;
				}
				elsif ($mode eq "subset") {
					if (grep /^$header$/, @samples) {
						push @strain_index, $header_count;
					}
				}
			}
			$header_count++;
		}
		$sample_count = scalar @strain_index;
	}
	unless ($line =~ /^#/) {
		my @cols = split(/\t/, $line);
		my $count = 0;
		my $A_count = 0;
		my $C_count = 0;
		my $G_count = 0;
		my $T_count = 0;
		my $ref = $cols[3];
		my $alt = $cols[4];
		foreach my $col (@cols) {
			no warnings 'numeric';
			if ( ($count > 8) and (grep /^$count$/, @strain_index) ) {
				my @calls = split(/:/, $col);
				if ($calls[0] eq "0") {
					if ($ref eq "A") {
						$A_count++;
					}
					elsif ($ref eq "C") {
						$C_count++;
					}
					elsif ($ref eq "G") {
						$G_count++;
					}
					elsif ($ref eq "T") {
						$T_count++;
					}
					else {
						die "ref base $ref at $cols[0]\t$cols[1]\n";
					}
				}
				elsif ($calls[0] eq "1") {
					if ($alt eq "A") {
						$A_count++;
					}
					elsif ($alt eq "C") {
						$C_count++;
					}
					elsif ($alt eq "G") {
						$G_count++;
					}
					elsif ($alt eq "T") {
						$T_count++;
					}
					else {
						die "alt base $alt at $cols[0]\t$cols[1]\n";
					}
				}
			}
			$count++;
		}
		my $total = $A_count + $C_count + $G_count + $T_count;
		unless ($total != $sample_count) {
			print OUT "$cols[0]\t$cols[1]\t$A_count:$C_count:$G_count:$T_count\n"
		}
	}
}


close IN;
close OUT;

exit;


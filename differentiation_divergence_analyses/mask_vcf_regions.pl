#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);
use Data::Dumper;

#script to mask regions from a VCF for a set of strains
#takes list of samples, and for each sample a bed file containing regions to mask
#will replace VCF record for that strain with the word "null"
#usage: perl mask_vcf_regions.pl --vcf in.vcf --samples sample1,sample2,sample3 --suffix .introgressed.bed --out masked.vcf

my $vcf;
my @samples;
my $suffix;
my $out;

GetOptions(
        'vcf=s' => \$vcf,
        'samples=s' => \@samples,
        'suffix=s' => \$suffix,
        'out=s' => \$out,
) or die "missing input\n";

@samples = split(/,/,join(',',@samples));

#make hash of hashes, sample->site

my %introgressed;

foreach my $sample (@samples) {
	my $file = "$sample$suffix";
	open (IN1, "$file") or next;
	while (my $line1 = <IN1>) {
		chomp $line1;
		my @cols1 = split(/\t/, $line1);
		my $chromosome = $cols1[0];
		my $start = $cols1[1] + 1; #convert 0-based to 1-based coordinate
		my $end = $cols1[2];
		my $count = $start;
		my $stop = $end + 1;
		until ($count == $stop) {
			my $mask = "$chromosome\t$count";
			$introgressed{$sample}{$mask} = 1; 
			$count++;
		}
	}
	print "read in $sample\n";
	close IN1;
}

open (IN2, "$vcf") or die;
open (OUT, ">$out") or die;

my %sample_index;
my $sample_count;

while (my $vcf_line = <IN2>) {
	chomp $vcf_line;
	if ($vcf_line =~ /^#CHROM/) {
		my @headers = split(/\t/, $vcf_line);
		my $header_count = 0;
		foreach my $header (@headers) {
			if ($header_count > 8) {
				foreach my $sample (@samples) {
					if ($header eq "$sample") {
						$sample_index{$header_count} = "$sample";
					}
				}
			}
			$header_count++;
		}
		$sample_count = scalar keys %sample_index;
		my $compare = scalar @samples;
		if ($sample_count != $compare) {
			print "WARNING: found fewer/more samples than expected in VCF\n";
		}
	}
	if ($vcf_line !~ /^#/) {
		my @vcf_cols = split(/\t/, $vcf_line);
		my $count = 0;
		my $site = "$vcf_cols[0]\t$vcf_cols[1]";
		foreach my $col (@vcf_cols) {
			no warnings 'numeric';
			if ( ($count > 8) and (exists $sample_index{$count}) ) {
				my $sample_check = $sample_index{$count};
				if (exists $introgressed{$sample_check} && exists $introgressed{$sample_check}{$site}) {
					$vcf_cols[$count] = "null";
				}
			}
			$count++;
		}
		print OUT join("\t", @vcf_cols),"\n";
	}
	else {
		print OUT "$vcf_line\n";
	}
}

close IN2;
close OUT;

exit;







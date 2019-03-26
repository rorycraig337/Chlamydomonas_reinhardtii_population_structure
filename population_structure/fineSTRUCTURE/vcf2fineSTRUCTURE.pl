#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#script converts VCF with no missing data to input format used by fineSTRUCTURE
#usage: perl vcf2fineSTRUCTURE.pl --vcf in.vcf --chr chromosomes.txt --suffix phase

my $vcf;
my $chr;
my $suffix;

GetOptions(
        'vcf=s' => \$vcf,
        'chr=s' => \$chr,
        'suffix=s' => \$suffix,
) or die "missing input\n";   

open (IN, "$vcf") or die;

my %snp_index;
my %haplo_index;
my %iso_index;
my $snp_count = 0;
my $total;
my $curr_chr = 0;

while (my $line = <IN>) {
	chomp $line;
	if ($line =~ /^#CHROM/) {
		my @head_cols = split(/\t/, $line);
		my $head_col_count = 0;
		foreach my $head_col (@head_cols) {
			if ($head_col_count > 8) {
				$iso_index{$head_col_count} = $head_col;
			}
			$head_col_count++;
		}
		$total = $head_col_count - 9;
	}
	elsif ($line !~ /#/) {
		my @cols = split(/\t/, $line);
		if ($cols[0] ne "$curr_chr") {
			$snp_count = 0;
		}
		$curr_chr = $cols[0];
		$snp_count++;
		$snp_index{$cols[0]}{$snp_count} = $cols[1];
		my $col_count = 0;
		foreach my $col (@cols) {
			if ($col_count > 8) {
				my @calls = split(/:/, $col);
				my $iso = $iso_index{$col_count};
				if ($calls[0] eq "0") {
					push(@{$haplo_index{$cols[0]}{$iso}}, 0);
				}
				elsif ($calls[0] eq "1") {
					push(@{$haplo_index{$cols[0]}{$iso}}, 1);
				}
				else {
					die;
				}
			}
			$col_count++;
		}
	}
}

close IN;

open (IN2, "$chr") or die;

while (my $chromosome = <IN2>) {
	chomp $chromosome;
	open (OUT, ">$chromosome.$suffix") or die;
	print OUT "$total\n";
	my @numbers = sort { $snp_index{$chromosome}{$a} <=> $snp_index{$chromosome}{$b} } keys %{$snp_index{$chromosome}};
	my $total_snps = $numbers[-1];
	print OUT "$total_snps\n";
	print OUT "P";
	foreach my $snp_pos (sort {$a <=> $b} keys %{$snp_index{$chromosome}}) {
		print OUT " $snp_index{$chromosome}{$snp_pos}";
	}
	print OUT "\n";
	foreach my $iso_haplo (sort keys %{$haplo_index{$chromosome}}) {
	    	foreach my $snp_call (@{$haplo_index{$chromosome}{$iso_haplo}}){
    			print OUT "$snp_call";
    		}
    		print OUT "\n";
	}
	close OUT;
}


close IN2;
close OUT;

exit;

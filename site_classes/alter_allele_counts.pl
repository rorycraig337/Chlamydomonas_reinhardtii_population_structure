#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#script to take Rob's annotation table, and replace alleles with an arbitrary SNP if one exists in VCF, and an invariant site if one doesn't
#does not actually add allele information, only a dummy SNP (0:10:0:10)
#VCF should contain only SNPs
#usage: perl alter_allele_counts.pl --vcf in.vcf --table in.tsv --out out.tsv

my $vcf;
my $table;
my $out;

GetOptions(
        'vcf=s' => \$vcf,
        'table=s' => \$table,
        'out=s' => \$out,
) or die "missing input\n";

open (IN1, "$vcf") or die;

my %snps;

while (my $vcf_line = <IN1>) {
	chomp $vcf_line;
	unless ($vcf_line =~ /^#/) {
		my @vcf_cols = split(/\t/, $vcf_line);
		my $match = "$vcf_cols[0]\t$vcf_cols[1]";
		$snps{$match} = 1;
	}
}

close IN1;

open (IN2, "$table") or die;
open (OUT, ">$out") or die;

while (my $line = <IN2>) {
	chomp $line;
	unless ($line =~ /^#/) {
		my @cols = split(/\t/, $line);
		my $remove = pop @cols;
		my $hit = "$cols[0]\t$cols[1]";
		if (exists $snps{$hit}) {
			push @cols, "0:10:0:10";
		}
		else {
			push @cols, "0:20:0:0";
		}
		print OUT join("\t", @cols),"\n";
	}	
	else {
		print OUT "$line\n";
	}
}

close IN2;
close OUT;

exit;

#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qq(GetOptions);
use Bio::AlignIO;
use Bio::SimpleAlign;

#takes a fasta alignment of SNPs for isolates and calculates consensus sequences for given clades
#percentage is the threshold required for a consensus base to be called
#assumes SNPs have been filtered and there are appropriate numbers of isolates at each site, as percentage is used regardless of count
#usage: perl fasta2consensus.pl --clade1 cladeA.txt --clade2 cladeB.txt --fasta snps.fa --percentage 60 --out consensus.txt 

my $clade1;
my $clade2;
my $fasta;
my $percentage;
my $out;

GetOptions(
	'clade1=s' => \$clade1,
	'clade2=s' => \$clade2,
	'percentage=i' => \$percentage,
	'fasta=s' => \$fasta,
	'out=s' => \$out,
) or die;

my @isolates1;
my @isolates2;

my $clade1_count = 0;
my $clade2_count = 0;

open (IN1, "$clade1") or die;
open (IN2, "$clade2") or die;

while (my $isolate1 = <IN1>) {
	chomp $isolate1;
	push @isolates1, $isolate1;
	$clade1_count++;
}

while (my $isolate2 = <IN2>) {
	chomp $isolate2;
	push @isolates2, $isolate2;
	$clade2_count++;
}

my $total_count = $clade1_count + $clade2_count;

close IN1;
close IN2;

open (OUT, ">$out") or die;

my $str = Bio::AlignIO->new(-file => "$fasta",
	-format => "fasta");
my $aln = $str->next_aln();
my $len = $aln->length();
my $pos = 1;
my $end = $len + 1;
until ($pos == $end) { #loop through alignment
	my @clade1_bases;
	my @clade2_bases;
	foreach my $clade1_seq (@isolates1) {
		my $seq = $aln->get_seq_by_id($clade1_seq);
		my $base = $seq->subseq($pos, $pos);
		push @clade1_bases, $base;
	}
	foreach my $clade2_seq (@isolates2) {
		my $seq = $aln->get_seq_by_id($clade2_seq);
		my $base = $seq->subseq($pos, $pos);
		push @clade2_bases, $base;
	}
	my $A1 = 0;
	my $C1 = 0;
	my $G1 = 0;
	my $T1 = 0;
	my $A2 = 0;
	my $C2 = 0;
	my $G2 = 0;
	my $T2 = 0;
	foreach my $base1 (@clade1_bases) {
		if ($base1 eq "A") {
			$A1++;
		}
		elsif ($base1 eq "C") {
			$C1++;
		}
		elsif ($base1 eq "G") {
			$G1++;
		}
		elsif ($base1 eq "T") {
			$T1++;
		}
	}
	foreach my $base2 (@clade2_bases) {
		if ($base2 eq "A") {
			$A2++;
		}
		elsif ($base2 eq "C") {
			$C2++;
		}
		elsif ($base2 eq "G") {
			$G2++;
		}
		elsif ($base2 eq "T") {
			$T2++;
		}
	}
	my $clade1_site_count = $A1 + $C1 + $G1 + $T1;
	my $clade2_site_count = $A2 + $C2 + $G2 + $T2;
	my $total_site_count = $clade1_site_count + $clade2_site_count;
	if ( ( ($A1 / $clade1_site_count) * 100) > $percentage) {
		print OUT "$pos\tA\t";
	}
	elsif ( ( ($C1 / $clade1_site_count) * 100) > $percentage) {
		print OUT "$pos\tC\t";
	}
	elsif ( ( ($G1 / $clade1_site_count) * 100) > $percentage) {
		print OUT "$pos\tG\t";
	}
	elsif ( ( ($T1 / $clade1_site_count) * 100) > $percentage) {
		print OUT "$pos\tT\t";
	}
	else {
		print OUT "$pos\tN\t";
	}
	if ( ( ($A2 / $clade2_site_count) * 100) > $percentage) {
		print OUT "A\n";
	}
	elsif ( ( ($C2 / $clade2_site_count) * 100) > $percentage) {
		print OUT "C\n";
	}
	elsif ( ( ($G2 / $clade2_site_count) * 100) > $percentage) {
		print OUT "G\n";
	}
	elsif ( ( ($T2 / $clade2_site_count) * 100) > $percentage) {
		print OUT "T\n";
	}
	else {
		print OUT "N\n";
	}
	if ($pos % 10000 == 0) {
		print "$pos\n";
	}
	$pos++;
}

close OUT;

exit;

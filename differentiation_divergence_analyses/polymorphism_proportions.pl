#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qq(GetOptions);
use Bio::AlignIO;
use Bio::SimpleAlign;

#script to compute number of fixed/shared/private alleles from a multiple fasta alignment containing two clades/comparison groups
#assumptions only work for bi-allelic sites
#outputs simple proportions (.txt file), and additionally detailed table of each SNP and its category for each clade (.tsv file)
#outputs changed by mode, simple is proportions only, detailed is proportions and SNP table
#usage: perl poymorphism_proportions.pl --clade1 list1.txt --clade2 list2.txt --fasta in.fasta --out prefix --mode simple/detailed

my $clade1;
my $clade2;
my $minA;
my $minB;
my $fasta;
my $out;
my $mode = "default";

GetOptions(
	'clade1=s' => \$clade1,
	'clade2=s' => \$clade2,
	'minA=i' => \$minA,
	'minB=i' => \$minB,
	'fasta=s' => \$fasta,
	'out=s' => \$out,
	'mode=s' => \$mode,
) or die;

if ($mode eq "default") {
	$mode = "simple";
}
elsif ($mode eq "detailed") {
	open (OUT1, ">$out.clade1.tsv") or die;
	open (OUT2, ">$out.clade2.tsv") or die;
}

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

my $str = Bio::AlignIO->new(-file => "$fasta",
	-format => "fasta");
my $aln = $str->next_aln();
my $len = $aln->length();
my $pos = 1;
my $fixed = 0;
my $shared = 0;
my $private1 = 0;
my $private2 = 0;
my $invariant = 0;
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
		elsif ($base1 ne "-") {
			print "$pos\t";
			print join(" ", @clade1_bases);
			die;
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
		elsif ($base2 ne "-") {
			print "$pos\t";
            		print join(" ", @clade2_bases);
            		die;
		}
	}
	my $clade1_site_count = $A1 + $C1 + $G1 + $T1;
	my $clade2_site_count = $A2 + $C2 + $G2 + $T2;
	my $total_site_count = $clade1_site_count + $clade2_site_count;
	if ( ($clade1_site_count >= $minA) and ($clade2_site_count >= $minB) ) {
		if ( ($A1 == $clade1_site_count) or ($C1 == $clade1_site_count) or ($G1 == $clade1_site_count) or ($T1 == $clade1_site_count) ) { #clade 1 is invariant
			if ( ($A2 == $clade2_site_count) or ($C2 == $clade2_site_count) or ($G2 == $clade2_site_count) or ($T2 == $clade2_site_count) ) { #clade 2 is also invariant
				if ( (($A1 + $A2) == $total_site_count) or (($C1 + $C2) == $total_site_count) or (($G1 + $G2) == $total_site_count) or (($T1 + $T2) == $total_site_count) ) { #same base in each clade
					$invariant++;
				}
				else { #different base in each clade
					$fixed++;
					if ($mode eq "detailed") {
						if ($A1 > 0) {
							print OUT1 "$pos\tA\tfixed\n";
						}
						elsif ($C1 > 0) {
							print OUT1 "$pos\tC\tfixed\n";
						}
						elsif ($G1 > 0) {
							print OUT1 "$pos\tG\tfixed\n";
						}
						elsif ($T1 > 0) {
							print OUT1 "$pos\tT\tfixed\n";
						}
						else {
							die "$pos\n";
						}
						if ($A2 > 0) {
							print OUT2 "$pos\tA\tfixed\n";
						}
						elsif ($C2 > 0) {
							print OUT2 "$pos\tC\tfixed\n";
						}
						elsif ($G2 > 0) {
							print OUT2 "$pos\tG\tfixed\n";
						}
						elsif ($T2 > 0) {
							print OUT2 "$pos\tT\tfixed\n";
						}
						else {
							die "$pos\n";
						}
					}
				}
			}
			else { #clade 2 is variant but clade 1 is not
				$private2++;
				if ($mode eq "detailed") {
					if ($A1 > 0) {
						if ( ($A2 > 0) and ($C2 > 0) ) {
							print OUT2 "$pos\tC\tprivate\n";
						}
						elsif ( ($A2 > 0) and ($G2 > 0) ) {
							print OUT2 "$pos\tG\tprivate\n";
						}
						elsif ( ($A2 > 0) and ($T2 > 0) ) {
							print OUT2 "$pos\tT\tprivate\n";
						}
						else {
							die "$pos\n";
						}
					}
					elsif ($C1 > 0) {
						if ( ($C2 > 0) and ($A2 > 0) ) {
							print OUT2 "$pos\tA\tprivate\n";
						}
						elsif ( ($C2 > 0) and ($G2 > 0) ) {
							print OUT2 "$pos\tG\tprivate\n";
						}
						elsif ( ($C2 > 0) and ($T2 > 0) ) {
							print OUT2 "$pos\tT\tprivate\n";
						}
						else {
							die "$pos\n";
						}
					}
					elsif ($G1 > 0) {
						if ( ($G2 > 0) and ($A2 > 0) ) {
							print OUT2 "$pos\tA\tprivate\n";
						}
						elsif ( ($G2 > 0) and ($C2 > 0) ) {
							print OUT2 "$pos\tC\tprivate\n";
						}
						elsif ( ($G2 > 0) and ($T2 > 0) ) {
							print OUT2 "$pos\tT\tprivate\n";
						}
						else {
							die "$pos\n";
						}
					}
					elsif ($T1 > 0) {
						if ( ($T2 > 0) and ($A2 > 0) ) {
							print OUT2 "$pos\tA\tprivate\n";
						}
						elsif ( ($T2 > 0) and ($C2 > 0) ) {
							print OUT2 "$pos\tC\tprivate\n";
						}
						elsif ( ($T2 > 0) and ($G2 > 0) ) {
							print OUT2 "$pos\tG\tprivate\n";
						}
						else {
							die "$pos\n";
						}
					}
					else {
						die "$pos\n";
					}
				}
			}
		}
		elsif ( ($A2 == $clade2_site_count) or ($C2 == $clade2_site_count) or ($G2 == $clade2_site_count) or ($T2 == $clade2_site_count) ) { #clade 2 is invariant, and clade 1 is variant
			$private1++;
			if ($mode eq "detailed") {
				if ($A2 > 0) {
					if ( ($A1 > 0) and ($C1 > 0) ) {
						print OUT1 "$pos\tC\tprivate\n";
					}
					elsif ( ($A1 > 0) and ($G1 > 0) ) {
						print OUT1 "$pos\tG\tprivate\n";
					}
					elsif ( ($A1 > 0) and ($T1 > 0) ) {
						print OUT1 "$pos\tT\tprivate\n";
					}
					else {
						die "$pos\n";
					}
				}
				if ($C2 > 0) {
					if ( ($C1 > 0) and ($A1 > 0) ) {
						print OUT1 "$pos\tA\tprivate\n";
					}
					elsif ( ($C1 > 0) and ($G1 > 0) ) {
						print OUT1 "$pos\tG\tprivate\n";
					}
					elsif ( ($C1 > 0) and ($T1 > 0) ) {
						print OUT1 "$pos\tT\tprivate\n";
					}
					else {
						die "$pos\n";
					}
				}
				if ($G2 > 0) {
					if ( ($G1 > 0) and ($A1 > 0) ) {
						print OUT1 "$pos\tA\tprivate\n";
					}
					elsif ( ($G1 > 0) and ($C1 > 0) ) {
						print OUT1 "$pos\tC\tprivate\n";
					}
					elsif ( ($G1 > 0) and ($T1 > 0) ) {
						print OUT1 "$pos\tT\tprivate\n";
					}
					else {
						die "$pos\n";
					}
				}
				if ($T2 > 0) {
					if ( ($T1 > 0) and ($A1 > 0) ) {
						print OUT1 "$pos\tA\tprivate\n";
					}
					elsif ( ($T1 > 0) and ($C1 > 0) ) {
						print OUT1 "$pos\tC\tprivate\n";
					}
					elsif ( ($T1 > 0) and ($G1 > 0) ) {
						print OUT1 "$pos\tG\tprivate\n";
					}
					else {
						die "$pos\n";
					}
				}
			}
		}
		elsif ( ($A1 != $clade1_site_count) and ($C1 != $clade1_site_count) and ($G1 != $clade1_site_count) and ($T1 != $clade1_site_count) ) { #clade 1 is variant, clade 2 is variant
			$shared++;
			if ($mode eq "detailed") {
				if ( ($A1 > 0) and ($C1 > 0) ) {
					print OUT1 "$pos\tA C\tshared\n";
					print OUT2 "$pos\tA C\tshared\n";
				}
				elsif ( ($A1 > 0) and ($G1 > 0) ) {
					print OUT1 "$pos\tA G\tshared\n";
					print OUT2 "$pos\tA G\tshared\n";
				}
				elsif ( ($A1 > 0) and ($T1 > 0) ) {
					print OUT1 "$pos\tA T\tshared\n";
					print OUT2 "$pos\tA T\tshared\n";
				}
				elsif ( ($C1 > 0) and ($G1 > 0) ) {
					print OUT1 "$pos\tC G\tshared\n";
					print OUT2 "$pos\tC G\tshared\n";
				}
				elsif ( ($C1 > 0) and ($T1 > 0) ) {
					print OUT1 "$pos\tC T\tshared\n";
					print OUT2 "$pos\tC T\tshared\n";
				}
				elsif ( ($G1 > 0) and ($T1 > 0) ) {
					print OUT1 "$pos\tG T\tshared\n";
					print OUT2 "$pos\tG T\tshared\n";
				}
				else {
					die "$pos\n";
				}
			}
		}
	}
	undef @clade1_bases;
	undef @clade2_bases;
	if ($pos % 10000 == 0) {
		print "$pos\n";
	}
	$pos++;
}


close OUT1;
close OUT2;

open (OUT3, ">$out.txt") or die;

print OUT3 "total\t$len\ninvariant\t$invariant\nfixed\t$fixed\nshared\t$shared\nprivate1\t$private1\nprivate2\t$private2\n";

close OUT3;

exit;




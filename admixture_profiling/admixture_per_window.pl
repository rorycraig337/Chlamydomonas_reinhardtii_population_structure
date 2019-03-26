#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#script calculates proportion of admixture per window/chromosome
#accounts for IBD by subtracting IBD tracts from cumulative length of windows
#admixture tracts (that have had IBD tracts subtracted) can then be matched to window for each isolate, and sum of admixture compared to cumulative window
#list contains isolates of interest
#usage: perl admixture_per_window.pl --windows windows.bed --IBD_suffix IBD.bed --admixture_suffix introgressed.IBD_sub.bed --list list.txt --out admixture_windows.tsv

my $windows;
my $IBD_suffix;
my $admixture_suffix;
my $list;
my $out;

GetOptions(
	'windows=s' => \$windows,
	'IBD_suffix=s' => \$IBD_suffix,
	'admixture_suffix=s' => \$admixture_suffix,
	'list=s' => \$list,
	'out=s' => \$out,
) or die "missing input\n";

my %IBD;
my %admixture;

open (IN1, "$list") or die;

my $iso_count = 0;

while (my $iso = <IN1>) {
	chomp $iso;
	open (IN2, "$iso.$IBD_suffix") or die;
	while (my $IBD_line = <IN2>) {
		chomp $IBD_line;
		push(@{$IBD{$iso}}, $IBD_line); #create hash with array of all IBD tracts
	}
	close IN2;
	open (IN3, "$iso.$admixture_suffix") or die;
	while (my $admix_line = <IN3>) {
		chomp $admix_line;
		push(@{$admixture{$iso}}, $admix_line); #create hash with array of all admixture tracts
	}
	close IN3;
	$iso_count++;
}

close IN1;

open (IN4, "$windows") or die;
open (OUT, ">$out") or die;

while (my $line = <IN4>) {
	chomp $line;
	my @cols = split(/\t/, $line);
	my $length = $cols[2] - $cols[1];
	my $cumulative_length = $length * $iso_count;
	my $cumulative_IBD = 0;
	my $cumulative_admixture = 0;
	foreach my $isolate (keys %IBD) {
		foreach my $IBD_tract (@{$IBD{$isolate}}) {
			my @tract_cols = split(/\t/, $IBD_tract);
			if ($cols[0] eq "$tract_cols[0]") { #same chromsome
				unless ( ($tract_cols[1] >= $cols[2]) or ($tract_cols[2] <= $cols[1]) ) { #tract is not out of range
					if ( ($tract_cols[1] < $cols[1]) and ($tract_cols[2] < $cols[2]) and ($tract_cols[2] > $cols[1]) ) { #tract right flank is within window
						$cumulative_IBD += $tract_cols[2] - $cols[1];
					}
					elsif ( ($tract_cols[1] > $cols[1]) and ($tract_cols[1] < $cols[2]) and ($tract_cols[2] > $cols[2]) ) { #tract left flank is within window
						$cumulative_IBD += $cols[2] - $tract_cols[1];
					}
					elsif ( ($tract_cols[1] >= $cols[1]) and ($tract_cols[1] < $cols[2]) and ($tract_cols[2] <= $cols[2]) and ($tract_cols[2] > $cols[1]) ) { #tract is entirely within window
						$cumulative_IBD += $tract_cols[2] - $tract_cols[1];
					}
					elsif ( ($tract_cols[1] <= $cols[1]) and ($tract_cols[2] >= $cols[2]) ) { #tract completely spans window
						$cumulative_IBD += $cols[2] - $cols[1];
					}
					else {
						die "what else is there!: $line\n$IBD_tract\n";
					}
				}
			}
		}
	}
	foreach my $isolate (keys %admixture) {
		foreach my $admix_tract (@{$admixture{$isolate}}) {
			my @tract_cols = split(/\t/, $admix_tract);
			if ($cols[0] eq "$tract_cols[0]") { #same chromsome
				unless ( ($tract_cols[1] >= $cols[2]) or ($tract_cols[2] <= $cols[1]) ) { #tract is not out of range
					if ( ($tract_cols[1] < $cols[1]) and ($tract_cols[2] < $cols[2]) and ($tract_cols[2] > $cols[1]) ) { #tract right flank is within window
						$cumulative_admixture += $tract_cols[2] - $cols[1];
					}
					elsif ( ($tract_cols[1] > $cols[1]) and ($tract_cols[1] < $cols[2]) and ($tract_cols[2] > $cols[2]) ) { #tract left flank is within window
						$cumulative_admixture += $cols[2] - $tract_cols[1];
					}
					elsif ( ($tract_cols[1] >= $cols[1]) and ($tract_cols[1] < $cols[2]) and ($tract_cols[2] <= $cols[2]) and ($tract_cols[2] > $cols[1]) ) { #tract is entirely within window
						$cumulative_admixture += $tract_cols[2] - $tract_cols[1];
					}
					elsif ( ($tract_cols[1] <= $cols[1]) and ($tract_cols[2] >= $cols[2]) ) { #tract completely spans window
						$cumulative_admixture += $cols[2] - $cols[1];
					}
					else {
						die "what else is there!: $line\n$admix_tract\n";
					}
				}
			}
		}
	}
	print "$line\t$cumulative_length\t$cumulative_IBD\t$cumulative_admixture\n";
	my $unique_length = $cumulative_length - $cumulative_IBD;
	my $admix_percent = ($cumulative_admixture / $unique_length) * 100;
	print OUT "$line\t$admix_percent\n";
}

close IN4;
close OUT;

exit;

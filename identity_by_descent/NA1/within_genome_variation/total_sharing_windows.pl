#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);
use Data::Dumper;

#script to calculate per window total sharing (i.e. mean IBD over all pairs)
#usage: perl total_sharing_windows.pl --pairs 325 --IBD ibd.txt --windows windows.txt --out total_sharing.txt

my $pairs;
my $IBD;
my $windows;
my $out;

GetOptions(
	'pairs=i' => \$pairs,
	'IBD=s' => \$IBD,
	'windows=s' => \$windows,
	'out=s' => \$out,
) or die "missing input\n"; 

my %index;

open (IN1, "$IBD") or die;

while (my $line = <IN1>) {
	chomp $line;
	my @cols = split(/\t/, $line);
	my $chr = "chromosome_$cols[2]";
	push @{$index{$chr}}, "$cols[2]\t$cols[3]\t$cols[4]";
}

close IN1;

print Dumper(\%index);

open (IN2, "$windows") or die;
open (OUT, ">$out") or die;

while (my $window = <IN2>) {
	chomp $window;
	my @wcols = split(/\t/, $window);
	my $length = $wcols[2] - $wcols[1];
	my $cumulative_length = $length * $pairs;
	my $cumulative_IBD = 0;
	foreach my $tract (@{$index{$wcols[0]}}){
		my @tcols = split (/\t/, $tract);
		unless ( ($tcols[1] >= $wcols[2]) or ($tcols[2] <= $wcols[1]) ) { #tract is not out of range
			if ( ($tcols[1] < $wcols[1]) and ($tcols[2] < $wcols[2]) and ($tcols[2] > $wcols[1]) ) { #tract right flank is within window
				$cumulative_IBD += $tcols[2] - $wcols[1];
			}
			elsif ( ($tcols[1] > $wcols[1]) and ($tcols[1] < $wcols[2]) and ($tcols[2] > $wcols[2]) ) { #tract left flank is within window
				$cumulative_IBD += $wcols[2] - $tcols[1];
			}
			elsif ( ($tcols[1] >= $wcols[1]) and ($tcols[1] < $wcols[2]) and ($tcols[2] <= $wcols[2]) and ($tcols[2] > $wcols[1]) ) { #tract is entirely within window
				$cumulative_IBD += $tcols[2] - $tcols[1];
			}
			elsif ( ($tcols[1] <= $wcols[1]) and ($tcols[2] >= $wcols[2]) ) { #tract completely spans window
				$cumulative_IBD += $wcols[2] - $wcols[1];
			}
			else {
				die "what else is there!: $window\n$tract\n";
			}
		}
	}
	my $IBD_percent = ($cumulative_IBD / $cumulative_length) * 100;
	print OUT "$window\t$IBD_percent\n";
}

close IN2;
close OUT;

exit;


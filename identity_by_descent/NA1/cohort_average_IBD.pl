#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);
use Data::Dumper;

#calculate cohort average IBD sharing for each isolate
#usage: perl cohort_average_IBD.pl --lineage NA1.txt --IBD NA1.IBD_tracts.100kb.txt --genome 106481163 --out NA1_cohort_average.txt

my $lineage;
my $IBD;
my $genome;
my $out;

GetOptions(
        'lineage=s' => \$lineage,
        'IBD=s' => \$IBD,
        'genome=i' => \$genome,
        'out=s' => \$out,
) or die "missing input\n";

open (IN1, "$lineage") or die;

my %index;

while (my $iso = <IN1>) {
	chomp $iso;
	open (IN2, "$IBD") or die;
	while (my $line = <IN2>) {
		chomp $line;
		my @cols = split(" ", $line);
		if ( ($cols[0] eq "$iso") or ($cols[1] eq "$iso") ) {
			my $length = $cols[4] - $cols[3] - 1;
			if ($cols[0] eq "$iso") {
				$index{$cols[0]}{$cols[1]} += $length;
			}
			elsif ($cols[1] eq "$iso") {
				$index{$cols[1]}{$cols[0]} += $length;
			}
		}
	}
	close IN2;
}

close IN1;

open (OUT, ">$out") or die;

foreach my $iso_key (sort keys %index) {
	my @totals;
	foreach my $match (sort keys %{$index{$iso_key}}) {
	    	my $per = ($index{$iso_key}{$match} / $genome) * 100;
		print "$iso_key\t$match\t$per\n";
    		push @totals, $per;
	}
	my $length = scalar @totals;
	print "$length\n";
	my $sum;
	foreach my $total (@totals) {
		$sum += $total;
	}
	my $final = $sum / $length;
	print OUT "$iso_key\t$final\n";
}

close OUT;

exit;




#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qq(GetOptions);
use Data::Dumper;

#script takes output of GATK DepthOfCoverage and caluclates mean coverage for each isolate
#only sites with a minimum of $min_DP reads will be included in calculation
#usage: perl mean_coverage.pl --cov in.cov --DP_min 3 --out means.txt

my $cov;
my $DP_min;
my $out;

GetOptions(
    'cov=s' => \$cov,
    'DP_min=i' => \$DP_min,
    'out=s' => \$out,
) or die "missing input\n";

open (IN, "$cov") or die;

my %isolate_index;
my %cov_tally;

my $header = <IN>;
chomp $header;
my @head_cols = split(/\t/, $header);
my $head_count = 0;
foreach my $head_col (@head_cols) {
	if ($head_count > 2) {
		my $isolate = substr($head_col, 10);
		$isolate_index{$head_count} = $isolate; #make index hash of isolate ID and associated column number
		$cov_tally{$isolate} = "0\t0"; #initialise hash containing two values, for total coverage and total sites
	}
	$head_count++;
}


while (my $line = <IN>) {
	chomp $line;
	my @cols = split(/\t/, $line);
	my $col_count = 0;
	foreach my $col (@cols) {
		if ($col_count > 2) {
			my $ID = $isolate_index{$col_count}; #get isolate ID from index hash
			if ($col >= $DP_min) {
				my @cov_info = split(/\t/, $cov_tally{$ID});
				$cov_info[0] += $col; #tally total coverage
				$cov_info[1] += 1; #tally number of sites
				$cov_tally{$ID} = "$cov_info[0]\t$cov_info[1]"; #store new values in hash
			}
		}
		$col_count++;
	}
}

close IN;

print Dumper(\%cov_tally);

open (OUT, ">$out") or die;

foreach my $iso (sort keys %cov_tally) {
	my @cov_results = split(/\t/, $cov_tally{$iso}); #access final tally of coverage and sites
	my $mean = $cov_results[0] / $cov_results[1];
	print OUT "$iso\t$mean\n";
}

close OUT;

exit;

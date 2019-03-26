#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#filters elements in bed file my size 'min_size'
#usage: perl filter_bed_by_size.pl --bed in.bed --min_size 50000 --out out.bed

my $bed;
my $min_size;
my $out;

GetOptions(
	'bed=s' => \$bed,
	'min_size=i' => \$min_size,
	'out=s' => \$out,
) or die "missing input\n";

open (IN, "$bed") or die;
open (OUT, ">$out") or die;

while (my $line = <IN>) {
	chomp $line;
	my @cols = split(/\t/, $line);
	my $length = $cols[2] - $cols[1];
	if ($length >= $min_size) {
		print OUT "$line\n";
	}
}

close IN;
close OUT;

exit;

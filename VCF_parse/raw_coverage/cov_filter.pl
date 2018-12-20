#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#filters cov file based on a bed file
#usage: perl cov_filter.pl --cov cov.txt --bed Ns.txt --out filtered_cov.txt

my $cov;
my $bed;
my $out;

GetOptions(
        'cov=s' => \$cov,
        'bed=s' => \$bed,
        'out=s' => \$out,
) or die "missing input\n";

open (IN1, "$bed") or die;

my %N_index;

while (my $pos = <IN1>) {
	chomp $pos;
	my @info = split(/\t/, $pos);
	for (my $i=$info[1]+1; $i <= $info[2]; $i++) {
		my $match = "$info[0]:$i";
		$N_index{$match} = 1;
	}
}

close IN1;

open (IN2, "$cov") or die;
open (OUT, ">$out") or die;

my $header = <IN2>;
print OUT "$header\n";

while (my $line = <IN2>) {
	chomp $line;
	my @cols = split(" ", $line);
	unless (exists $N_index{$cols[0]}) {
		print OUT "$line\n";
	}
}

close IN2;
close OUT;

exit;

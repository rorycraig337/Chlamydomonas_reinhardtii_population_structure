#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#validates output of vcf2hmmIBD.pl, removing non-SNP lines
#usage: perl validate_hmmIBD.pl --in in.txt --out out.txt

my $in;
my $out;

GetOptions(
        'in=s' => \$in,
        'out=s' => \$out,
) or die "missing input\n";

open (IN, "$in") or die;
open (OUT, ">$out") or die;

my $header = <IN>;
print OUT "$header";

while (my $line = <IN>) {
	chomp $line;
	my @cols = split(/\t/, $line);
	my $total = 0;
	my $ref = 0;
	my $alt = 0;
	my $col_count = 0;
	foreach my $col (@cols) {
		if ($col_count > 1) {
			if ( ($col == 1) or ($col == 0) ) {
				$total++;
				if ($col == 1) {
					$alt++;
				}
				elsif ($col == 0) {
					$ref++;
				}
			}
		}
		$col_count++;
	}
	if ( ($total > 0) and ($alt != $total) and ($ref != $total) ) {
		print OUT "$line\n";
	}
}

close IN;
close OUT;

exit;

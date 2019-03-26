#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#given a fasta file computes GC content
#usage: perl get_GC.pl --fasta in.fasta

my $fasta;

GetOptions(
        'fasta=s' => \$fasta,
) or die "missing input\n";

my $length = 0;
my $gc_count = 0;

open (IN, "$fasta") or die;

while (my $line = <IN>) {
	chomp $line;
	if ($line !~ /^>/) {
		$length += length $line;
		$gc_count += $line =~ tr/GCgc//;
	}
}

my $percent = ($gc_count / $length) * 100;

print "GC content of $fasta is $percent %\n";

close IN;

exit;

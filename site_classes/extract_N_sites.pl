#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#uses the big chlamy annotation file to produce a bed file of all Ns (should be merged afterwards)
#usage: perl extract_N_sites.pl --annotation annotation_table.txt --Nbed Ns.bed --Nsite Ns.txt

my $annotation;
my $Nbed;
my $Nsite;

GetOptions(
        'annotation=s' => \$annotation,
        'Nbed=s' => \$Nbed,
	'Nsite=s' => \$Nsite,
) or die "missing input\n";

open (IN, "$annotation") or die;
open (OUT1, ">$Nbed") or die;
open (OUT2, ">$Nsite") or die;

while (my $line = <IN>) {
	unless ($line =~ /^#/) {
		chomp $line;
		my @cols = split(/\t/, $line);
		if ( ($cols[2] eq "N") or ($cols[2] eq "n") ) {
			my $start = $cols[1] - 1;
			print OUT1 "$cols[0]\t$start\t$cols[1]\n";
			print OUT2 "$cols[0]\t$cols[1]\n";
		}
	}
}

close IN;
close OUT1;
close OUT2;

exit;

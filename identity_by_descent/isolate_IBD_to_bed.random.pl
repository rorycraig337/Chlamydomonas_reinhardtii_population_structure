#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#script takes IBD tracts file and file of isolates, producing per isolate BED file of IBD tracts for masking (one pair is chosen randomly)
#usage: perl isolate_IBD_to_bed.random.pl --IBD tracts.txt --temp temp.txt --isolates list.txt --suffix IBD.bed

my $IBD;
my $temp;
my $isolates;
my $suffix;

GetOptions(
        'IBD=s' => \$IBD,
 	'temp=s' => \$temp,
        'isolates=s' => \$isolates,
	'suffix=s' => \$suffix,
) or die "missing input\n";

open (IBD, "$IBD") or die;
open (TEMP, ">$temp") or die;

while (my $ibd = <IBD>) {
	chomp $ibd;
	my @ibd_cols = split(/\t/, $ibd);
	my $int = int(rand(2)); #randomly selects 0 or 1
	print TEMP "$ibd_cols[$int]\t$ibd_cols[2]\t$ibd_cols[3]\t$ibd_cols[4]\t$ibd_cols[5]\t$ibd_cols[6]\n";
}

close IBD;
close TEMP;

open (IN1, "$isolates") or die;

while (my $isolate = <IN1>) {
	chomp $isolate;
	open (IN2, "$temp") or die;
	open (OUT, ">$isolate.$suffix") or die;
	while (my $line = <IN2>) {
		chomp $line;
		my @cols = split(/\t/, $line);
		if ($cols[0] eq "$isolate") {
			my $start = $cols[2] - 1;
			print OUT "chromosome_$cols[1]\t$start\t$cols[3]\n";
		}
	} 
	close IN2;
	close OUT;
}

close IN1;

system("rm $temp");

exit;



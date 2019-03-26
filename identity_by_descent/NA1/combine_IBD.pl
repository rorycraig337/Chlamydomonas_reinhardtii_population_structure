#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);
use Data::Dumper;

#combine hmmIBD estimate, and 100 kb and 500 kb estimates of proportion of genome IBD
#usage: perl combine_IBD.pl --genome 106481163 --pairs NA1_pairs.txt --fract NA1_fract.txt --one NA1_cumulative.100kb.txt --five NA1_cumulative.500kb.txt --out NA1_IBD_summary.txt

my $genome;
my $pairs;
my $fract;
my $one;
my $five;
my $out;

GetOptions(
	'genome=i' => \$genome,
        'pairs=s' => \$pairs,
        'fract=s' => \$fract,
        'one=s' => \$one,
        'five=s' => \$five,
        'out=s' => \$out,
) or die "missing input\n";

open (IN1, "$fract") or die;

my %f_index;

while (my $line1 = <IN1>) {
	chomp $line1;
	my @cols1 = split(" ", $line1);
	unless ($line1 =~ /^sample/) {
		my $pair = "$cols1[0]\t$cols1[1]";
		my $f = $cols1[9] * 100; 
		$f_index{$pair} = $f; 
	}
}

close IN1;

open (IN2, "$one") or die;

my %one_index;

while (my $line2 = <IN2>) {
	chomp $line2;
	my @cols2 = split(" ", $line2);
	my $pair = "$cols2[0]\t$cols2[1]";
	my $f = $cols2[2] / $genome * 100;
	$one_index{$pair} = $f;
}

close IN2;

open (IN3, "$five") or die;

my %five_index;

while (my $line3 = <IN3>) {
	chomp $line3;
	my @cols3 = split(" ", $line3);
	my $pair = "$cols3[0]\t$cols3[1]";
	my $f =	$cols3[2] / $genome * 100;
        $five_index{$pair} = $f;
}

close IN3;

open (IN, "$pairs") or die;
open (OUT, ">$out") or die;


while (my $pair = <IN>) {
	chomp $pair;
	my @cols = split(" ", $pair);
	my $i1 = $f_index{"$cols[0]\t$cols[1]"};
	my $i2 = $one_index{"$cols[0]\t$cols[1]"};
	my $i3 = $five_index{"$cols[0]\t$cols[1]"};
	print OUT "$cols[0]\t$cols[1]\t$i1\t$i2\t$i3\n";
}

close IN;
close OUT;

exit;





exit;

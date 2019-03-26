#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);
use Data::Dumper;

#takes a matrix output by MEGA, and a list of isolates, and converts to csv file for use with R
#list of isolates needs clade ID in second column and sampling site ID in third column (tab-separated)
#distances is tab separated file containing distance between sampling sites
#perl matrix_converter.pl --matrix in.matrix --isolates all.txt --out out.csv

my $matrix;
my $isolates;
my $distances;
my $out;

GetOptions(
        'matrix=s' => \$matrix,
        'isolates=s' => \$isolates,
        'distances=s' => \$distances,
        'out=s' => \$out,
) or die "missing input\n";

my %iso_index;

open (IN1, "$isolates") or die;

my $iso_count = 0;

while (my $iso = <IN1>) {
	chomp $iso;
	$iso_count++;
	$iso_index{$iso_count} = $iso;
}

close IN1;

my %distance_index;

open (IN2, "$distances") or die;

while (my $pair = <IN2>) {
	chomp $pair;
	my @pair_cols = split(" ", $pair);
	my $site = "$pair_cols[0]\t$pair_cols[1]";
	$distance_index{$site} = $pair_cols[2];
}

close IN2;

open (IN3, "$matrix") or die;
open (OUT, ">$out") or die;

my $header = <IN3>;
my $count = 0;

while (my $line = <IN3>) {
	chomp $line;
	$count++; #isolate 1 in matrix
	my @cols = split(" ", $line);
	if ( ($count < 10) and ($count != 1) ) { #these cols have the isolate number as [ *], so remove first part element of @cols
		shift @cols;
		my $id = substr $cols[0], 0, -1; #remove other "]"
		my @isolate1_info = split(/\t/, $iso_index{$count}); #get isolate and site in array
		shift @cols; #now remove id col
		my $col_number = scalar @cols;
		my $col_count = 0;
		foreach my $gen_distance (@cols) {
			$col_count++; #col count is now isolate 2 in matrix
			my @isolate2_info = split(/\t/, $iso_index{$col_count});
			my $iso1 = $isolate1_info[0];
			my $iso2 = $isolate2_info[0];
			my $site1 = $isolate1_info[2];
			my $site2 = $isolate2_info[2];
			my $clade1 = $isolate1_info[1];
			my $clade2 = $isolate2_info[1];
			my $geo_distance = "null";
			foreach my $distance_pair (keys %distance_index) {
				my @distance_pairs = split(/\t/, $distance_pair);
				if ( ( ($site1 eq "$distance_pairs[0]") and ($site2 eq "$distance_pairs[1]") ) or ( ($site1 eq "$distance_pairs[1]") and ($site2 eq "$distance_pairs[0]") ) ) {
					$geo_distance = $distance_index{$distance_pair};
				}
			}
			if ($geo_distance eq "null") {
				die "error: couldn't find distance\n";
			}
			my $type;
			if ( ($clade1 eq "A") and ($clade2 eq "A") ) {
				$type = "0,0";
			}
			elsif ( ( ($clade1 eq "A") and ($clade2 eq "B") ) or ($clade2 eq "A") and ($clade1 eq "B") ) {
				$type = "0,1";
			}
			elsif ( ($clade1 eq "B") and ($clade2 eq "B") ) {
				$type = "1,1";
			}
			else {
				die "error: couldn't find clade\n";
			}
			print OUT "$iso1,$iso2,$gen_distance,$geo_distance,\"$type\"\n";
		}
	}
	elsif ($count != 1) {
		my $temp = substr $cols[0], 1; #remove first [
		my $id = substr $temp, 0, -1; #remove last ]
		my @isolate1_info = split(/\t/, $iso_index{$count}); #get isolate and site in array
		shift @cols; #now remove id col
		my $col_number = scalar @cols;
		my $col_count = 0;
		foreach my $gen_distance (@cols) {
			$col_count++; #col count is now isolate 2 in matrix
			my @isolate2_info = split(/\t/, $iso_index{$col_count});
			my $iso1 = $isolate1_info[0];
			my $iso2 = $isolate2_info[0];
			my $site1 = $isolate1_info[2];
			my $site2 = $isolate2_info[2];
			my $clade1 = $isolate1_info[1];
			my $clade2 = $isolate2_info[1];
			my $geo_distance = "null";
			foreach my $distance_pair (keys %distance_index) {
				my @distance_pairs = split(/\t/, $distance_pair);
				if ( ( ($site1 eq "$distance_pairs[0]") and ($site2 eq "$distance_pairs[1]") ) or ( ($site1 eq "$distance_pairs[1]") and ($site2 eq "$distance_pairs[0]") ) ) {
					$geo_distance = $distance_index{$distance_pair};
				}
			}
			if ($geo_distance eq "null") {
				die "error: couldn't find distance\n";
			}
			my $type;
			if ( ($clade1 eq "A") and ($clade2 eq "A") ) {
				$type = "0,0";
			}
			elsif ( ( ($clade1 eq "A") and ($clade2 eq "B") ) or ($clade2 eq "A") and ($clade1 eq "B") ) {
				$type = "0,1";
			}
			elsif ( ($clade1 eq "B") and ($clade2 eq "B") ) {
				$type = "1,1";
			}
			else {
				die "error: couldn't find clade\n";
			}
			print OUT "$iso1,$iso2,$gen_distance,$geo_distance,\"$type\"\n";
		}
	}
}

close IN3;
close OUT;

exit;

#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qq(GetOptions);
use Data::Dumper;

#formats admixture data for R heatmap
#isolates is list of isolates in original order
#isoaltes_new is list of isolates as it appears in phylogeny 
#outputs per chromosome file
#usage: perl format_admixture_for_R_heat.pl --isolates all.txt --isolates_new all_phylo.txt --in 20kb.R.tsv --suffix 20kb.R.csv

my $isolates;
my $isolates_new;
my $in;
my $suffix;

GetOptions(
	'isolates=s' => \$isolates,
	'isolates_new=s' => \$isolates_new,
	'in=s' => \$in,
	'suffix=s' => \$suffix,
) or die;

my %isolate_index;
my $counter1 = 0;

open (IN1, "$isolates") or die;

while (my $isolate = <IN1>) {
	chomp $isolate;
	$counter1++;
	$isolate_index{old}{$counter1} = "$isolate";
}

close IN1;

my $counter2 = 0;

open (IN2, "$isolates_new") or die;

while (my $isolate = <IN2>) {
	chomp $isolate;
	$counter2++;
	$isolate_index{new}{$counter2} = "$isolate";
}

close IN2;

my %index;
my %chr_index;
my $line_count = 0;

open (IN3, "$in") or die;

while (my $line = <IN3>) {
	chomp $line;
	my @cols = split(/\t/, $line);
	my $col_count = 0;
	if ($line_count == 0) {
		foreach my $col (@cols) {
			$col_count++;
			my @win_info = split(" ", $col);
			my $chr = $win_info[0];
			push @{ $chr_index{$chr} }, $col_count;
			$index{windows}{$col_count} = $col;
		}
	}
	else {
		$col_count = 0;
		my $id = $isolate_index{old}{$line_count};
		foreach my $col (@cols) {
			$col_count++;
			$index{$id}{$col_count} = $col;
		}
	}
	$line_count++;
}

close IN3;

foreach my $chrom (sort keys %chr_index) {
	open (OUT, ">$chrom.$suffix") or die;
	print OUT "isolates";
	foreach my $window_num (@{$chr_index{$chrom}}){
		my $window = $index{windows}{$window_num};
		print OUT ",$window";
	}
	print OUT "\n";
	foreach my $isolate_num (sort { $a <=> $b } keys %{ $isolate_index{new} }) {
		my $isolate_id = $isolate_index{new}{$isolate_num};
		print OUT "$isolate_id";
		foreach my $col_num (@{$chr_index{$chrom}}){
			my $value = $index{$isolate_id}{$col_num};
			print OUT ",$value";
		}
		print OUT "\n";
	}
	close OUT;
}

exit;


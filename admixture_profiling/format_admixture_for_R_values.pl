#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qq(GetOptions);
use Data::Dumper;

#formats admixture data for R, can be used as input for format_admixture_for_R_heat.pl
#isolates is list of isolates as they appear in admixutre output in col1, and clade in cols2 (tab separated)
#usage: perl format_admixture_for_R_values.pl --isolates all_info.txt --in f_20kb.introgression.tsv --out 20kb.R.tsv

my $isolates;
my $in;
my $out;

GetOptions(
	'isolates=s' => \$isolates,
	'in=s' => \$in,
	'out=s' => \$out,
) or die;

my %isolate_index;

open (IN1, "$isolates") or die;

my $counter = 0;

while (my $isolate = <IN1>) {
	chomp $isolate;
	my @i_info = split(/\t/, $isolate);
	$counter++;
	if ($i_info[1] eq "A") {
		$isolate_index{$counter} = "$isolate";
	}
	elsif ($i_info[1] eq "B") {
		$isolate_index{$counter} = "$isolate";
	}
	else {
		die "error in isolates file\n";
	}
}

close IN1;

my %index;
my $line_count = 0;

open (IN2, "$in") or die;

while (my $line = <IN2>) {
	chomp $line;
	$line_count++;
	my @cols = split(/\t/, $line);
	my $window = "$cols[0] $cols[1] $cols[2]";
	$index{windows}{$line_count} = $window;
	my $col_count = 0;
	foreach my $col (@cols) {
		$col_count++;
		if ($col_count > 3) {
			my $isolate_number = $col_count - 3;
			my $id = $isolate_index{$isolate_number};
			my @id_info = split(/\t/, $id);
			my $name = $id_info[0];
			my $clade = $id_info[1];
			if ( ($col ne "NA") ) {
				if ($clade eq "A") {
					my $A_value = $col / 100;
					$index{$name}{$line_count} = $A_value;
				}
				elsif ($clade eq "B") {
					my $B_value = 1 - ($col / 100);
					$index{$name}{$line_count} = $B_value;
				}
			}
			else {
				if ($clade eq "A") {
					$index{$name}{$line_count} = "NA";
				}
				elsif ($clade eq "B") {
					$index{$name}{$line_count} = "NA";
				}
			}
		}
	}
}

close IN2;

open (OUT, ">$out") or die;

foreach my $win_count (sort { $a <=> $b } keys %{ $index{windows} } ) {
		unless ($win_count == 1) {
			print OUT "\t";
		}
		my $win = $index{windows}{$win_count};
		print OUT "$win";
}

print OUT "\n";

foreach my $key (sort keys %index) {
	unless ($key eq "windows") {
		foreach my $key_count (sort { $a <=> $b } keys %{ $index{$key} } ) {
			unless ($key_count == 1) {
				print OUT "\t";
			}
			my $call = $index{$key}{$key_count};
			print OUT "$call";		
		}
	}
	print OUT "\n";
}


close OUT;

exit;



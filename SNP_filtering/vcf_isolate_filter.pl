#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#script to remove certain isolates from a VCF (a list given as --exclude)
#usage: perl vcf_isolate_filter.pl --vcf in.vcf --exclude exclude.txt --out out.vcf

my $in;
my $exclude;
my $out;

GetOptions(
        'in=s' => \$in,
        'exclude=s' => \$exclude,
        'out=s' => \$out,
) or die "missing input\n";

my %exc1;

open (LIST, "$exclude") or die;

while (my $iso = <LIST>) {
	chomp $iso;
	$exc1{$iso} = 1; #store isolate IDs to be excluded
}

close LIST;

my %exc2;

open (IN, "$in") or die;
open (OUT, ">$out") or die;

while (my $line = <IN>) {
	chomp $line;
	if ($line =~ /^#CHROM/) {
		my @head_cols = split(/\t/, $line);
		my $head_col_count = 0;
		foreach my $head_col (@head_cols) {
			if ($head_col_count > 8) {
				if (exists $exc1{$head_col}) {
					$exc2{$head_col_count} = 1; #store column numbers to be excluded
				}
				else {
					print OUT "\t$head_col"; #print header of isolate IDs to be retained
				}	
			}
                        elsif ($head_col_count == 0) {
                               	print OUT "$head_col";
                        }
                        else {
                                print OUT "\t$head_col";
                        }
			$head_col_count++;
		}
		print OUT "\n";
	}
	elsif ($line =~ /^#/) {
		print OUT "$line\n";
	}
	unless ($line =~ /^#/) {
		my @cols = split(/\t/, $line);
		my $col_count = 0;
		foreach my $col (@cols) {
			if ($col_count > 8) {
				unless (exists $exc2{$col_count}) {
					print OUT "\t$col";
				}
			}
			elsif ($col_count == 0) {
				print OUT "$col";
			}
			else {
				print OUT "\t$col";
			}
			$col_count++;
		}
		print OUT "\n";	
	}
}

close IN;
close OUT;

exit;

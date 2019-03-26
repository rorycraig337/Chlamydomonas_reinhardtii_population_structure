#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#takes a CDS table as output by extract_site_classes.pl, and removes sense-antisense or out of frame CDS overlaps
#usage: perl CDS_overlap_filter.pl --table CDS_table.txt --out CDS_no_overlaps_table.txt

my $table;
my $out;

GetOptions(
        'table=s' => \$table,
        'out=s' => \$out,
) or die "missing input\n";

open (IN, "$table") or die;
open (OUT, ">$out") or die;

while (my $line = <IN>) {
        chomp $line;
        unless ($line =~ /^#/) {
                my $pass = 0;
                my @cols = split(/\t/, $line);
                my @strands = split(/:/, $cols[21]); #strand (+,-) present at site
                my $strand_count = scalar @strands;
                my @frames = split(/:/, $cols[22]); #reading frames (0,1,2) present at site
                if ($strand_count > 1) { #possible sense-antisense overlap
                        my $plus = 0;
                        my $minus = 0;
                        foreach my $strand_elem (@strands) {
                                if ($strand_elem eq "+") {
                                        $plus = 1;
                                }
                                if ($strand_elem eq "-") {
                                        $minus = 1;
                                }
                        }
                        if ( ($plus == 1) and ($minus == 1) ) {
                                $pass = 1; #sense-antisense overlap
                        	print "sense-antisense: $line\n";
			}
                        else { #same strand overlap
                                my $initial = $frames[0];
                                my $in_frame = 0;
                                foreach my $frame_elem (@frames) {
                                        if ($frame_elem != "$initial") { #break in frame overlap
                                                $in_frame = 1;
					}
                                }
                                if ($in_frame == 1) {
                                        $pass = 1; #overlap same strand, but not same frame
                                	print "frame-break: $line\n";
				}
                        }
                }
                if ($pass == 0) {
                        print OUT "$line\n";
                }
        }
}

close IN;
close OUT;

exit;

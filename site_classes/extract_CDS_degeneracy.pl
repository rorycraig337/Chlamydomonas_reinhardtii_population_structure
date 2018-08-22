#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#script to take CDS table with no overlaps or multi-allelic codons, and extract 0D, 2D and 4D sites to bed files
#optional unique flag will output only sites that are uniquely CDS (i.e. with no overlap with UTR/introns etc.)
#usage: perl extract_CDS_degeneracy.pl --table CDS_table_final.txt --unique

my $table;
my $unique = 0;

GetOptions(
        'table=s' => \$table,
        'unique' => \$unique,
) or die "missing input\n";

open (IN, "$table") or die;
open (OUT1, ">4D.bed") or die;
open (OUT2, ">2D.bed") or die;
open (OUT3, ">0D.bed") or die;


while (my $line = <IN>) {
        chomp $line;
        my @cols = split(/\t/, $line);
        if ($unique == 1) {
                if ( ($cols[10] == 1) and ($cols[9] == 0) and ($cols[11] == 0) and ($cols[12] == 0) and ($cols[3] == 1) and ($cols[4] == 1) and ($cols[13] == 1) and ($cols[14] == 1) and ($cols[5] == 0) and ($cols[6] == 0) and ($cols[7] == 0) and ($cols[8] == 0) and ($cols[15] == 0) and ($cols[16] == 0) ) {
                        my $start = $cols[1] - 1;
                        print OUT1 "$cols[0]\t$start\t$cols[1]\t$cols[19]\t.\t$cols[21]\n";
                }
                if ( ($cols[11] == 1) and ($cols[9] == 0) and ($cols[10] == 0) and ($cols[12] == 0) and ($cols[3] == 1) and ($cols[4] == 1) and ($cols[13] == 1) and ($cols[14] == 1) and ($cols[5] == 0) and ($cols[6] == 0) and ($cols[7] == 0) and ($cols[8] == 0) and ($cols[15] == 0) and ($cols[16] == 0) ) {
                        my $start = $cols[1] -1;
                        print OUT2 "$cols[0]\t$start\t$cols[1]\t$cols[19]\t.\t$cols[21]\n";
                }
                if ( ($cols[9] == 1) and ($cols[10] == 0) and ($cols[11] == 0) and ($cols[12] == 0) and ($cols[3] == 1) and ($cols[4] == 1) and ($cols[13] == 1) and ($cols[14] == 1) and ($cols[5] == 0) and ($cols[6] == 0) and ($cols[7] == 0) and ($cols[8] == 0) and ($cols[15] == 0) and ($cols[16] == 0) ) {
                        my $start = $cols[1] -1;
                        print OUT3 "$cols[0]\t$start\t$cols[1]\t$cols[19]\t.\t$cols[21]\n";
                }
        }
        elsif ($unique == 0) {
                if ( ($cols[10] == 1) and ($cols[9] == 0) and ($cols[11] == 0) and ($cols[12] == 0) and ($cols[3] == 1) and ($cols[4] == 1) and ($cols[13] == 1) and ($cols[14] == 1) ) {
                        my $start = $cols[1] - 1;
                        print OUT1 "$cols[0]\t$start\t$cols[1]\t$cols[19]\t.\t$cols[21]\n";
                }
                if ( ($cols[11] == 1) and ($cols[9] == 0) and ($cols[10] == 0) and ($cols[12] == 0) and ($cols[3] == 1) and ($cols[4] == 1) and ($cols[13] == 1) and ($cols[14] == 1) ) {
                        my $start = $cols[1] -1;
                        print OUT2 "$cols[0]\t$start\t$cols[1]\t$cols[19]\t.\t$cols[21]\n";
                }
                if ( ($cols[9] == 1) and ($cols[10] == 0) and ($cols[11] == 0) and ($cols[12] == 0) and ($cols[3] == 1) and ($cols[4] == 1) and ($cols[13] == 1) and ($cols[14] == 1) ) {
                        my $start = $cols[1] -1;
                        print OUT3 "$cols[0]\t$start\t$cols[1]\t$cols[19]\t.\t$cols[21]\n";
                }
        }
}

close IN;
close OUT1;
close OUT2;
close OUT3;

exit;

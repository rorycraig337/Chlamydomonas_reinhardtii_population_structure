#!/usr/bin/perl

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# perl chlamy_popgen.pl <input> <table> <prefix> [options]      # 
# ============================================================= #
# calculates population genetics statistics from Rob's          #
# annotation table format (i.e. allele frequencies A:C:G:T)     #
# ============================================================= #
# output is to <prefix>.<mode>.txt                              #
# ============================================================= #
# --input_format: tsv or vcf (default tsv)                      #
# tsv is standard annotation table with allele frequencies,     #
# alternatively provide a VCF file from which a temporary table #
# will be produced (temp.tsv). This must contain only invariant #
# sites and bi-allelic SNPs, with all INDELs etc. filtered      #
# ============================================================= #
# --bed_format : BED3 or BED6 (default BED3)                    #
# standard 0-based BED format, if BED6 option to report         #
# output per feature (e.g. gene ID) from column 4               #
# ============================================================= #
# --mode : feature and/or overall (comma separated)             #
# defines output, statistics per feature or for all sites	#
# ============================================================= #
# --output : pi, theta, TajD, folded-SFS, Fst and/or dxy        #
# defines statistics, folded-SFS is for max copy number only,   #
# Fst and dxy are for two populations only                      #
# ============================================================= #
# --max_copies: integer                                         #
# number of alleles if all individuals have good base calls,    #
# for 2 populations give 2 values separated by comma            #
# ============================================================= #
# --min_copies: integer (default is max_copies)                 #   
# for 2 populations give 2 values separated by comma            #
# will calculate statistics for range of copy numbers,          #
# pi (and Fst/dxy) can be calculated directly for all copy      # 
# numbers (default), or can be taken as a weighted average.     #
# Theta and Tajima's D are always taken as weighted averages    #
# in min_copies is less than max_copies                         #
# ============================================================= #
# --pi_weighted_average                                         #
# pi will be calculated as a weighted average for each copy     #
# number (see above), this flag affects Fst and dxy if they     #
# are being calculated                                          #
# ============================================================= #
# --features: by_line or by_annotation                          #
# will calculate statistics per line (e.g. a 5 kb window)       #
# or per annotation (BED6 only) in column 4 (e.g. gene ID)      #
# only compatible with --mode feature                           #
# ============================================================= #
# --populations: 1 or 2 (default is 1)                          #
# can provide table where last two columns contain allele       #
# frequencies for two populations, enabling calculation of      #
# Fst and dxy between populations                               #
# ============================================================= #
# --pop_names: name1,name2 (default is pop1,pop2)               #
# set population names if two population mode is used           #
# ============================================================= #
# --pop_list: file.txt                                          #
# for use with --input vcf and --populations 2, a two line text #
# file with a list of comma separated individuals, one line per #
# population                                                    #
# ============================================================= #
# --codon_usage                                                 #
# additionally outputs SFS for P->P, U->U and P->U/U->P         #
# segregating sites (based on G/C preferred codons)             #
# only compatible with --output folded-SFS                      #
# ============================================================= #
# example:                                                      #
# perl chlamy_popgen.pl 4D.bed alleles.tsv 4D --bed_format \    #
# BED3 --mode feature,overall --output pi,theta,folded-SFS  \   #
# --max_copies 18 --min_copies 12 --features by_line            #
# ============================================================= #
# Rory Craig, May 2018                                          #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#per feature with codon usage currently not implemented

use warnings;
use strict;
use Getopt::Long qw(GetOptions);
use Data::Dumper;

my $bed = $ARGV[0];
my $table = $ARGV[1];
my $prefix = $ARGV[2];

my $max_copies;
my $min_copies = "default";
my $input_format = "tsv";
my @modes;
my @outputs;
my $pi_weighted_average = 0;
my $format = "BED3";
my $feature_method = "none";
my $populations = 1;
my @pop_names = ("pop1", "pop2");
my $pop_list = "none";
my $codon_usage = 0;

GetOptions(
	'max_copies=s' => \$max_copies,
	'min_copies=s' => \$min_copies,
	'input_format=s' => \$input_format,
	'mode=s' => \@modes, 
	'output=s' => \@outputs,
	'pi_weighted_average' => \$pi_weighted_average,
	'bed_format=s' => \$format,
	'features=s' => \$feature_method,
	'populations=i' => \$populations,
	'pop_names=s' => \@pop_names,
	'pop_list=s' => \$pop_list,
	'codon_usage' => \$codon_usage,  
) or die "missing input\n";

my @min_copies;
my @max_copies = split(/,/,join(',',$max_copies));

if ($min_copies eq "default") {
	@min_copies = @max_copies;
}
else {
	@min_copies = split(/,/,join(',',$min_copies));
}

if ( ($populations != 1) and ($populations != 2) ) {
	print "WARNING: populations was set as $populations, only 1 or 2 is allowed, setting as 1...\n";
	$populations = 1;
}

if ($populations == 1) {
	if (scalar @min_copies > 1) {
		print "WARNING: populations set to 1, but 2 values given for min_copies. The first value will be used, is this correct?\n";
	}
	if (scalar @max_copies > 1) {
		print "WARNING: populations set to 1, but 2 values given for max_copies. The first value will be used, is this correct?\n";
	}
}

if (scalar @pop_names > 2) {
	shift @pop_names;
	shift @pop_names;
	@pop_names = split(/,/,join(',',@pop_names));
	if (scalar @pop_names != 2) {
		print "WARNING: pop_names not given two names, default names will be used\n";
		@pop_names = ("pop1", "pop2");
	}
}

################ section 1: parse VCF if used as input file ################

if ($input_format eq "vcf") {
	print "input format is vcf, parsing vcf for allele frequencies...\n";
	if ( ($populations == 2) and ($pop_list eq "none") ) {
		print "WARNING: vcf file provided and populations set as 2, but populations file is missing. Running as a single population\n";
		$populations = 1;
	}
	my %pop_hash;
	my %pop_hash_index;
	if ( ($populations == 2) and ($pop_list ne "none") ) {
		open (POP, "$pop_list") or die "pop list file name invalid\n";
		my $pop1_inds = <POP>;
		my $pop2_inds = <POP>;
		chomp $pop1_inds;
		chomp $pop2_inds;
		my @pop1_list = split(/,/, $pop1_inds);
		my @pop2_list = split(/,/, $pop2_inds);
		foreach my $ind1 (@pop1_list) {
			$pop_hash{pop1}{$ind1} = 1;
		}
		foreach my $ind2 (@pop2_list) {
			$pop_hash{pop2}{$ind2} = 1;
		}
		close POP;	
	}
	open (VCF, "$table") or die "VCF file name invalid\n";
	open (TEMP, ">temp.tsv") or die;
	my $ind_index_check = 0;
	while (my $vcf_line = <VCF>) {
		chomp $vcf_line;
		if ($vcf_line =~ /^#CHROM/) {
			my @head_cols = split(/\t/, $vcf_line);
			my $head_col_count = 0;
			foreach my $head_col (@head_cols) {
				if ($head_col_count > 8) {
					if (exists $pop_hash{pop1}{$head_col}) {
						$pop_hash_index{pop1}{$head_col_count} = 1;
					}
					elsif (exists $pop_hash{pop2}{$head_col}) {
						$pop_hash_index{pop2}{$head_col_count} = 1;
					}
				}
				$head_col_count++;
			}
			$ind_index_check = 1;
		}
		unless ($vcf_line =~ /^#/) {
			if ($ind_index_check == 0) {
				die "VCF file does not contain header with valid information\n";
			}
			my @vcf_cols = split(/\t/, $vcf_line);
			if ($populations == 1) {
				my $vcf_count = 0;
				my $A = 0;
				my $C = 0;
				my $G = 0;
				my $T = 0;
				my $ref = $vcf_cols[3];
				my $alt = $vcf_cols[4];
				foreach my $vcf_col (@vcf_cols) {
					if ($vcf_count > 8) {
						my @calls = split(/:/, $vcf_col);
						if ($calls[0] eq "0") {
							if ($ref eq "A") {
								$A++;
							}
							elsif ($ref eq "C") {
								$C++;
							}
							elsif ($ref eq "G") {
								$G++;
							}
							elsif ($ref eq "T") {
								$T++;
							}
							else {
								die "VCF reference base error, line: $vcf_line\n";
							}
						}
						elsif ($calls[0] eq "1") {
							if ($alt eq "A") {
								$A++;
							}
							elsif ($alt eq "C") {
								$C++;
							}
							elsif ($alt eq "G") {
								$G++;
							}
							elsif ($alt eq "T") {
								$T++;
							}
							else {
								die "VCF alternative base error, line: $vcf_line\n";
							}	
						}
					}
					$vcf_count++;
				}
				print TEMP "$vcf_cols[0]\t$vcf_cols[1]\t$A:$C:$G:$T\n";
			}
			if ($populations == 2) {
				my $vcf_count = 0;
				my $A1 = 0;
				my $A2 = 0;
				my $C1 = 0;
				my $C2 = 0;
				my $G1 = 0;
				my $G2 = 0;
				my $T1 = 0;
				my $T2 = 0;
				my $ref = $vcf_cols[3];
				my $alt = $vcf_cols[4];
				foreach my $vcf_col (@vcf_cols) {
					if ( ($vcf_count > 8) and (exists $pop_hash_index{pop1}{$vcf_count}) ) {
						my @calls = split(/:/, $vcf_col);
						if ($calls[0] eq "0") {
							if ($ref eq "A") {
								$A1++;
							}
							elsif ($ref eq "C") {
								$C1++;
							}
							elsif ($ref eq "G") {
								$G1++;
							}
							elsif ($ref eq "T") {
								$T1++;
							}
							else {
								die "VCF reference base error, line: $vcf_line\n";
							}
						}
						elsif ($calls[0] eq "1") {
							if ($alt eq "A") {
								$A1++;
							}
							elsif ($alt eq "C") {
								$C1++;
							}
							elsif ($alt eq "G") {
								$G1++;
							}
							elsif ($alt eq "T") {
								$T1++;
							}
							else {
								die "VCF alternative base error, line: $vcf_line\n";
							}	
						}	
					}
					elsif ( ($vcf_count > 8) and (exists $pop_hash_index{pop2}{$vcf_count}) ) {
						my @calls = split(/:/, $vcf_col);
						if ($calls[0] eq "0") {
							if ($ref eq "A") {
								$A2++;
							}
							elsif ($ref eq "C") {
								$C2++;
							}
							elsif ($ref eq "G") {
								$G2++;
							}
							elsif ($ref eq "T") {
								$T2++;
							}
							else {
								die "VCF reference base error, line: $vcf_line\n";
							}
						}
						elsif ($calls[0] eq "1") {
							if ($alt eq "A") {
								$A2++;
							}
							elsif ($alt eq "C") {
								$C2++;
							}
							elsif ($alt eq "G") {
								$G2++;
							}
							elsif ($alt eq "T") {
								$T2++;
							}
							else {
								die "VCF alternative base error, line: $vcf_line\n";
							}	
						}	
					}
					$vcf_count++;
				}		
				print TEMP "$vcf_cols[0]\t$vcf_cols[1]\t$A1:$C1:$G1:$T1\t$A2:$C2:$G2:$T2\n";
			}
		}
	}
	close VCF;
	close TEMP;
	$table = "temp.tsv";
	print "...successfully parsed\n";
}


################ section 2: read input commands and setup output files ################

#turn off/on (0/1) different output modes (F=feature, O=overall)
my $modeF = 0;
my $modeO = 0;

#parse mode arguments 
@modes = split(/,/,join(',',@modes));
my @sorted_modes = sort { lc($a) cmp lc($b) } @modes;
my $argument_number = @sorted_modes;
if ( ($argument_number == 0) or ($argument_number > 2) ) {
	die "error in mode selection, can be overall, feature, or overall,feature\n";
}
else {
	if ( ($sorted_modes[0] eq "overall") and ($argument_number == 1) ) {
		$modeO = 1;
		print "calculating overall statistics for all sites...\n";
	}
	elsif ( ($sorted_modes[0] eq "feature") and ($argument_number == 1) ) {
		$modeF = 1;
		print "calculating statistics per feature...\n";
	}
	elsif ( ($sorted_modes[0] eq "feature") and ($sorted_modes[1] eq "overall") and ($argument_number == 2) ) {
		$modeF = 1;
		$modeO = 1;
		print "calculating statistics per feature, and overall for all sites...\n";
	}
	else {
		die "ERROR: mode not set correctly, feature or overall?\n";
	}
} 

#open output files with user specified prefix
if ( $modeO == 1 ) {
	open (OUT1, ">$prefix.overall.txt") or die; 
}
if ( $modeF == 1 ) {
	open (OUT2, ">$prefix.per_feature.txt") or die; 
}

#parse output arguments
@outputs = split(/,/,join(',',@outputs));
my $pi = 0;
my $theta = 0;
my $TajD = 0;
my $foldedSFS = 0;
my $Fst = 0;
my $dxy = 0;
foreach my $output_arguments (@outputs) {
	chomp $output_arguments;
	if ($output_arguments eq "pi") {
		$pi = 1;
	}
	elsif ($output_arguments eq "theta") {
		$theta = 1;
	}
	elsif ($output_arguments eq "TajD") {
		$TajD = 1;
	}
	elsif ($output_arguments eq "folded-SFS") {
    	$foldedSFS = 1;
	}
	elsif ($output_arguments eq "Fst") {
		if ($populations == 2) {
			$Fst = 1;
		}
		else {
			print "WARNING: to calculate Fst the populations flag must be set to 2. Fst will not be calculated...\n";
		}
	}
	elsif ($output_arguments eq "dxy") {
		if ($populations == 2) {
			$dxy = 1;
		}
		else {
			print "WARNING: to calculate dxy the populations flag must be set to 2. dxy will not be calculated...\n";
		}
	}
	else {
		print "WARNING: unknown output argument: $output_arguments...\n";
	}
}

#print appropriate headers for by feature output file
if ( $modeF == 1 ) { #if mode is set to 'feature'
	if ($format eq "BED3") {
		print OUT2 "chromosome\tstart\tend";
	}
	elsif ($format eq "BED6") {
		print OUT2 "chromosome\tstart\tend\tname\tscore\tstrand";
	}
	else {
		die "incorrect argument: bed_format, BED3 or BED6?";
	}
	if ($feature_method eq "by_line") {
		if ($populations == 1) {
			if ($pi == 1) {
				print OUT2 "\tpi";
			}
			if ($theta == 1) {
				print OUT2 "\ttheta";
			}
			if ($TajD == 1) {
				print OUT2 "\tTajima's_D";
			}
			if ($foldedSFS == 1) {
				print OUT2 "\tfolded-SFS($max_copies[0] copies)";
			}
			print OUT2 "\n";
		}
		elsif ($populations == 2) {
			foreach my $pop_name (@pop_names) {
				if ($pi == 1) {
					print OUT2 "\t$pop_name pi";
				}
				if ($theta == 1) {
					print OUT2 "\t$pop_name theta";
				}
				if ($TajD == 1) {
					print OUT2 "\t$pop_name Tajima's_D";
				}
				if ($foldedSFS == 1) {
					print OUT2 "\t$pop_name folded-SFS";
				}
			}
			if ($Fst == 1) {
				print OUT2 "\tFst";
			}
			if ($dxy == 1) {
				print OUT2 "\tdxy";
			}
			print OUT2 "\n";
		}
	}
	elsif ($feature_method eq "by_annotation") {
		print OUT2 "\tfeature_ID";
		if ($populations == 1) {
			if ($pi == 1) {
				print OUT2 "\tpi";
			}
			if ($theta == 1) {
				print OUT2 "\ttheta";
			}
			if ($TajD == 1) {
				print OUT2 "\tTajima's_D";
			}
			if ($foldedSFS == 1) {
				print OUT2 "\tfolded-SFS($max_copies[0] copies)";
			}
			print OUT2 "\n";
		}
		elsif ($populations == 2) {
			foreach my $pop_name (@pop_names) {
				if ($pi == 1) {
					print OUT2 "\t$pop_name pi";
				}
				if ($theta == 1) {
					print OUT2 "\t$pop_name theta";
				}
				if ($TajD == 1) {
					print OUT2 "\t$pop_name Tajima's_D";
				}
				if ($foldedSFS == 1) {
					print OUT2 "\t$pop_name folded-SFS";
				}
			}
			if ($Fst == 1) {
				print OUT2 "\tFst";
			}
			if ($dxy == 1) {
				print OUT2 "\tdxy";
			}
			print OUT2 "\n";
		}
	}
	else {
		die "features option error, by_line or by_annotation?";
	}
}

print "...input options successfully read\n";

################ section 3: read input file and make coordinate hashes ################

open (IN1, "$bed") or die;

print "now parsing input file...";

my %feature;
my %feature_index;
my $feature_count = 0;

while (my $line1 = <IN1>) {
	chomp $line1;
	$feature_count++;
	my @cols1 = split(/\t/, $line1);
	my $chromosome = $cols1[0];
	my $start = $cols1[1] + 1; #convert 0-based to 1-based coordinate
	my $end = $cols1[2];
	my $count = $start;
	until ($count == ($end + 1) ) {
		if ( ( $modeF == 1 ) and ($feature_method eq "by_line") ) {
			my $key_feature = "$chromosome\t$count";
			my $feature_id = "f$feature_count";
			$feature{$key_feature} = $feature_id; #add feature ID as value if mode is feature, and feature is assigned by line
			$feature_index{$feature_id} = $line1; #add original bed line value with feature_id key
		}
		if ( ( $modeF == 1 ) and ($feature_method eq "by_annotation") ) {
			my $key_feature = "$chromosome\t$count";
			$feature{$key_feature} = $cols1[3]; #add feature ID as value if mode is feature, and feature is assigned by annotation (must be in column 4 of bed file!)
		}
		elsif ($feature_method eq "none") { 
			my $key_feature = "$chromosome\t$count";
			$feature{$key_feature} = 1; #output not by feature, so null value of 1
		}
		$count++;
	}
}

close IN1;

print "successfully completed\n";

################ section 4: read annotation table and store SFS for calculating popgen statistics ################

open (IN2, "$table") or die;

my @copy_counts1;
my $loop_count1 = $min_copies[0];

#make an array of accepted number of copies at a locus
until ($loop_count1 == ($max_copies[0] + 1) ) {
	push @copy_counts1, $loop_count1;
	$loop_count1++;
}

my @copy_counts2;
my $loop_count2 = $min_copies[1];

if ($populations == 2) {
	until ($loop_count2 == ($max_copies[1] + 1) ) {
		push @copy_counts2, $loop_count2;
		$loop_count2++;
	}	
}

my @polymorphic_sites1;
my $loop_count3 = 0;

#make an array of zeros for possible number of polymorphic sites
my $allele_match1 = ( ($max_copies[0] / 2) + 1);
if ($max_copies[0] % 2 == 1) { #number is odd
	$allele_match1 = $allele_match1 - 0.5;
}
until ($loop_count3 == $allele_match1 ) {
	push @polymorphic_sites1, 0;
	$loop_count3++;
}

my @polymorphic_sites2;
my $loop_count4 = 0;

if ($populations == 2) {
	my $allele_match2 = ( ($max_copies[1] / 2) + 1);
	if ($max_copies[1] % 2 == 1) { #number is odd
		$allele_match2 = $allele_match2 - 0.5;
	}
	until ($loop_count4 == $allele_match2 ) {
		push @polymorphic_sites2, 0;
		$loop_count4++;
	}
}

#if mode is overall, make hashes to contain: number of copies -> SFS
my %SFS;
my %SFS_PPUU;
my %SFS_PU;
my %SFS_UP;
my %pop_freqs;

#make hash of arrays containing SFS for each possible number of copies
if ($modeO == 1) {
	if ($populations == 1) {
		foreach my $copy_number (@copy_counts1) {
			$SFS{$copy_number} = [@polymorphic_sites1];
			if ($codon_usage == 1) { #make additional SFS hases for each of the 3 codon categories
				$SFS_PPUU{$copy_number} = [@polymorphic_sites1];
				$SFS_PU{$copy_number} = [@polymorphic_sites1];
				$SFS_UP{$copy_number} = [@polymorphic_sites1];
			}	
		}
	}
	elsif ($populations == 2) { #make hashses with additional hash layer for each population: number of copies -> population -> SFS
		my $copy_counter = 1;
		foreach my $pop_name (@pop_names) {
			if ($copy_counter == 1) {
				foreach my $copy_number1 (@copy_counts1) {
					$SFS{$pop_name}{$copy_number1} = [@polymorphic_sites1];
					if ($codon_usage == 1) { #make additional SFS hases for each of the 3 codon categories
						$SFS_PPUU{$pop_name}{$copy_number1} = [@polymorphic_sites1];
						$SFS_PU{$pop_name}{$copy_number1} = [@polymorphic_sites1];
						$SFS_UP{$pop_name}{$copy_number1} = [@polymorphic_sites1];
					}				
				}
			}
			elsif ($copy_counter == 2) {
				foreach my $copy_number2 (@copy_counts2) {
					$SFS{$pop_name}{$copy_number2} = [@polymorphic_sites2];
					if ($codon_usage == 1) { #make additional SFS hases for each of the 3 codon categories
						$SFS_PPUU{$pop_name}{$copy_number2} = [@polymorphic_sites2];
						$SFS_PU{$pop_name}{$copy_number2} = [@polymorphic_sites2];
						$SFS_UP{$pop_name}{$copy_number2} = [@polymorphic_sites2];
					}				
				}
			}
			$copy_counter++;
		}
	}
}

#if mode is features, make hashes of hashes to contain: feature_id -> number of copies -> SFS
my %SFS_features;
my %SFS_PPUU_features;
my %SFS_PU_features;
my %SFS_UP_features;
my %pop_freqs_features;


#make hash of arrays containing SFS for each possible number of copies
if ( $modeF == 1 ) {
	foreach my $feature_ids (values %feature) { #feature IDs, either unique ID per line or per annotation
		if ($populations == 1) {
			foreach my $copy_number (@copy_counts1) {
				$SFS_features{$feature_ids}{$copy_number} = [@polymorphic_sites1];
				if ($codon_usage == 1) {
					$SFS_PPUU_features{$feature_ids}{$copy_number} = [@polymorphic_sites1];
					$SFS_PU_features{$feature_ids}{$copy_number} = [@polymorphic_sites1];
					$SFS_UP_features{$feature_ids}{$copy_number} = [@polymorphic_sites1];
				}
			}
		}
		elsif ($populations == 2) { #make hashses with additional hash layer for each population: feature_id -> population -> number of copies -> SFS
			my $copy_counter = 1;
			foreach my $pop_name (@pop_names) {
				if ($copy_counter == 1) {
					foreach my $copy_number1 (@copy_counts1) {
						$SFS_features{$feature_ids}{$pop_name}{$copy_number1} = [@polymorphic_sites1];
						if ($codon_usage == 1) {
							$SFS_PPUU_features{$feature_ids}{$pop_name}{$copy_number1} = [@polymorphic_sites1];
							$SFS_PU_features{$feature_ids}{$pop_name}{$copy_number1} = [@polymorphic_sites1];
							$SFS_UP_features{$feature_ids}{$pop_name}{$copy_number1} = [@polymorphic_sites1];
						}
					}
				}
				elsif ($copy_counter == 2) {
					foreach my $copy_number2 (@copy_counts2) {
						$SFS_features{$feature_ids}{$pop_name}{$copy_number2} = [@polymorphic_sites2];
						if ($codon_usage == 1) {
							$SFS_PPUU_features{$feature_ids}{$pop_name}{$copy_number2} = [@polymorphic_sites2];
							$SFS_PU_features{$feature_ids}{$pop_name}{$copy_number2} = [@polymorphic_sites2];
							$SFS_UP_features{$feature_ids}{$pop_name}{$copy_number2} = [@polymorphic_sites2];
						}
					}		
				}
				$copy_counter++;
			}
		}
	}
}

print "now transversing the big table...";

my $line_count = 0;

#loop through table and match sites within features
while (my $line2 = <IN2>) {
	$line_count++;
	if ($line2 !~ /^#/) {
		chomp $line2;
		my @cols2 = split(/\t/, $line2);
		my $key_match = "$cols2[0]\t$cols2[1]";
		if (exists $feature{$key_match}) { #match between bed file and table
			my $bases2 = pop @cols2; #assumes alleles are last column in line
			my @allele_counts;
			push @allele_counts, $bases2;
			if ($populations == 2) {
				my $bases1 = pop @cols2;
				push @allele_counts, $bases1; #if two populations, assumes pop1 alleles are second last column in line, and pop2 alleles are last column
			}
			my @rev_allele_counts = reverse @allele_counts;
			my @allele_freqs;
			my $pop_count = 0;
			foreach my $allele_count (@rev_allele_counts) {
				my @alleles = split(/:/, $allele_count);
				my $A_count = $alleles[0];
				my $C_count = $alleles[1];
				my $G_count = $alleles[2];
				my $T_count = $alleles[3];
				my @spectrum = ($alleles[0], $alleles[1], $alleles[2], $alleles[3]);
				my $number_seqs = $A_count + $C_count + $G_count + $T_count;
				my $number_alleles = 0;
				if ($A_count > 0) {
					$number_alleles++;
				}
				if ($C_count > 0) {
					$number_alleles++;
				}
				if ($G_count > 0) {
					$number_alleles++;
				}
				if ($T_count > 0) {
					$number_alleles++;
				}
				if ($number_alleles == 2) { #ignore sites with more than 2 alleles
					if ($number_seqs >= $min_copies[$pop_count]) { #ignore sites with insufficient number of called individuals 
						if ($number_seqs > $max_copies[$pop_count]) {
							die "allele copies exceeds declared max_copies at table line $line_count\n$number_seqs copies and max copies set to $max_copies[$pop_count]\n";
						}
						my @segregating_alleles;
						foreach my $base (@spectrum) {
							unless ($base == 0) {
								push @segregating_alleles, $base; #collect frequency of two alleles in array
							} 
						}
						my @sorted_spectrum = sort { $a <=> $b } @segregating_alleles;
						my $maf = $sorted_spectrum[0]; #sort array numerically, and take lowest number as minor allele frequency
						if ($populations == 1) {
							$SFS{$number_seqs}[$maf]++; #add one to SFS category for specific number of sequences
							if ($codon_usage == 1) { #get SFS for each of the the 3 codon categories
								if ( ( ($G_count > 0) and ($C_count > 0) ) or ( ($A_count > 0) and ($T_count > 0) ) ) { #all SNPs where either P (G/C) alleles are present together, or U (A/T)
									$SFS_PPUU{$number_seqs}[$maf]++;
								}
								if ( ( ($G_count > 0) and ($A_count > 0) and ($G_count > $A_count) ) or ( ($G_count > 0) and ($T_count > 0) and ($G_count > $T_count) ) or ( ($C_count > 0) and ($A_count > 0) and ($C_count > $A_count) ) or ( ($C_count > 0) and ($T_count > 0) and ($C_count > $T_count) ) ) { #assumes minor allele is derived 
									$SFS_PU{$number_seqs}[$maf]++;
								}
								if ( ( ($A_count > 0) and ($G_count > 0) and ($A_count > $G_count) ) or ( ($T_count > 0) and ($G_count > 0) and ($T_count > $G_count) ) or ( ($A_count > 0) and ($C_count > 0) and ($A_count > $C_count) ) or ( ($T_count > 0) and ($C_count > 0) and ($T_count > $C_count) ) ) { #assumes minor allele is derived
									$SFS_UP{$number_seqs}[$maf]++;
								}
							}
							if ( $modeF == 1 ) { #if mode is set to "feature"
								$SFS_features{$feature{$key_match}}{$number_seqs}[$maf]++; #add one to SFS category for specific number of sequences for specific feature ID
								if ($codon_usage == 1) { #get SFS for each of the the 3 codon categories
									if ( ( ($G_count > 0) and ($C_count > 0) ) or ( ($A_count > 0) and ($T_count > 0) ) ) { #all SNPs where either P (G/C) alleles are present together, or U (A/T)
										$SFS_PPUU_features{$feature{$key_match}}{$number_seqs}[$maf]++;
									}
									if ( ( ($G_count > 0) and ($A_count > 0) and ($G_count > $A_count) ) or ( ($G_count > 0) and ($T_count > 0) and ($G_count > $T_count) ) or ( ($C_count > 0) and ($A_count > 0) and ($C_count > $A_count) ) or ( ($C_count > 0) and ($T_count > 0) and ($C_count > $T_count) ) ) { #assumes minor allele is derived 
										$SFS_PU_features{$feature{$key_match}}{$number_seqs}[$maf]++;
									}
									if ( ( ($A_count > 0) and ($G_count > 0) and ($A_count > $G_count) ) or ( ($T_count > 0) and ($G_count > 0) and ($T_count > $G_count) ) or ( ($A_count > 0) and ($C_count > 0) and ($A_count > $C_count) ) or ( ($T_count > 0) and ($C_count > 0) and ($T_count > $C_count) ) ) { #assumes minor allele is derived
										$SFS_UP_features{$feature{$key_match}}{$number_seqs}[$maf]++;
									}
								}
							}
						}
						elsif ($populations == 2) {
							$SFS{$pop_names[$pop_count]}{$number_seqs}[$maf]++; #add one to SFS category for specific number of sequences
							if ($codon_usage == 1) { #get SFS for each of the the 3 codon categories
								if ( ( ($G_count > 0) and ($C_count > 0) ) or ( ($A_count > 0) and ($T_count > 0) ) ) { #all SNPs where either P (G/C) alleles are present together, or U (A/T)
									$SFS_PPUU{$pop_names[$pop_count]}{$number_seqs}[$maf]++;
								}
								if ( ( ($G_count > 0) and ($A_count > 0) and ($G_count > $A_count) ) or ( ($G_count > 0) and ($T_count > 0) and ($G_count > $T_count) ) or ( ($C_count > 0) and ($A_count > 0) and ($C_count > $A_count) ) or ( ($C_count > 0) and ($T_count > 0) and ($C_count > $T_count) ) ) { #assumes minor allele is derived 
									$SFS_PU{$pop_names[$pop_count]}{$number_seqs}[$maf]++;
								}
								if ( ( ($A_count > 0) and ($G_count > 0) and ($A_count > $G_count) ) or ( ($T_count > 0) and ($G_count > 0) and ($T_count > $G_count) ) or ( ($A_count > 0) and ($C_count > 0) and ($A_count > $C_count) ) or ( ($T_count > 0) and ($C_count > 0) and ($T_count > $C_count) ) ) { #assumes minor allele is derived
									$SFS_UP{$pop_names[$pop_count]}{$number_seqs}[$maf]++;
								}
							}
							if ( $modeF == 1 ) { #if mode is set to "feature"
								$SFS_features{$feature{$key_match}}{$pop_names[$pop_count]}{$number_seqs}[$maf]++; #add one to SFS category for specific number of sequences for specific feature ID
								if ($codon_usage == 1) { #get SFS for each of the the 3 codon categories
									if ( ( ($G_count > 0) and ($C_count > 0) ) or ( ($A_count > 0) and ($T_count > 0) ) ) { #all SNPs where either P (G/C) alleles are present together, or U (A/T)
										$SFS_PPUU_features{$feature{$key_match}}{$pop_names[$pop_count]}{$number_seqs}[$maf]++;
									}
									if ( ( ($G_count > 0) and ($A_count > 0) and ($G_count > $A_count) ) or ( ($G_count > 0) and ($T_count > 0) and ($G_count > $T_count) ) or ( ($C_count > 0) and ($A_count > 0) and ($C_count > $A_count) ) or ( ($C_count > 0) and ($T_count > 0) and ($C_count > $T_count) ) ) { #assumes minor allele is derived 
										$SFS_PU_features{$feature{$key_match}}{$pop_names[$pop_count]}{$number_seqs}[$maf]++;
									}
									if ( ( ($A_count > 0) and ($G_count > 0) and ($A_count > $G_count) ) or ( ($T_count > 0) and ($G_count > 0) and ($T_count > $G_count) ) or ( ($A_count > 0) and ($C_count > 0) and ($A_count > $C_count) ) or ( ($T_count > 0) and ($C_count > 0) and ($T_count > $C_count) ) ) { #assumes minor allele is derived
										$SFS_UP_features{$feature{$key_match}}{$pop_names[$pop_count]}{$number_seqs}[$maf]++;
									}
								}
							}
						}
					}
				}	
				elsif ($number_alleles == 1) { #add invariant sites
					if ($number_seqs >= $min_copies[$pop_count]) {
						if ($populations == 1) {
							$SFS{$number_seqs}[0]++; #add one to SFS for invariant sites for specific number of sequences 
							if ($modeF == 1) { #if mode is set to "feature"
								$SFS_features{$feature{$key_match}}{$number_seqs}[0]++; #add one to invariant sites for specific number of sequences for specific feature ID
							}
							if ($codon_usage == 1) {
								$SFS_PPUU{$number_seqs}[0]++;
								if ( ($G_count > 0) or ($C_count > 0) ) {
									$SFS_PU{$number_seqs}[0]++;
								}
								if ( ($A_count > 0) or ($T_count > 0) ) {
									$SFS_UP{$number_seqs}[0]++;
								}            
							}
						}
						elsif ($populations == 2) {
							$SFS{$pop_names[$pop_count]}{$number_seqs}[0]++; #add one to SFS for invariant sites for specific number of sequences 
							if ($modeF == 1) { #if mode is set to "feature"
								$SFS_features{$feature{$key_match}}{$pop_names[$pop_count]}{$number_seqs}[0]++; #add one to invariant sites for specific number of sequences for specific feature ID
							}
							if ($codon_usage == 1) {
								$SFS_PPUU{$pop_names[$pop_count]}{$number_seqs}[0]++;
								if ( ($G_count > 0) or ($C_count > 0) ) {
									$SFS_PU{$pop_names[$pop_count]}{$number_seqs}[0]++;
								}
								if ( ($A_count > 0) or ($T_count > 0) ) {
									$SFS_UP{$pop_names[$pop_count]}{$number_seqs}[0]++;
								}            
							}
						}

					}
				}
				$pop_count++;
			}
			if ( ($populations == 2) and ( ($Fst == 1) or ($dxy == 1) ) ) {
				my @alleles1 = split(/:/, $rev_allele_counts[0]);
				my @alleles2 = split(/:/, $rev_allele_counts[1]);
				my $A_count1 = $alleles1[0];
				my $C_count1 = $alleles1[1];
				my $G_count1 = $alleles1[2];
				my $T_count1 = $alleles1[3];
				my $A_count2 = $alleles2[0];
				my $C_count2 = $alleles2[1];
				my $G_count2 = $alleles2[2];
				my $T_count2 = $alleles2[3];
				my $number_seqs1 = $A_count1 + $C_count1 + $G_count1 + $T_count1;
				my $number_seqs2 = $A_count2 + $C_count2 + $G_count2 + $T_count2;
				my $total_count = $number_seqs1 + $number_seqs2;
				if ( ($number_seqs1 >= $min_copies[0]) and ($number_seqs2 >= $min_copies[1])) {
					if ( ($A_count1 == $number_seqs1) or ($C_count1 == $number_seqs1) or ($G_count1 == $number_seqs1) or ($T_count1 == $number_seqs1) ) { #pop 1 is invariant
						if ( ($A_count2 == $number_seqs2) or ($C_count2 == $number_seqs2) or ($G_count2 == $number_seqs2) or ($T_count2 == $number_seqs2) ) { #pop 1 is invariant	
							if ( (($A_count1 + $A_count2) == $total_count) or (($C_count1 + $C_count2) == $total_count) or (($G_count1 + $G_count2) == $total_count) or (($T_count1 + $T_count2) == $total_count) ) { #site is invariant overall
								$pop_freqs{"$number_seqs1 $number_seqs2"}{"1 0"}{"1 0"}++; #add one to frequency 1,0 - 1,0 for a given number of sequences combination
								if ($modeF == 1) { #if mode is set to "feature"
									$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"1 0"}{"1 0"}++; #add one to frequency 1,0 - 1,0 for given feature and given number of sequences combination
								}
							}
							else { #site is fixed between pops
								$pop_freqs{"$number_seqs1 $number_seqs2"}{"1 0"}{"0 1"}++; #add one to frequency 1,0 for a given number of sequences combination
								if ($modeF == 1) { #if mode is set to "feature"
									$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"1 0"}{"0 1"}++; #add one to frequency 1,0 - 0,1 for given feature and given number of sequences combination
								}
							}
						}
						else { #site is invariant in pop 1 but variant in pop 2
							if ($A_count1 > 0) { #A is fixed in pop 1
								if ( ($A_count2 > 0) and ($C_count2 > 0) ) {
									my $Afreq2 = $A_count2 / $number_seqs2;
									my $Cfreq2 = $C_count2 / $number_seqs2;
									$pop_freqs{"$number_seqs1 $number_seqs2"}{"1 0"}{"$Afreq2 $Cfreq2"}++;
									if ($modeF == 1) { #if mode is set to "feature"
										$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"1 0"}{"$Afreq2 $Cfreq2"}++; 
									}
								}
								elsif ( ($A_count2 > 0) and ($G_count2 > 0) ) {
									my $Afreq2 = $A_count2 / $number_seqs2;
									my $Gfreq2 = $G_count2 / $number_seqs2;
									$pop_freqs{"$number_seqs1 $number_seqs2"}{"1 0"}{"$Afreq2 $Gfreq2"}++;
									if ($modeF == 1) { #if mode is set to "feature"
										$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"1 0"}{"$Afreq2 $Gfreq2"}++; 
									}
								}
								elsif ( ($A_count2 > 0) and ($T_count2 > 0) ) {
									my $Afreq2 = $A_count2 / $number_seqs2;
									my $Tfreq2 = $T_count2 / $number_seqs2;
									$pop_freqs{"$number_seqs1 $number_seqs2"}{"1 0"}{"$Afreq2 $Tfreq2"}++;
									if ($modeF == 1) { #if mode is set to "feature"
										$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"1 0"}{"$Afreq2 $Tfreq2"}++; 
									}
								}
							}						
							elsif ($C_count1 > 0) { #C is fixed in pop 1
								if ( ($C_count2 > 0) and ($A_count2 > 0) ) {
									my $Cfreq2 = $C_count2 / $number_seqs2;
									my $Afreq2 = $A_count2 / $number_seqs2;
									$pop_freqs{"$number_seqs1 $number_seqs2"}{"1 0"}{"$Cfreq2 $Afreq2"}++;
									if ($modeF == 1) { #if mode is set to "feature"
										$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"1 0"}{"$Cfreq2 $Afreq2"}++; 
									}
								}
								elsif ( ($C_count2 > 0) and ($G_count2 > 0) ) {
									my $Cfreq2 = $C_count2 / $number_seqs2;
									my $Gfreq2 = $G_count2 / $number_seqs2;
									$pop_freqs{"$number_seqs1 $number_seqs2"}{"1 0"}{"$Cfreq2 $Gfreq2"}++;
									if ($modeF == 1) { #if mode is set to "feature"
										$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"1 0"}{"$Cfreq2 $Gfreq2"}++; 
									}
								}						
								elsif ( ($C_count2 > 0) and ($T_count2 > 0) ) {
									my $Cfreq2 = $C_count2 / $number_seqs2;
									my $Tfreq2 = $T_count2 / $number_seqs2;
									$pop_freqs{"$number_seqs1 $number_seqs2"}{"1 0"}{"$Cfreq2 $Tfreq2"}++;
									if ($modeF == 1) { #if mode is set to "feature"
										$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"1 0"}{"$Cfreq2 $Tfreq2"}++; 
									}
								}	
							}
							elsif ($G_count1 > 0) { #G is fixed in pop 1
								if ( ($G_count2 > 0) and ($A_count2 > 0) ) {
									my $Gfreq2 = $G_count2 / $number_seqs2;
									my $Afreq2 = $A_count2 / $number_seqs2;
									$pop_freqs{"$number_seqs1 $number_seqs2"}{"1 0"}{"$Gfreq2 $Afreq2"}++;
									if ($modeF == 1) { #if mode is set to "feature"
										$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"1 0"}{"$Gfreq2 $Afreq2"}++; 
									}
								}
								elsif ( ($G_count2 > 0) and ($C_count2 > 0) ) {
									my $Gfreq2 = $G_count2 / $number_seqs2;
									my $Cfreq2 = $C_count2 / $number_seqs2;
									$pop_freqs{"$number_seqs1 $number_seqs2"}{"1 0"}{"$Gfreq2 $Cfreq2"}++;
									if ($modeF == 1) { #if mode is set to "feature"
										$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"1 0"}{"$Gfreq2 $Cfreq2"}++; 
									}
								}
								elsif ( ($G_count2 > 0) and ($T_count2 > 0) ) {
									my $Gfreq2 = $G_count2 / $number_seqs2;
									my $Tfreq2 = $T_count2 / $number_seqs2;
									$pop_freqs{"$number_seqs1 $number_seqs2"}{"1 0"}{"$Gfreq2 $Tfreq2"}++;
									if ($modeF == 1) { #if mode is set to "feature"
										$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"1 0"}{"$Gfreq2 $Tfreq2"}++; 
									}
								}
							}
							elsif ($T_count1 > 0) {
								if ( ($T_count2 > 0) and ($A_count2 > 0) ) {
									my $Tfreq2 = $T_count2 / $number_seqs2;
									my $Afreq2 = $A_count2 / $number_seqs2;
									$pop_freqs{"$number_seqs1 $number_seqs2"}{"1 0"}{"$Tfreq2 $Afreq2"}++;
									if ($modeF == 1) { #if mode is set to "feature"
										$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"1 0"}{"$Tfreq2 $Afreq2"}++; 
									}	
								}
								elsif ( ($T_count2 > 0) and ($C_count2 > 0) ) {
									my $Tfreq2 = $T_count2 / $number_seqs2;
									my $Cfreq2 = $C_count2 / $number_seqs2;
									$pop_freqs{"$number_seqs1 $number_seqs2"}{"1 0"}{"$Tfreq2 $Cfreq2"}++;
									if ($modeF == 1) { #if mode is set to "feature"
										$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"1 0"}{"$Tfreq2 $Cfreq2"}++; 
									}	
								}
								elsif ( ($T_count2 > 0) and ($G_count2 > 0) ) {
									my $Tfreq2 = $T_count2 / $number_seqs2;
									my $Gfreq2 = $G_count2 / $number_seqs2;
									$pop_freqs{"$number_seqs1 $number_seqs2"}{"1 0"}{"$Tfreq2 $Gfreq2"}++;
									if ($modeF == 1) { #if mode is set to "feature"
										$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"1 0"}{"$Tfreq2 $Gfreq2"}++; 
									}	
								}
							}
						}
					}
					elsif ( ($A_count2 == $number_seqs2) or ($C_count2 == $number_seqs2) or ($G_count2 == $number_seqs2) or ($T_count2 == $number_seqs2) ) { #pop2 is invariant, pop1 is variant
						if ($A_count2 > 0) { #A is fixed in pop 2
							if ( ($A_count1 > 0) and ($C_count1 > 0) ) {
								my $Afreq1 = $A_count1 / $number_seqs1;
								my $Cfreq1 = $C_count1 / $number_seqs1;
								$pop_freqs{"$number_seqs1 $number_seqs2"}{"$Afreq1 $Cfreq1"}{"1 0"}++;
								if ($modeF == 1) { #if mode is set to "feature"
									$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"$Afreq1 $Cfreq1"}{"1 0"}++; 
								}
							}
							elsif ( ($A_count1 > 0) and ($G_count1 > 0) ) {
								my $Afreq1 = $A_count1 / $number_seqs1;
								my $Gfreq1 = $G_count1 / $number_seqs1;
								$pop_freqs{"$number_seqs1 $number_seqs2"}{"$Afreq1 $Gfreq1"}{"1 0"}++;
								if ($modeF == 1) { #if mode is set to "feature"
									$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"$Afreq1 $Gfreq1"}{"1 0"}++; 
								}
							}
							elsif ( ($A_count1 > 0) and ($T_count1 > 0) ) {
								my $Afreq1 = $A_count1 / $number_seqs1;
								my $Tfreq1 = $T_count1 / $number_seqs1;
								$pop_freqs{"$number_seqs1 $number_seqs2"}{"$Afreq1 $Tfreq1"}{"1 0"}++;
								if ($modeF == 1) { #if mode is set to "feature"
									$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"$Afreq1 $Tfreq1"}{"1 0"}++; 
								}
							}
						}			
						elsif ($C_count2 > 0) { #C is fixed in pop 2
							if ( ($C_count1 > 0) and ($A_count1 > 0) ) {
								my $Cfreq1 = $C_count1 / $number_seqs1;
								my $Afreq1 = $A_count1 / $number_seqs1;
								$pop_freqs{"$number_seqs1 $number_seqs2"}{"$Cfreq1 $Afreq1"}{"1 0"}++;
								if ($modeF == 1) { #if mode is set to "feature"
									$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"$Cfreq1 $Afreq1"}{"1 0"}++; 
								}
							}
							elsif ( ($C_count1 > 0) and ($G_count1 > 0) ) {
								my $Cfreq1 = $C_count1 / $number_seqs1;
								my $Gfreq1 = $G_count1 / $number_seqs1;
								$pop_freqs{"$number_seqs1 $number_seqs2"}{"$Cfreq1 $Gfreq1"}{"1 0"}++;
								if ($modeF == 1) { #if mode is set to "feature"
									$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"$Cfreq1 $Gfreq1"}{"1 0"}++; 
								}
							}						
							elsif ( ($C_count1 > 0) and ($T_count1 > 0) ) {
								my $Cfreq1 = $C_count1 / $number_seqs1;
								my $Tfreq1 = $T_count1 / $number_seqs1;
								$pop_freqs{"$number_seqs1 $number_seqs2"}{"$Cfreq1 $Tfreq1"}{"1 0"}++;
								if ($modeF == 1) { #if mode is set to "feature"
									$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"$Cfreq1 $Tfreq1"}{"1 0"}++; 
								}
							}	
						}
						elsif ($G_count2 > 0) { #G is fixed in pop 2
							if ( ($G_count1 > 0) and ($A_count1 > 0) ) {
								my $Gfreq1 = $G_count1 / $number_seqs1;
								my $Afreq1 = $A_count1 / $number_seqs1;
								$pop_freqs{"$number_seqs1 $number_seqs2"}{"$Gfreq1 $Afreq1"}{"1 0"}++;
								if ($modeF == 1) { #if mode is set to "feature"
									$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"$Gfreq1 $Afreq1"}{"1 0"}++; 
								}
							}
							elsif ( ($G_count1 > 0) and ($C_count1 > 0) ) {
								my $Gfreq1 = $G_count1 / $number_seqs1;
								my $Cfreq1 = $C_count1 / $number_seqs1;
								$pop_freqs{"$number_seqs1 $number_seqs2"}{"$Gfreq1 $Cfreq1"}{"1 0"}++;
								if ($modeF == 1) { #if mode is set to "feature"
									$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"$Gfreq1 $Cfreq1"}{"1 0"}++; 
								}
							}
							elsif ( ($G_count1 > 0) and ($T_count1 > 0) ) {
								my $Gfreq1 = $G_count1 / $number_seqs1;
								my $Tfreq1 = $T_count1 / $number_seqs1;
								$pop_freqs{"$number_seqs1 $number_seqs2"}{"$Gfreq1 $Tfreq1"}{"1 0"}++;
								if ($modeF == 1) { #if mode is set to "feature"
									$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"$Gfreq1 $Tfreq1"}{"1 0"}++; 
								}
							}
						}
						elsif ($T_count2 > 0) {
							if ( ($T_count1 > 0) and ($A_count1 > 0) ) {
								my $Tfreq1 = $T_count1 / $number_seqs1;
								my $Afreq1 = $A_count1 / $number_seqs1;
								$pop_freqs{"$number_seqs1 $number_seqs2"}{"$Tfreq1 $Afreq1"}{"1 0"}++;
								if ($modeF == 1) { #if mode is set to "feature"
									$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"$Tfreq1 $Afreq1"}{"1 0"}++; 
								}	
							}
							elsif ( ($T_count1 > 0) and ($C_count1 > 0) ) {
								my $Tfreq1 = $T_count1 / $number_seqs1;
								my $Cfreq1 = $C_count1 / $number_seqs1;
								$pop_freqs{"$number_seqs1 $number_seqs2"}{"$Tfreq1 $Cfreq1"}{"1 0"}++;
								if ($modeF == 1) { #if mode is set to "feature"
									$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"$Tfreq1 $Cfreq1"}{"1 0"}++; 
								}	
							}
							elsif ( ($T_count1 > 0) and ($G_count1 > 0) ) {
								my $Tfreq1 = $T_count1 / $number_seqs1;
								my $Gfreq1 = $G_count1 / $number_seqs1;
								$pop_freqs{"$number_seqs1 $number_seqs2"}{"$Tfreq1 $Gfreq1"}{"1 0"}++;
								if ($modeF == 1) { #if mode is set to "feature"
									$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"$Tfreq1 $Gfreq1"}{"1 0"}++; 
								}	
							}
						}
					}
					elsif ( ($A_count1 != $number_seqs1) and ($C_count1 != $number_seqs1) and ($G_count1 != $number_seqs1) and ($T_count1 != $number_seqs1) ) { #both pops are variant, assumes bi-allelic SNPs between pops!	
						if ( ( ($A_count1 > 0) and ($C_count1 > 0) ) and  ( ($A_count2 > 0) and ($C_count2 > 0) ) ) {
							my $Afreq1 = $A_count1 / $number_seqs1;
							my $Cfreq1 = $C_count1 / $number_seqs1;
							my $Afreq2 = $A_count2 / $number_seqs2;
							my $Cfreq2 = $C_count2 / $number_seqs2;
							$pop_freqs{"$number_seqs1 $number_seqs2"}{"$Afreq1 $Cfreq1"}{"$Afreq2 $Cfreq2"}++;
							if ($modeF == 1) { #if mode is set to "feature"
								$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"$Afreq1 $Cfreq1"}{"$Afreq2 $Cfreq2"}++; 
							}
						}
						elsif ( ( ($A_count1 > 0) and ($G_count1 > 0) ) and  ( ($A_count2 > 0) and ($G_count2 > 0) ) ) {
							my $Afreq1 = $A_count1 / $number_seqs1;
							my $Gfreq1 = $G_count1 / $number_seqs1;
							my $Afreq2 = $A_count2 / $number_seqs2;
							my $Gfreq2 = $G_count2 / $number_seqs2;
							$pop_freqs{"$number_seqs1 $number_seqs2"}{"$Afreq1 $Gfreq1"}{"$Afreq2 $Gfreq2"}++;
							if ($modeF == 1) { #if mode is set to "feature"
								$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"$Afreq1 $Gfreq1"}{"$Afreq2 $Gfreq2"}++; 
							}
						}
						elsif ( ( ($A_count1 > 0) and ($T_count1 > 0) ) and  ( ($A_count2 > 0) and ($T_count2 > 0) ) ) {
							my $Afreq1 = $A_count1 / $number_seqs1;
							my $Tfreq1 = $T_count1 / $number_seqs1;
							my $Afreq2 = $A_count2 / $number_seqs2;
							my $Tfreq2 = $T_count2 / $number_seqs2;
							$pop_freqs{"$number_seqs1 $number_seqs2"}{"$Afreq1 $Tfreq1"}{"$Afreq2 $Tfreq2"}++;
							if ($modeF == 1) { #if mode is set to "feature"
								$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"$Afreq1 $Tfreq1"}{"$Afreq2 $Tfreq2"}++; 
							}
						}
						elsif ( ( ($C_count1 > 0) and ($G_count1 > 0) ) and  ( ($C_count2 > 0) and ($G_count2 > 0) ) ) {
							my $Cfreq1 = $C_count1 / $number_seqs1;
							my $Gfreq1 = $G_count1 / $number_seqs1;
							my $Cfreq2 = $C_count2 / $number_seqs2;
							my $Gfreq2 = $G_count2 / $number_seqs2;
							$pop_freqs{"$number_seqs1 $number_seqs2"}{"$Cfreq1 $Gfreq1"}{"$Cfreq2 $Gfreq2"}++;
							if ($modeF == 1) { #if mode is set to "feature"
								$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"$Cfreq1 $Gfreq1"}{"$Cfreq2 $Gfreq2"}++; 
							}
						}
						elsif ( ( ($C_count1 > 0) and ($T_count1 > 0) ) and  ( ($C_count2 > 0) and ($T_count2 > 0) ) ) {
							my $Cfreq1 = $C_count1 / $number_seqs1;
							my $Tfreq1 = $T_count1 / $number_seqs1;
							my $Cfreq2 = $C_count2 / $number_seqs2;
							my $Tfreq2 = $T_count2 / $number_seqs2;
							$pop_freqs{"$number_seqs1 $number_seqs2"}{"$Cfreq1 $Tfreq1"}{"$Cfreq2 $Tfreq2"}++;
							if ($modeF == 1) { #if mode is set to "feature"
								$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"$Cfreq1 $Tfreq1"}{"$Cfreq2 $Tfreq2"}++; 
							}
						}
						elsif ( ( ($G_count1 > 0) and ($T_count1 > 0) ) and  ( ($G_count2 > 0) and ($T_count2 > 0) ) ) {
							my $Gfreq1 = $G_count1 / $number_seqs1;
							my $Tfreq1 = $T_count1 / $number_seqs1;
							my $Gfreq2 = $G_count2 / $number_seqs2;
							my $Tfreq2 = $T_count2 / $number_seqs2;
							$pop_freqs{"$number_seqs1 $number_seqs2"}{"$Gfreq1 $Tfreq1"}{"$Gfreq2 $Tfreq2"}++;
							if ($modeF == 1) { #if mode is set to "feature"
								$pop_freqs_features{$feature{$key_match}}{"$number_seqs1 $number_seqs2"}{"$Gfreq1 $Tfreq1"}{"$Gfreq2 $Tfreq2"}++; 
							}
						}
						else {
							$Fst = 0;
							$dxy = 0;
							print "WARNING: error in assignment of allele frequencies between populations, Fst and dxy will not be calculated. Are all SNPs bi-allelic even between populations?\n";
						}
					}
				}
			}
		}
	}
}

close IN2;

if ($input_format eq "vcf") {
	system("rm temp.tsv");
}

print "finished!\n";
print "calculating popgen statistics...";


################ section 5: calculating popgen statistics ################

if ($modeO == 1) { #calculate overall statistics
	if ( ( ($pi == 1) or ($theta == 1) or ($TajD == 1) ) and ($populations == 1) ) {
		my %pi_hash;
		my %theta_hash;
		my %D_hash;
		my $total_number_sites_all_copies = 0;
		my $total_pi_all_copies = 0;
		my $invalid_sites = 0;
		foreach my $key (sort { $a <=> $b } keys %SFS) { #keys are different numbers of copies, values contain SFS for each
			my ($denominator, $numerator, $number_of_sites, $number_of_variable_sites) = pi ($key, \@{$SFS{$key}});
			$total_number_sites_all_copies += $number_of_sites;
			if ($numerator == 0) {
				if ($pi == 1) {
					if ($number_of_sites > 0) {
						print OUT1 "pi for $key copies is 0 and number of sites is $number_of_sites\n";
					}
					elsif ($number_of_sites == 0) {
						print OUT1 "pi for $key copies is NA, number of valid sites is $number_of_sites\n";
					}	
				}
				if ($theta == 1) {
					if ($number_of_sites > 0) {
						print OUT1 "theta for $key copies is 0 and number of sites is $number_of_sites\n";
					}
					elsif ($number_of_sites == 0) {
						print OUT1 "theta for $key copies is NA, number of valid sites is $number_of_sites\n";
					}	
				}
				if ($TajD == 1) {
					if ($number_of_sites > 0) {
						print OUT1 "Tajima's D for $key copies is NA and number of sites is $number_of_sites\n";
					}
					elsif ($number_of_sites == 0) {
						print OUT1 "Tajima's D for $key copies is NA, number of valid sites is $number_of_sites\n";
					}
				}
				$pi_hash{$key} = "0\t$number_of_sites"; #if no variable sites store 0 for pi, no variable sites
				$theta_hash{$key} = "0\t$number_of_sites"; #if no variable sites store 0 for theta
				$invalid_sites += $number_of_sites;
			}
			else {
				my $total_pi = $numerator / $denominator;
				$total_pi_all_copies += $total_pi;
				my $per_site_pi = ($numerator / $denominator) / $number_of_sites;
				$pi_hash{$key} = "$per_site_pi\t$number_of_sites";
				if ($pi == 1) {
					print OUT1 "pi for $key copies is $per_site_pi and number of sites is $number_of_sites\n";
				}
				my $total_theta = theta ( $key, $number_of_variable_sites );
				my $per_site_theta = $total_theta / $number_of_sites;
				$theta_hash{$key} = "$per_site_theta\t$number_of_sites";
				if ($theta == 1) {
					print OUT1 "theta for $key copies is $per_site_theta and number of sites is $number_of_sites\n";
				}
				my $D = tajima( $key, $total_pi, $number_of_variable_sites );
				unless ($D eq "NA") {
					$D_hash{$key} = "$D\t$number_of_sites";
				}
				if ($TajD == 1) {
					print OUT1 "Tajima's D for $key copies is $D and number of sites is $number_of_sites\n";
				}
			}
		}
		#calculate overall pi
		if ( ($pi_weighted_average == 1) and ($pi == 1 ) ) { #calculate overall pi as a weighted average
			my $sum_of_weights_pi = 0;
			my $weighted_average_pi = 0;
			foreach my $pis (keys %pi_hash) { #keys are copy numbers
				my @values = split(/\t/, $pi_hash{$pis});
				$sum_of_weights_pi += $values[0] * $values[1]; #per site pi * number of sites
			}
			if ($sum_of_weights_pi > 0) {
				$weighted_average_pi = $sum_of_weights_pi / $total_number_sites_all_copies;
				print OUT1 "overall weighted average pi is: $weighted_average_pi\ntotal number of sites is: $total_number_sites_all_copies\n";
			}
			elsif ($total_number_sites_all_copies > 0) {
		    	print OUT1 "overall weighted average pi is: 0\ttotal number of sites is: $total_number_sites_all_copies\n";
			}
			elsif ($total_number_sites_all_copies == 0) {
				print OUT1 "no sites with sufficient number of allele copies for regions of interest, pi cannot be calculated\n";
			} 
		}
		elsif ($pi == 1) { #calculate overall pi as sum of all pi values / total number of sites for all copies
			if ($total_pi_all_copies > 0) {
				my $overall_per_site_pi = $total_pi_all_copies / $total_number_sites_all_copies;
				print OUT1 "overall pi is: $overall_per_site_pi\ntotal number of sites is: $total_number_sites_all_copies\n";
			}
			elsif ($total_number_sites_all_copies > 0) {
		    	print OUT1 "overall pi is: 0\ttotal number of sites is: $total_number_sites_all_copies\n";
			}
			elsif ($total_number_sites_all_copies == 0) {
				print OUT1 "no sites with sufficient number of allele copies for regions of interest, pi cannot be calculated\n";
			} 
		}
		#take weighted average for overall theta
		if ($theta == 1) {
			my $sum_of_weights_theta = 0;
			my $weighted_average_theta = 0;
			foreach my $thetas (keys %theta_hash) { #keys are copy numbers
				my @values = split(/\t/, $theta_hash{$thetas});
				$sum_of_weights_theta += $values[0] * $values[1]; #per site theta * number of sites
			}
			if ($sum_of_weights_theta > 0) {
				$weighted_average_theta = $sum_of_weights_theta / $total_number_sites_all_copies;
				print OUT1 "overall weighted average theta is: $weighted_average_theta\ntotal number of sites is: $total_number_sites_all_copies\n";
			}
			elsif ($total_number_sites_all_copies > 0) {
		    	print OUT1 "overall weighted average theta is: 0\ntotal number of sites is: $total_number_sites_all_copies\n";
			}
			elsif ($total_number_sites_all_copies == 0) {
				print OUT1 "no sites with sufficient number of alleles for regions of interest, theta cannot be calculated\n";
			}
		}
		#take weighted average for overall Tajima's D
		if ($TajD == 1) {
			my $sum_of_weights_D = 0;
			my $weighted_average_D = 0;
			my $valid_sites = $total_number_sites_all_copies - $invalid_sites;
			unless (!keys %D_hash) {
				foreach my $Ds (keys %D_hash) { #keys are copy numbers
					my @values = split(/\t/, $D_hash{$Ds});
					$sum_of_weights_D += $values[0] * $values[1]; #Tajima'S D * number of sites
				}
				if ( ($sum_of_weights_D > 0) or ($sum_of_weights_D < 0) ) {
					$weighted_average_D = $sum_of_weights_D / ($valid_sites);
					print OUT1 "overall weighted average Tajima's D is: $weighted_average_D\ntotal number of sites is: $valid_sites\n";
				}
				elsif ($total_number_sites_all_copies > 0) {
			    	print OUT1 "overall weighted average Tajima's D is: NA\ntotal number of sites is: $valid_sites\n";
				}
				elsif ($total_number_sites_all_copies == 0) {
					print OUT1 "no sites with sufficient number of alleles for regions of interest, Tajima's D cannot be calculated\n";
				}
			}
			else {
				print OUT1 "overall weighted average Tajima's D is: NA\ntotal number of sites is: $valid_sites\n";
			}
		}
	}
	#output SFS for max allele copies
	if ( ($foldedSFS == 1) and ($populations == 1) ) {    
		print OUT1 "\nfolded-SFS for $max_copies allele copies is: @{$SFS{$max_copies[0]}}\n";
		if ($codon_usage == 1) {
			print OUT1 "\noverall SFS for codon categories is:\n";
			print OUT1 "P->P/U->U: @{$SFS_PPUU{$max_copies[0]}}\n";
			print OUT1 "P->U: @{$SFS_PU{$max_copies[0]}}\n";
			print OUT1 "U->P: @{$SFS_UP{$max_copies[0]}}\n";
		}
	}
	if ($populations == 2) {
		my $pop_count = 0;
		my @pop_pis;
		foreach my $pop_name (@pop_names) {
			if ( ($pi == 1) or ($theta == 1) or ($TajD == 1) or ($Fst == 1) ) {
				my %pi_hash;
				my %theta_hash;
				my %D_hash;
				my $total_number_sites_all_copies = 0;
				my $total_pi_all_copies = 0;
				my $invalid_sites = 0;
				foreach my $sub_key (sort { $a <=> $b } keys %{ $SFS{$pop_name} } ) {
					my ($denominator, $numerator, $number_of_sites, $number_of_variable_sites) = pi ($sub_key, \@{$SFS{$pop_name}{$sub_key}});
					$total_number_sites_all_copies += $number_of_sites;	
					if ($numerator == 0) {
						if ($pi == 1) {
							if ($number_of_sites > 0) {
								print OUT1 "pi for $pop_name and $sub_key copies is 0 and number of sites is $number_of_sites\n";
							}
							elsif ($number_of_sites == 0) {
								print OUT1 "pi for $pop_name and $sub_key copies is NA, number of valid sites is $number_of_sites\n";
							}	
						}
						if ($theta == 1) {
							if ($number_of_sites > 0) {
								print OUT1 "theta for $pop_name and $sub_key copies is 0 and number of sites is $number_of_sites\n";
							}
							elsif ($number_of_sites == 0) {
								print OUT1 "theta for $pop_name and $sub_key copies is NA, number of valid sites is $number_of_sites\n";
							}	
						}
						if ($TajD == 1) {
							if ($number_of_sites > 0) {
								print OUT1 "Tajima's D for $pop_name and $sub_key copies is NA and number of sites is $number_of_sites\n";
							}
							elsif ($number_of_sites == 0) {
								print OUT1 "Tajima's D for $pop_name and $sub_key copies is NA, number of valid sites is $number_of_sites\n";
							}
						}
						$pi_hash{$sub_key} = "0\t$number_of_sites"; #if no variable sites store 0 for pi, no variable sites
						$theta_hash{$sub_key} = "0\t$number_of_sites"; #if no variable sites store 0 for theta
						$invalid_sites += $number_of_sites;
					}
					else {
						my $total_pi = $numerator / $denominator;
						$total_pi_all_copies += $total_pi;
						my $per_site_pi = ($numerator / $denominator) / $number_of_sites;
						$pi_hash{$sub_key} = "$per_site_pi\t$number_of_sites";
						if ($pi == 1) {
							print OUT1 "pi for $pop_name and $sub_key copies is $per_site_pi and number of sites is $number_of_sites\n";
						}
						my $total_theta = theta ( $sub_key, $number_of_variable_sites );
						my $per_site_theta = $total_theta / $number_of_sites;
						$theta_hash{$sub_key} = "$per_site_theta\t$number_of_sites";
						if ($theta == 1) {
							print OUT1 "theta for $pop_name and $sub_key copies is $per_site_theta and number of sites is $number_of_sites\n";
						}
						my $D = tajima( $sub_key, $total_pi, $number_of_variable_sites );
						unless ($D eq "NA") {
							$D_hash{$sub_key} = "$D\t$number_of_sites";
						}
						if ($TajD == 1) {
							print OUT1 "Tajima's D for $pop_name and $sub_key copies is $D and number of sites is $number_of_sites\n";
						}
					}
				}
				#calculate overall pi
				if ( ($pi_weighted_average == 1) and ($pi == 1 ) ) { #calculate overall pi as a weighted average
					my $sum_of_weights_pi = 0;
					my $weighted_average_pi = 0;
					foreach my $pis (keys %pi_hash) { #keys are copy numbers
						my @values = split(/\t/, $pi_hash{$pis});
						$sum_of_weights_pi += $values[0] * $values[1]; #per site pi * number of sites
					}
					if ($sum_of_weights_pi > 0) {
						$weighted_average_pi = $sum_of_weights_pi / $total_number_sites_all_copies;
						print OUT1 "overall weighted average pi for $pop_name is: $weighted_average_pi and total number of sites is: $total_number_sites_all_copies\n";
						push @pop_pis, $weighted_average_pi;
					}
					elsif ($total_number_sites_all_copies > 0) {
				    	print OUT1 "overall weighted average pi for $pop_name is: 0 and total number of sites is: $total_number_sites_all_copies\n";
						push @pop_pis, 0;
					}
					elsif ($total_number_sites_all_copies == 0) {
						print OUT1 "no sites with sufficient number of allele copies for regions of interest, pi for $pop_name cannot be calculated\n";
						push @pop_pis, "NA";
					} 
				}
				elsif ($pi == 1) { #calculate overall pi as sum of all pi values / total number of sites for all copies
					if ($total_pi_all_copies > 0) {
						my $overall_per_site_pi = $total_pi_all_copies / $total_number_sites_all_copies;
						print OUT1 "overall pi for $pop_name is: $overall_per_site_pi and total number of sites is: $total_number_sites_all_copies\n";
						push @pop_pis, $overall_per_site_pi;
					}
					elsif ($total_number_sites_all_copies > 0) {
				    	print OUT1 "overall pi for $pop_name is: 0 and total number of sites is: $total_number_sites_all_copies\n";
						push @pop_pis, 0;
					}
					elsif ($total_number_sites_all_copies == 0) {
						print OUT1 "no sites with sufficient number of allele copies for regions of interest, pi for $pop_name cannot be calculated\n";
						push @pop_pis, "NA";
					} 
				}
				#take weighted average for overall theta
				if ($theta == 1) {
					my $sum_of_weights_theta = 0;
					my $weighted_average_theta = 0;
					foreach my $thetas (keys %theta_hash) { #keys are copy numbers
						my @values = split(/\t/, $theta_hash{$thetas});
						$sum_of_weights_theta += $values[0] * $values[1]; #per site theta * number of sites
					}
					if ($sum_of_weights_theta > 0) {
						$weighted_average_theta = $sum_of_weights_theta / $total_number_sites_all_copies;
						print OUT1 "overall weighted average theta for $pop_name is: $weighted_average_theta and total number of sites is: $total_number_sites_all_copies\n";
					}
					elsif ($total_number_sites_all_copies > 0) {
				    	print OUT1 "overall weighted average theta for $pop_name is: 0 and total number of sites is: $total_number_sites_all_copies\n";
					}
					elsif ($total_number_sites_all_copies == 0) {
						print OUT1 "no sites with sufficient number of alleles for regions of interest, theta for $pop_name cannot be calculated\n";
					}
				}
				#take weighted average for overall Tajima's D
				my $sum_of_weights_D = 0;
				my $weighted_average_D = 0;
				my $valid_sites = $total_number_sites_all_copies - $invalid_sites;
				if ($TajD == 1) {
					unless (!keys %D_hash) {
						foreach my $Ds (keys %D_hash) { #keys are copy numbers
							my @values = split(/\t/, $D_hash{$Ds});
							$sum_of_weights_D += $values[0] * $values[1]; #Tajima'S D * number of sites
						}
						if ( ($sum_of_weights_D > 0) or ($sum_of_weights_D < 0) ) {
							$weighted_average_D = $sum_of_weights_D / $valid_sites;
							print OUT1 "overall weighted average Tajima's D for $pop_name is: $weighted_average_D and total number of sites is: $valid_sites\n";
						}
						elsif ($total_number_sites_all_copies > 0) {
					    	print OUT1 "overall weighted average Tajima's D for $pop_name is: NA and total number of sites is: $valid_sites\n";
						}
						elsif ($total_number_sites_all_copies == 0) {
							print OUT1 "no sites with sufficient number of alleles for regions of interest, Tajima's D for $pop_name cannot be calculated\n";
						}
					}
					else {
						print OUT1 "overall weighted average Tajima's D for $pop_name is: NA and total number of sites is: $valid_sites\n";
					}	
				}
			}
			if ($foldedSFS == 1) {
				print OUT1 "\nfolded-SFS for $pop_name and $max_copies[$pop_count] allele copies is: @{$SFS{$pop_name}{$max_copies[$pop_count]}}\n";
				if ($codon_usage == 1) {
					print OUT1 "\noverall SFS for codon categories is:\n";
					print OUT1 "P->P/U->U: @{$SFS_PPUU{$pop_name}{$max_copies[$pop_count]}}\n";
					print OUT1 "P->U: @{$SFS_PU{$pop_name}{$max_copies[$pop_count]}}\n";
					print OUT1 "U->P: @{$SFS_UP{$pop_name}{$max_copies[$pop_count]}}\n";
				}
			}
			$pop_count++;
			print OUT1 "\n";
		}
		if ( ($Fst == 1) or ($dxy == 1) ) {
			my $pi_within;
			if ( ($pop_pis[0] eq "NA") or ($pop_pis[1] eq "NA") ) {
				$pi_within = "NA";
			} 
			else {
				$pi_within = ($pop_pis[0] + $pop_pis[1]) / 2;
			}
			my $total_pi_between_all_copies = 0;
			my $total_number_sites_all_copies = 0;
			my %pi_between_hash;
			foreach my $copies_key (keys %pop_freqs) {
				my $total_pi_between = 0;
				my $number_of_sites = 0;
				my @pop_copies = split (" ", $copies_key);
				my $pop1_copies = $pop_copies[0];
				my $pop2_copies = $pop_copies[1];
				foreach my $pop1_key (keys %{ $pop_freqs{$copies_key} } ) {
					my @pop1_freqs = split (" ", $pop1_key);
					my $p1 = $pop1_freqs[0];
					my $q1 = $pop1_freqs[1];
					foreach my $pop2_key (keys %{ $pop_freqs{$copies_key}{$pop1_key} } ) {
						my @pop2_freqs = split (" ", $pop2_key);
						my $p2 = $pop2_freqs[0];
						my $q2 = $pop2_freqs[1];
						$total_pi_between += ( ($p1 * $q2) + ($p2 * $q1) ) * $pop_freqs{$copies_key}{$pop1_key}{$pop2_key};
						$total_pi_between_all_copies += ( ($p1 * $q2) + ($p2 * $q1) ) * $pop_freqs{$copies_key}{$pop1_key}{$pop2_key};
						$number_of_sites += $pop_freqs{$copies_key}{$pop1_key}{$pop2_key};
						$total_number_sites_all_copies += $pop_freqs{$copies_key}{$pop1_key}{$pop2_key};
					}
				} 
				my $per_site_pi_between = $total_pi_between / $number_of_sites;
				$pi_between_hash{$copies_key} = "$per_site_pi_between\t$number_of_sites";
			}
			if ($pi_weighted_average == 1) {
				print OUT1 "\n";
				my $sum_of_weights_pi_between = 0;
				my $weighted_average_pi_between = 0;
				foreach my $pis_between (keys %pi_between_hash) {
					my @values = split(/\t/, $pi_between_hash{$pis_between});
					$sum_of_weights_pi_between += $values[0] * $values[1]; #per site pi * number of sites
				}
				if ($sum_of_weights_pi_between > 0) {
					$weighted_average_pi_between = $sum_of_weights_pi_between / $total_number_sites_all_copies;
					if ( ($Fst == 1) and ($pi_within ne "NA") ) {
						my $Fst_overall = 1 - ($pi_within / $weighted_average_pi_between);
						print OUT1 "Fst calculated from weighted average pi within and between populations is $Fst_overall and total number of sites is: $total_number_sites_all_copies\n";
					}
					elsif ( ($Fst == 1) and ($pi_within eq "NA") ) {
						print OUT1 "Fst is NA and total number of sites is: $total_number_sites_all_copies\n";
					}
					if ($dxy == 1) {
						print OUT1 "dxy calculated as a weighted average is $weighted_average_pi_between and total number of sites is: $total_number_sites_all_copies\n";
					}
				}
				elsif ($total_number_sites_all_copies > 0) {
					if ($Fst == 1) {
						print OUT1 "pi between populations is 0, Fst is NA and total number of sites is: $total_number_sites_all_copies\n";
					}
					if ($dxy == 1) {
						print OUT1 "dxy calculated as a weighted average is 0 and total number of sites is: $total_number_sites_all_copies\n";
					}
				}
				elsif ($total_number_sites_all_copies == 0) {
					if ($Fst == 1) {
						print OUT1 "pi between populations is NA, Fst is NA and total number of sites is 0\n";
					}
					if ($dxy == 1) {
						print OUT1 "dxy is NA and total number of sites is 0\n";
					}			
				}
			}
			else {
				print OUT1 "\n";
				if ($total_pi_between_all_copies > 0) {
					my $overall_per_site_pi_between = $total_pi_between_all_copies / $total_number_sites_all_copies;
					if ( ($Fst == 1) and ($pi_within ne "NA") ) {
						my $Fst_overall = 1 - ($pi_within / $overall_per_site_pi_between);
						print OUT1 "Fst calculated from pi within and between populations is $Fst_overall and total number of sites is: $total_number_sites_all_copies\n";
					}
					elsif ( ($Fst == 1) and ($pi_within eq "NA") ) {
						print OUT1 "Fst is NA and total number of sites is: $total_number_sites_all_copies\n";
					}
					if ($dxy == 1) {
						print OUT1 "dxy is $overall_per_site_pi_between and total number of sites is: $total_number_sites_all_copies\n";
					}
				}
				elsif ($total_number_sites_all_copies > 0) {
					if ($Fst == 1) {
						print OUT1 "pi between populations is 0, Fst is NA and total number of sites is: $total_number_sites_all_copies\n";
					}
					if ($dxy == 1) {
						print OUT1 "dxy is 0 and total number of sites is: $total_number_sites_all_copies\n";
					}
				}
				elsif ($total_number_sites_all_copies == 0) {
					if ($Fst == 1) {
						print OUT1 "pi between populations is NA, Fst is NA and total number of sites is 0\n";
					}
					if ($dxy == 1) {
						print OUT1 "dxy is NA and total number of sites is 0\n";
					}			
				}
			}
		}
	}
	close OUT1;
}

if ($modeF == 1) { #calculate overall statistics
	foreach my $feature_key (sort keys %SFS_features) { #loop through top dimension of hash - feature IDs
		if ($feature_method eq "by_line") {
			print OUT2 "$feature_index{$feature_key}"; #recall line coordinates from hash of feature IDs and feature coordinates
		}
		elsif ($feature_method eq "by_annotation") {
			print OUT2 "$feature_key"; #print out feature ID
		}
		if ( ( ($pi == 1) or ($theta == 1) or ($TajD == 1) ) and ($populations == 1) ) {
			my %pi_hash_features;
			my %theta_hash_features;
			my %D_hash_features;
			my $total_number_sites_all_copies_features = 0;
			my $total_pi_all_copies_features = 0;
			my $invalid_sites_features = 0;
			foreach my $feature_sub_key (sort { $a <=> $b } keys %{ $SFS_features{$feature_key} } ) { #loop through second dimension of hash - number of copies
				my ($denominator, $numerator, $number_of_sites, $number_of_variable_sites) = pi ($feature_sub_key, \@{$SFS_features{$feature_key}{$feature_sub_key}});
				$total_number_sites_all_copies_features += $number_of_sites;
				if ($numerator == 0) {
					$pi_hash_features{$feature_sub_key} = "0\t$number_of_sites"; #if undefined store 0 for pi
					$theta_hash_features{$feature_sub_key} = "0\t$number_of_sites"; #if undefined store 0 for theta
					$invalid_sites_features += $number_of_sites;
				}
				else {
					my $total_pi_features = $numerator / $denominator;
					$total_pi_all_copies_features += $total_pi_features;
					my $per_site_pi_features = ($numerator / $denominator) / $number_of_sites;
					$pi_hash_features{$feature_sub_key} = "$per_site_pi_features\t$number_of_sites";
					my $total_theta_features = theta ( $feature_sub_key, $number_of_variable_sites);
					my $per_site_theta_features = $total_theta_features / $number_of_sites;
					$theta_hash_features{$feature_sub_key} = "$per_site_theta_features\t$number_of_sites";
					my $D_features = tajima( $feature_sub_key, $total_pi_features, $number_of_variable_sites );
					unless ($D_features eq "NA") {
						$D_hash_features{$feature_sub_key} = "$D_features\t$number_of_sites";
					}
				}
			}
			if ( ($pi_weighted_average == 1) and ($pi == 1) ) {
				my $sum_of_weights_pi_features = 0;
				my $weighted_average_pi_features = 0;
				foreach my $pis_features (keys %pi_hash_features) { #keys are copy numbers
					my @values_features = split(/\t/, $pi_hash_features{$pis_features});
					$sum_of_weights_pi_features += $values_features[0] * $values_features[1]; #per site pi * number of sites
				}
				if ($sum_of_weights_pi_features > 0) {
					$weighted_average_pi_features = $sum_of_weights_pi_features / $total_number_sites_all_copies_features;
					print OUT2 "\t$weighted_average_pi_features $total_number_sites_all_copies_features";
				}
				elsif ($total_number_sites_all_copies_features > 0) {
					print OUT2 "\t0 $total_number_sites_all_copies_features";
				}
				elsif ($total_number_sites_all_copies_features == 0) {
					print OUT2 "\tNA 0";
				} 
			}
			elsif ($pi == 1) { #calculate overall pi as sum of all pi values / total number of sites for all copies
				if ($total_pi_all_copies_features > 0) {
					my $overall_per_site_pi_features = $total_pi_all_copies_features / $total_number_sites_all_copies_features;
					print OUT2 "\t$overall_per_site_pi_features";
				}
				elsif ($total_number_sites_all_copies_features > 0) {
			    	print OUT2 "\t0";
				}
				elsif ($total_number_sites_all_copies_features == 0) {
					print OUT2 "\tNA";
				} 
			}
			if ($theta == 1) {
				my $sum_of_weights_theta_features = 0;
				my $weighted_average_theta_features = 0;
				foreach my $thetas_features (keys %theta_hash_features) { #keys are copy numbers
					my @values_features = split(/\t/, $theta_hash_features{$thetas_features});
					$sum_of_weights_theta_features += $values_features[0] * $values_features[1]; #per site theta * number of sites
				}
				if ($sum_of_weights_theta_features > 0) {
					$weighted_average_theta_features = $sum_of_weights_theta_features / $total_number_sites_all_copies_features;
					print OUT2 "\t$weighted_average_theta_features";
				}
				elsif ($total_number_sites_all_copies_features > 0) {
					print OUT2 "\t0";
				}
				elsif ($total_number_sites_all_copies_features == 0) {
					print OUT2 "\tNA";
				}
			}
			if ($TajD == 1) {
				my $sum_of_weights_D_features = 0;
				my $weighted_average_D_features = 0;
				my $valid_sites_features = $total_number_sites_all_copies_features - $invalid_sites_features;
				unless (!keys %D_hash_features) {
					foreach my $Ds_features (keys %D_hash_features) { #keys are copy numbers
						my @values_features = split(/\t/, $D_hash_features{$Ds_features});
						$sum_of_weights_D_features += $values_features[0] * $values_features[1]; #Tajima'S D * number of sites
					}
					if ( ($sum_of_weights_D_features > 0) or ($sum_of_weights_D_features < 0) ) {
						$weighted_average_D_features = $sum_of_weights_D_features / $total_number_sites_all_copies_features;
						print OUT2 "\t$weighted_average_D_features";
					}
					elsif ($total_number_sites_all_copies_features > 0) {
						print OUT2 "\t0";
					}
					elsif ($total_number_sites_all_copies_features == 0) {
						print OUT2 "\tNA";
					}
				}
				else {
					print OUT2 "\tNA";
				}
			}
		}
		if ( ($foldedSFS == 1) and ($populations == 1) ) { #print SFS for max copy number
			print OUT2 "\t@{$SFS_features{$feature_key}{$max_copies[0]} }"
		}
		if ($populations == 2) {
			my $pop_count_features = 0;
			my @pop_pis_features;
			foreach my $pop_name (@pop_names) {
				if ( ($pi == 1) or ($theta == 1) or ($TajD == 1) or ($Fst == 1) ) {
					my %pi_hash_features;
					my %theta_hash_features;
					my %D_hash_features;
					my $total_number_sites_all_copies_features = 0;
					my $total_pi_all_copies_features = 0;
					my $invalid_sites_features = 0;	
					foreach my $feature_sub_key (sort { $a <=> $b } keys %{ $SFS_features{$feature_key}{$pop_name} } ) {
						my ($denominator, $numerator, $number_of_sites, $number_of_variable_sites) = pi ($feature_sub_key, \@{$SFS_features{$feature_key}{$pop_name}{$feature_sub_key}});
						$total_number_sites_all_copies_features += $number_of_sites;
						if ($numerator == 0) {
							$pi_hash_features{$feature_sub_key} = "0\t$number_of_sites"; #if undefined store 0 for pi
							$theta_hash_features{$feature_sub_key} = "0\t$number_of_sites"; #if undefined store 0 for theta
							$invalid_sites_features += $number_of_sites;
						}
						else {
							my $total_pi_features = $numerator / $denominator;
							$total_pi_all_copies_features += $total_pi_features;
							my $per_site_pi_features = ($numerator / $denominator) / $number_of_sites;
							$pi_hash_features{$feature_sub_key} = "$per_site_pi_features\t$number_of_sites";
							my $total_theta_features = theta ( $feature_sub_key, $number_of_variable_sites);
							my $per_site_theta_features = $total_theta_features / $number_of_sites;
							$theta_hash_features{$feature_sub_key} = "$per_site_theta_features\t$number_of_sites";
							my $D_features = tajima( $feature_sub_key, $total_pi_features, $number_of_variable_sites );
							unless ($D_features eq "NA") {
								$D_hash_features{$feature_sub_key} = "$D_features\t$number_of_sites";
							}
						}
					}
					if ( ($pi_weighted_average == 1) and ($pi == 1) ) {
						my $sum_of_weights_pi_features = 0;
						my $weighted_average_pi_features = 0;
						foreach my $pis_features (keys %pi_hash_features) { #keys are copy numbers
							my @values_features = split(/\t/, $pi_hash_features{$pis_features});
							$sum_of_weights_pi_features += $values_features[0] * $values_features[1]; #per site pi * number of sites
						}
						if ($sum_of_weights_pi_features > 0) {
							$weighted_average_pi_features = $sum_of_weights_pi_features / $total_number_sites_all_copies_features;
							print OUT2 "\t$weighted_average_pi_features";
							push @pop_pis_features, $weighted_average_pi_features;
						}
						elsif ($total_number_sites_all_copies_features > 0) {
							print OUT2 "\t0";
							push @pop_pis_features, 0;
						}
						elsif ($total_number_sites_all_copies_features == 0) {
							print OUT2 "\tNA";
							push @pop_pis_features, "NA";
						} 
					}
					elsif ($pi == 1) { #calculate overall pi as sum of all pi values / total number of sites for all copies
						if ($total_pi_all_copies_features > 0) {
							my $overall_per_site_pi_features = $total_pi_all_copies_features / $total_number_sites_all_copies_features;
							print OUT2 "\t$overall_per_site_pi_features";
							push @pop_pis_features, $overall_per_site_pi_features;
						}
						elsif ($total_number_sites_all_copies_features > 0) {
					    	print OUT2 "\t0";
					    	push @pop_pis_features, 0;
						}
						elsif ($total_number_sites_all_copies_features == 0) {
							print OUT2 "\tNA";
							push @pop_pis_features, "NA";
						} 
					}
					if ($theta == 1) {
						my $sum_of_weights_theta_features = 0;
						my $weighted_average_theta_features = 0;
						foreach my $thetas_features (keys %theta_hash_features) { #keys are copy numbers
							my @values_features = split(/\t/, $theta_hash_features{$thetas_features});
							$sum_of_weights_theta_features += $values_features[0] * $values_features[1]; #per site theta * number of sites
						}
						if ($sum_of_weights_theta_features > 0) {
							$weighted_average_theta_features = $sum_of_weights_theta_features / $total_number_sites_all_copies_features;
							print OUT2 "\t$weighted_average_theta_features";
						}
						elsif ($total_number_sites_all_copies_features > 0) {
							print OUT2 "\t0";
						}
						elsif ($total_number_sites_all_copies_features == 0) {
							print OUT2 "\tNA";
						}
					}
					if ($TajD == 1) {
						my $sum_of_weights_D_features = 0;
						my $weighted_average_D_features = 0;
						my $valid_sites_features = $total_number_sites_all_copies_features - $invalid_sites_features;
						unless (!keys %D_hash_features) {
							foreach my $Ds_features (keys %D_hash_features) { #keys are copy numbers
								my @values_features = split(/\t/, $D_hash_features{$Ds_features});
								$sum_of_weights_D_features += $values_features[0] * $values_features[1]; #Tajima'S D * number of sites
							}
							if ( ($sum_of_weights_D_features > 0) or ($sum_of_weights_D_features < 0) ) {
								$weighted_average_D_features = $sum_of_weights_D_features / $total_number_sites_all_copies_features;
								print OUT2 "\t$weighted_average_D_features";
							}
							elsif ($total_number_sites_all_copies_features > 0) {
								print OUT2 "\t0";
							}
							elsif ($total_number_sites_all_copies_features == 0) {
								print OUT2 "\tNA";
							}
						}
						else {
							print OUT2 "\tNA";
						}
					}
					if ($foldedSFS == 1) {
						print OUT2 "\t@{$SFS_features{$feature_key}{$pop_name}{$max_copies[$pop_count_features]} }";
					}
				}
				$pop_count_features++;
			}
			if ( ($Fst == 1) or ($dxy == 1) ) {
				my $pi_within_features;
				if ( ($pop_pis_features[0] eq "NA") or ($pop_pis_features[1] eq "NA") ) {
					$pi_within_features = "NA";
				} 
				else {
					$pi_within_features = ($pop_pis_features[0] + $pop_pis_features[1]) / 2;
				}
				my $total_pi_between_all_copies_features = 0;
				my $total_number_sites_all_copies_features = 0;
				my %pi_between_hash_features;
				foreach my $copies_key_features (keys %{ $pop_freqs_features{$feature_key} } ) {
					my $total_pi_between_features = 0;
					my $number_of_sites_features = 0;
					my @pop_copies_features = split (" ", $copies_key_features);
					my $pop1_copies_features = $pop_copies_features[0];
					my $pop2_copies_features = $pop_copies_features[1];
					foreach my $pop1_key_features (keys %{ $pop_freqs_features{$feature_key}{$copies_key_features} } ) {
						my @pop1_freqs_features = split (" ", $pop1_key_features);
						my $p1 = $pop1_freqs_features[0];
						my $q1 = $pop1_freqs_features[1];
						foreach my $pop2_key_features (keys %{ $pop_freqs_features{$feature_key}{$copies_key_features}{$pop1_key_features} } ) {
							my @pop2_freqs_features = split (" ", $pop2_key_features);
							my $p2 = $pop2_freqs_features[0];
							my $q2 = $pop2_freqs_features[1];
							$total_pi_between_features += ( ($p1 * $q2) + ($p2 * $q1) ) * $pop_freqs_features{$feature_key}{$copies_key_features}{$pop1_key_features}{$pop2_key_features};
							$total_pi_between_all_copies_features += ( ($p1 * $q2) + ($p2 * $q1) ) * $pop_freqs_features{$feature_key}{$copies_key_features}{$pop1_key_features}{$pop2_key_features};
							$number_of_sites_features += $pop_freqs_features{$feature_key}{$copies_key_features}{$pop1_key_features}{$pop2_key_features};
							$total_number_sites_all_copies_features += $pop_freqs_features{$feature_key}{$copies_key_features}{$pop1_key_features}{$pop2_key_features};
						}
					} 
					my $per_site_pi_between_features = $total_pi_between_features / $number_of_sites_features;
					$pi_between_hash_features{$copies_key_features} = "$per_site_pi_between_features\t$number_of_sites_features";
				}
				if ($pi_weighted_average == 1) {
					my $sum_of_weights_pi_between_features = 0;
					my $weighted_average_pi_between_features = 0;
					foreach my $pis_between_features (keys %pi_between_hash_features) {
						my @values_features = split(/\t/, $pi_between_hash_features{$pis_between_features});
						$sum_of_weights_pi_between_features += $values_features[0] * $values_features[1]; #per site pi * number of sites
					}
					if ($sum_of_weights_pi_between_features > 0) {
						$weighted_average_pi_between_features = $sum_of_weights_pi_between_features / $total_number_sites_all_copies_features;
						if ( ($Fst == 1) and ($pi_within_features ne "NA") ) {
							my $Fst_overall = 1 - ($pi_within_features / $weighted_average_pi_between_features);
							print OUT2 "\t$Fst_overall";
						}
						elsif ( ($Fst == 1) and ($pi_within_features eq "NA") ) {
							print OUT2 "\tNA";
						}
						if ($dxy == 1) {
							print OUT2 "\t$weighted_average_pi_between_features";
						}
					}
					elsif ($total_number_sites_all_copies_features > 0) {
						if ($Fst == 1) {
							print OUT2 "\tNA";
						}
						if ($dxy == 1) {
							print OUT2 "\t0";
						}
					}
					elsif ($total_number_sites_all_copies_features == 0) {
						if ($Fst == 1) {
							print OUT2 "\tNA";
						}
						if ($dxy == 1) {
							print OUT2 "\tNA";
						}			
					}
				}
				else {
					if ($total_pi_between_all_copies_features > 0) {
						my $overall_per_site_pi_between_features = $total_pi_between_all_copies_features / $total_number_sites_all_copies_features;
						if ( ($Fst == 1) and ($pi_within_features ne "NA") ) {
							my $Fst_overall = 1 - ($pi_within_features / $overall_per_site_pi_between_features);
							print OUT2 "\t$Fst_overall";
						}
						elsif ( ($Fst == 1) and ($pi_within_features eq "NA") ) {
							print OUT2 "\tNA";
						}
						if ($dxy == 1) {
							print OUT2 "\t$overall_per_site_pi_between_features";
						}
					}
					elsif ($total_number_sites_all_copies_features > 0) {
						if ($Fst == 1) {
							print OUT2 "\tNA";
						}
						if ($dxy == 1) {
							print OUT2 "\t0";
						}
					}
					elsif ($total_number_sites_all_copies_features == 0) {
						if ($Fst == 1) {
							print OUT2 "\tNA";
						}
						if ($dxy == 1) {
							print OUT2 "\tNA";
						}			
					}
				}
			}
		}
		print OUT2 "\n";
	}
	close OUT2;
}

print "all done!\n";

exit;

################ subroutines ################

sub pi {
	my $k = $_[0];
	my @sfs = @{$_[1]};
	my $denominator = ( $k * ($k - 1)) / 2; #in pi equation: k(k-1)/2
	my $numerator = 0;
	my $number_of_sites = 0;
	my $number_of_variable_sites = 0;
	my $sfs_count = 0;
	foreach my $sfs_elem (@sfs) {
		$number_of_sites += $sfs_elem;
		if ( ($sfs_count != 0) and ($sfs_elem != 0) ) { #do not consider invariant sites or nonexistant polymorphism counts
			$number_of_variable_sites += $sfs_elem;
			$numerator += ($sfs_count * ($k - $sfs_count) * $sfs_elem); #sfs_count is singleton (==1), doubleton(==2), tripleton (==3) etc.
		}
		$sfs_count++;
	}
	return ($denominator, $numerator, $number_of_sites, $number_of_variable_sites);
}

sub theta {
	my $n = $_[0];
	my $S = $_[1];
	my $a1;
	my $i = 1;
	until ($i == $n) {
		$a1 += (1 / $i);
		$i++;		
	}
	my $theta = $S / $a1;
	return $theta;
}

sub tajima {
	my $n = $_[0];
	my $pi = $_[1];
	my $S = $_[2];
	my $a1;
	my $a2;
	my $i = 1;
	until ($i == $n) {
		$a1 += (1 / $i);
		$a2 += (1 / ($i * $i) );
		$i++;		
	}
	my $b1 = ($n + 1) / (3 * ($n - 1) );
	my $b2 = (2 * ( ($n * $n) + $n + 3) ) / (9 * $n * ($n - 1) );
	my $c1 = $b1 - (1 / $a1);
	my $c2 = $b2 - ( ($n + 2) / ($a1 * $n) ) + ($a2 / ($a1 * $a1) );
	my $e1 = $c1 / $a1;
	my $e2 = $c2 / ( ($a1 * $a1) + $a2 );
	if ( ( ($pi - ($S / $a1) ) ) == 0) {
		return "NA";
	}
	else {
		my $D = ($pi - ($S / $a1) ) / sqrt ( ($e1 * $S) + ($e2 * $S * ($S - 1) ) );
		return $D;		
	}
}


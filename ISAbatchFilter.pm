#!/usr/bin/perl

# ISAbatchFilter.pm
# Dr. Xiaolin Wu & Daria W. Wells

# This module contains the subroutines that process the demultiplexed and pair-matched 
# ISA data and generates fasta files for alignment.

package ISAbatchFilter; 

use strict;
use warnings;

use Exporter qw(import);
use IO::Zlib;

our @EXPORT_OK = qw(batch_filter_1_internal batch_filter_2_length20 batch_filter_3_q20 make_unique_pair_fasta);




		
# If an internal sequence is defined, this subroutine reads the copied original data file and
# remove the reads that contain that sequence immediately after the LTR junction.
	
sub batch_filter_1_internal {
	
	# Initialize variables
	my ($filter_seq,$LTR, $data_file,$folder)= @_;
	my $filter_match = my $others = my $all_reads = 0;
	
	print "Filtering by sequence $filter_seq...\n";
	
	# Open the file handle to read the data. First check if the data file is a .gz
	# file and open appropriately.
	my $data_in; # Initialize the import file handle
	if ($data_file =~ m/\.gz$/) {
		open($data_in, "gzip -dc $data_file |") or die "Unable to open the data file in filter 1 $!\n";
	}else {
		open $data_in, ">", "$folder/PEread_LTR_NonInternal_linker.txt" or die "Unable to open the data file in filter 1 $!\n";
	}
	open my $data_out, ">", "$folder/PEread_LTR_NonInternal_linker_long.txt" or die "Unable print to new file in filter 1: internal. $!\n";
	
		
	# Look at each line of the file. Print the reads that don't contain the internal
	# sequence to a new file. 	
	while (my $newline = <$data_in>) {
		$all_reads++; # Count all pairs in the demultiplex results
		my ($post_LTR) = (split /\t/, $newline)[3]; # Column 3 contains the potential post-junction host sequence
		if ($post_LTR =~/^$filter_seq/) { # Check the post-junction target sequence
			$filter_match++; # Count the number of times the target sequence is found
		}
		else {
			print $data_out $newline; # Print the lines that don't contain the target sequence
			$others++; # Count the lines that don't match the target sequence
		}
	}
	close $data_in;
	close $data_out;
	
	# Report the number of internal and non-internal reads to the log.
	open my $log_fh, ">>", "$folder/blat_log.txt" or die "Unable to print to log file in batch filter 1: internal. $!";
	print $log_fh "Filtering read pairs by sequence\nPairs analyzed:\t$all_reads\n\tFiltered pairs:\t$filter_match\n\tPassed pairs:\t$others\n\n";
	close $log_fh;	
	
	# Delete interim report.
	unlink "$folder/PEread_LTR_Linker.txt";	

	return;
}

##########################################################################################################################

# This subroutine filters out the reads where the post-LTR sequence is shorter than 20 bases.

sub batch_filter_2_length20 {

	# Initialize variables
	my ($data_file,$LTR,$folder) = @_;	
	my $all_reads = my $ing2short = my $ing2long = 0;
	
	# Obtain the LTR junction sequence from the primer file generated in the demultiplex step.
	# Change it to the reverse complement for matching to read 2.
	my $jct_file = "$LTR/$LTR"."_primer.txt";
	open my $jct_in, "<", $jct_file or die "$!: $jct_file";
	<$jct_in>;
	my $junction = <$jct_in>;
	close $jct_in;
	$junction = reverse $junction;
	$junction =~ tr/ATGC/TACG/;
	
	my $data_in;
	
	
	# Open the file handle to read the data. First check if the data file is a .gz
	# file and open appropriately.
	if ($data_file =~ m/\.gz$/) {
		open($data_in, "gzip -dc $data_file |") or die "Unable to open the data file in filter 2: length >= 20 bp. $data_file $!\n";
	}else {
		open $data_in, "<", "$folder/$data_file" or die "Unable to open the data file in filter 2: length >= 20 bp. $!\n";
	}
	open my $data_out, ">", "$folder/PEread_LTR_NonInternal_linker_long.txt" or die "Unable print to new file in filter 2: length >= 20 bp. $! $data_file\n";
	
	print "Removing post-LTR sequences shorter than 20 bases...\n";
	
	
	# Test each for the length of the post LTR sequence. Print the lines to a new file if
	# the post LTR sequence is 20 bases or longer.
	while (my $newline = <$data_in>) {
		$all_reads++; # Count all of the reads imported.
		# Obtain the post-LTR sequence (read 1) and the post-linker sequence (read 2)
		my ($R1insert,$R2insert) = (split/\t/, $newline)[3,8];
		
		# First remove any reverse complement linker sequence from read 1 and any reverse
		# complement LTR sequence from read 2.
		$R1insert =~ s/(.*)AGTCCTCTAA.*/$1/; # First remove any reverse complement linker sequence from read 1
		$R2insert =~ s/(.*)$junction.*/$1/; # Remove any reverse complement LTR sequence from read 2

		# Only accept pairs where both the post-LTR sequence and the post-Linker sequence are at least 20 bp.
		if(length $R1insert<20 or length $R2insert<20) { # If the sequence is too short, skip it and track it.
		
			$ing2short++; # Count the number of rejected pairs.
			
		}else { # If the sequence is long enough, print it and track it.
		
			print $data_out $newline; # Print the lines that are at least 20 bases
			$ing2long++; # Count the long enough read pairs
		}
	}
	close $data_in;
	close $data_out;
	
	# Report the number of too short and long enough reads to the log.
	open my $log_fh, ">>", "$folder/blat_log.txt" or die "Unable to print to log file in batch filter 2: length >= 20. $!";
	print $log_fh "Filtering read pairs by length\nPairs analyzed:\t$all_reads\n\tFiltered pairs <20bp:\t$ing2short\n\tPassed pairs >=20bp:\t$ing2long\n\n";
	close $log_fh;

	# Delete interim report.
	unlink "$folder/$data_file";
	return;
}

##########################################################################################################################

# This subroutine will filter out the low quality reads, defined as reads whose q20
# score is less than 20.

sub batch_filter_3_q20 {
	
	my ($LTR,$folder) = @_;
	
	open my $data_in, "<$folder/PEread_LTR_NonInternal_linker_long.txt" or die "Unable to open the data file in filter 3: q20 > 20. $!\n";
	open my $data_out, ">$folder/PEread_LTR_NonInternal_linker_long_q20.txt" or die "Unable print to new file in filter 3: q20 > 20. $!\n";

	print "Removing low quality sequences...\n";

	# Initialize variables
	my $qcutoff = 20;
	my $totalin = my $totalout = my $filtered = 0;

	# Iterate through the data. 
	while (my $newline = <$data_in>) {
	
		# Compare each read separately
		my $sumL=0;
		my $sumR=0;
		
		$totalin++; # Count the number of reads imported
		my @elesnew = split(/\t/, $newline); # Split the input line into its component parts
		my $LTR = $elesnew[2]; # Obtain the HIV LTR sequence from read 1
		my $insertL = $elesnew[3];  # Obtain the post-LTR sequence from read 1
		my $qscoreL = $elesnew[4];  # Obtain the quality score for the post-LTR sequence
		
		# Calculate the quality score of the first 20 bp of the post-LTR sequence
		my $LTRlen = length($LTR); # Obtain the length of the HIV LTR sequence
		# Obtain the quality score of the first 20 post-LTR sequence by skipping bases equal to the HIV LTR length
		my $insertqscore20bpL = substr $qscoreL, $LTRlen, 20;
		
		# Convert the read1 quality score to an array
		my @qscoresL = split(//,$insertqscore20bpL);
		# Convert each ascii to numerical value and calculate average, A=65
		foreach(@qscoresL) {
			my $qL=ord($_);
			$sumL += $qL;
		}
		my $averageqscore20bpL = $sumL/20-33;  #start from ASCII 33 in Illumina new version
		# Source: wikipedia Fastq_format, older version -64
		
		# Get the quality score of post-linker sequence
		my $insertR = $elesnew[8]; # Obtain the post-linker sequence from read 2
		# Obtain the quality score of the first 20 post-linker sequence by skipping bases equal to the linker length
		my $qscoreR = substr $elesnew[9], 35, 20;
		# Convert the read2 quality score to an array
		my @qscoresR = split(//, $qscoreR);
		# Convert each ascii to numerical value and calculate average, A=65
		foreach(@qscoresR) {
			my $qR = ord($_);
			$sumR += $qR;
		}
		my $averagescore20bpR = $sumR/20-33; #start from ASCII 33 in Illumina new version

		# Determine if both reads have acceptable Q scores		
		if($averageqscore20bpL>$qcutoff and $averagescore20bpR>$qcutoff) {
			$totalout++; # Count the passing reads
			print $data_out $newline; # Print the passing reads
		}else{$filtered++;} # Count the removed reads
	}

	close $data_in;
	close $data_out;
	
	# Report the number of reads remaining
	open my $log_fh, ">>$folder/blat_log.txt" or die "Unable to print to log file in filter 3: q20 > 20. $!\n";
	print $log_fh "Filtering read pairs by quality score\nPairs analyzed:\t$totalin\n\tFiltered low-quality pairs:\t$filtered\n\tHigh-quality long pairs:\t$totalout\n\n";
	close $log_fh;
	
	unlink "$folder/PEread_LTR_NonInternal_linker_long.txt";

	# Delete interim report.
	return;
}

##########################################################################################################################

# This subroutine prepares fasta files for BLAT alignment. The first 15 bases of the post-LTR
# sequence and the post-linker sequence are used to create a unique "tag" to identify 
# identical PCR products and the number of times each unique pair is seen counted. A forward
# and a reverse fasta file are created.

sub make_unique_pair_fasta {
	
	my ($LTR,$folder) = @_;
	open my $data_in, "<$folder/PEread_LTR_NonInternal_linker_long_q20.txt" or die "Unable to open the data file in filter 4: generate fasta $!\n";
	
	print "Identifying unique inserts...\n";
	
	# Initialize variables
	my $total_unique = my $all_reads = 0; # Count the total number of read pairs analyzed and how many of them are unique
	
	# Store the unique pairs and their data in a hash
	# %unique hash structure:
	# $unique{UNIQUE_PAIR_ID} = [TIMES PAIR IS SEEN, [@INPUT_LINE_COMPONENT_PARTS], $LTR_seq, $Linker_seq];
	my %unique; 
	
	# The LTR junction from the LTR primer file is required
	my $jct_file = "$LTR/$LTR"."_primer.txt"; 
	open my $jct_in, "<", $jct_file or die "$!: $jct_file";
	<$jct_in>;
	my $junction = <$jct_in>;
	close $jct_in;
	
	# Obtain the reverse complement of the LTR junction
	$junction = reverse $junction;
	$junction =~ tr/ATGC/TACG/;
	
	# Iterate through the data to search for unique pairs
	while (my $line = <$data_in>) {
		$all_reads++; # Count all pairs analyzed
		chomp($line);
		my @elements=split /\t/, $line; # Split the input line into it's component parts
		
		# Obtain the relevant parts and modify as necessary
		my $seq_ID = $elements[0]; # Obtain the sequence ID
		my $mid = $elements[10]; # Obtain the 10 base molecular identifier
		$seq_ID = $seq_ID."#".$mid; # Append the 10 base molecular identifier to the sequence ID
		my $LTR_seq = $elements[3]; # Obtain the post-LTR read 1 sequence
		my $Linker_seq = $elements[8]; # Obtain the post-linker read 2 sequence
		
		# Create an identifier ($tag) for each unique read by combining the first 15 bases each of $LTR_seq_end and $Linker_seq_end
		my $LTR_seq_end = substr($LTR_seq, 0, 15);
		my $Linker_seq_end = substr($Linker_seq, 0, 15);
		my $tag="$LTR_seq_end\t$Linker_seq_end";
		
		# If this is the first time seeing a unique pair, store it in the hash
		unless (exists $unique{$tag}) {
			$LTR_seq =~ s/(.*)AGTCCTCTAA.*/$1/; # Remove any linker sequence from read 1 -> AGTCCTCTAA
			$Linker_seq =~ s/(.*)$junction.*/$1/; # Remove any HIV LTR sequence from read 2
			$unique{$tag} = [1, [@elements], $LTR_seq, $Linker_seq]; # Set the uniqe pair count to 1. Store the relevant sequence data.
			$total_unique++; # Count the number of unique pairs
			
		# If the pair has been seen before, count it.	
		}else {$unique{$tag}[0]++;}
	}
	close $data_in;
	
	# Delete interim report.
	unlink "$folder/PEread_LTR_NonInternal_linker_long_q20.txt";
	
	# Print the unique pairs without their LTR and linker sequences
	open my $out_counts, ">$folder/PEread_LTR_NonInternal_linker_long_UniquePairWct.txt" or die "Unable to print the unique pair counts file in filter 4: $!\n";
	open my $outL, ">$folder/PEread_LTR_NonInternal_linker_long_UniquePairWct_F.fa" or die "Unable to print the forward fasta file in filter 4: $!\n";
	open my $outR, ">$folder/PEread_LTR_NonInternal_linker_long_UniquePairWct_R.fa" or die "Unable to print the reverse fasta file in filter 4: $!\n";	
	
	print "Generating fasta files...\n";
	
	foreach my $tag (keys %unique) {
		$unique{$tag}[1][0] .= "#$unique{$tag}[0]"; # Append the unique pair count to the sequence ID w/ MID
		my $newline = join("\t", @{$unique{$tag}[1]}); # Join the line elements with the modified sequence ID to make the output line
		print $out_counts "$newline\n"; # Print the new line to an interim report
		print $outL ">$unique{$tag}[1][0]\n$unique{$tag}[2]\n"; # Print the read 1 sequences to a fasta file
		print $outR ">$unique{$tag}[1][0]\n$unique{$tag}[3]\n"; # Print the read 2 sequences to a fasta file
	}
	
	close $out_counts;
	close $outL;
	close $outR;
	
	open my $log_fh, ">>$folder/blat_log.txt" or die "Unable to print to log file in filter 4: generate fasta. $!\n";
	print $log_fh "Count and merge unique read pairs for alignment\nPairs analyzed:\t$all_reads\n\tUnique, long, high-quality pairs: $total_unique\n\n";
	close $log_fh;
	
	return;
}

1;
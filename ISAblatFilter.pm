#!/usr/bin/perl

# ISAblatFilter.pm
# Dr. Xiaolin Wu & Daria W. Wells

# This module contains the subroutines that process the demultiplexed and pair-matched 
# ISA data and generates fasta files for alignment.

package ISAblatFilter; 

use strict;
use warnings;
use Exporter qw(import);

our @EXPORT_OK = qw(blat_filter1 blat_filter2 blat_filter3 blat_filter4 blat_filter5 blat_filter6 blat_filter7);


#Blat output format
#psLayout version 3
#
#match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts
#     	match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count
#---------------------------------------------------------------------------------------------------------------------------------------------------------------
#80	0	0	0	0	0	0	0	+	FTF2IS401C9OTF_LTR	123	42	122	chr1	197195432	43516678	43516758	1	80,	42,	43516678,
#71	0	0	0	0	0	0	0	-	FTF2IS401A43VI_LTR	77	0	71	chr1	197195432	82334892	82334963	1	71,	6,	82334892,


# This subroutine filters out alignments from the BLAT output that start more than 3 bases
# into the ISA sequence and are shorter than 20 bases.
sub blat_filter1 {
	# Create the input and output file names.
	my ($lib, $folder) = @_;
	my $libname = "$lib"."_blat.txt";
	my $file2open="$libname";
	my $file2out="$libname\_filter1";

	my @lines; # The passing alignments will be stored in @lines.

	open my $data_in, "<$file2open" or die "Unable $libname in BLAT filter 1: $!\n";

	my $num_passed_reads = my $num_rejected_reads = 0;
	# Test each line.
	while (my $line=<$data_in>) {
		next unless $line =~ /^\d/;  # Skip lines that don't start with a digit.
		my @eles=split(/\t/, $line);
		#start <4
		#Qgap[5] <10
		#Tgap[7] <10
		#identity >95% [0]/([12]-[11])
		if($eles[11]<4 and $eles[5]<10 and $eles[7]<10 and $eles[0]>19 and $eles[0]/($eles[12]-$eles[11])>0.95) {
			push(@lines, $line);
			$num_passed_reads++;
		}else{$num_rejected_reads++;}
	}
	close $data_in;

	open my $data_out, ">$file2out" or die "Unable to print data in BLAT filter 1: $!\n";

	# Sort @lines to organize the matches for the next filter step.
	my @sorted_lines = sort {
		my @a_fields = split /\t/, $a;
		my @b_fields = split /\t/, $b;

		$a_fields[9] cmp $b_fields[9]  # sort by seq ID
		||
		$b_fields[0] <=> $a_fields[0]  # then sort by match length
	} @lines;
	
	foreach my $line (@sorted_lines) {
		print $data_out "$line";
	}
	print $data_out "\n";
	close $data_out;
}


##########################################################################################

# This subroutine iterates through the sorted data from sub blat_filter1.  The alignments
# are in order by sequence ID. The subroutine tracks the highest score for each sequence
# ID and prints the reads that match only 1 time and the reads that have matches within 10
# bases of each other to different files.

sub blat_filter2 {
	# Create the input file name.
	my ($lib,$folder) = @_;
	my $libname = "$lib"."_blat.txt";
	my $file2open="$libname\_filter1";


	# Create the output file names.
	my $filename_outRepeat = $file2open;
	$filename_outRepeat =~s/filter1/filter2_repeatMatch/;
	my $filename_outBest = $file2open;
	$filename_outBest =~s/filter1/filter2_bestMatch/;
	
	
	open my $data_in, "$file2open" or die "Unable to read data in BLAT filter 2: $!\n";

	# Initialize the variables that will be used to track the scores from
	# the first line of the input file.
	chomp(my $line = <$data_in>);
	my @eles = split /\t/, $line;

	my ($sequence_id) = (split /\t/, $line)[9];
	
	# @seq_ID_entries will temporarily store the lines associated with each sequence ID
	# while checking for matche to repetitive elements. It's reset when the script moves on
	# to another sequence ID.
	my @seq_ID_entries = ($line); # First, make the first line of the file the first entry in @seq_ID_entries.
	my @out_best; # An array to store the results for _outBest
	my @out_repeat;# An array to store the results for _outRepeat
	# Iterate through the high quality blat results, searching for reads that mapped to multiple
	# places in the genome
	while (my $new_line = <$data_in>) {

		
		chomp($new_line);
		# Obtain the sequence ID of the next entry. Make $newsequence_ID a junk ID if we're at the end of the file.
		my ($newsequence_ID) = $new_line =~ m/\w/ ? (split /\t/, $new_line)[9] : "<<EOF>>";
		# If the next sequence ID is the same as the first one, append it to @seq_ID_entries and move on to the next comparison
 		if ($sequence_id eq $newsequence_ID) {
 			push @seq_ID_entries, $new_line;
 		}else{  # If the next sequence ID is different from the first, determine wich, if any, entries
 				# are considered multi-mapped.
 			my $tag_count = @seq_ID_entries;
 			if($tag_count == 1) { # If there is only one instance of a particular sequence ID, print it to _outBest
 				push @out_best, "$tag_count\t$line";
 			}else { # If there are more than 1 entry, check for multi-mapping

 				my $first_entry = shift @seq_ID_entries; # The top scoring entry is the first one to be compared
 				my @multi_mapped = ($first_entry); # Create an array for multi-mapped reads.
 				my ($first_entry_matches,$site) = (split /\t/, $first_entry)[0,15];
 				
 				while (my $next_entry = shift @seq_ID_entries) { # Check each potentially multi-mapped read in order of descending match score
 					
 					my ($next_entry_matches,$next_site) = (split /\t/, $next_entry)[0,15];
 					# Reads are considered multi-mapped if the match score difference of at least one other entry
 					# and the top entry core are within 10 matches of one another. 
 					# If a sequence is multi-mapped, add it to @multi_mapped then check if there are more remaining to be checked
 					if ($first_entry_matches - $next_entry_matches < 10) {
 						push @multi_mapped, $next_entry;
 						# If that was the last element to be checked, add the entries to @out_repeat. If not, move on to the next element
 						if (@seq_ID_entries == 0) {
 							my $multi_count = @multi_mapped;
 							push @out_repeat, "$multi_count\t$_" foreach @multi_mapped;
 						}
					# If the next entry isn't within 10 matches of the first, the multi-map check is complete
 					}else {
 						# Count the number of mapped reads
 						my $multi_count = @multi_mapped;
 						# If only one read maps, add it to @out_best and end the loop
 						if ($multi_count == 1) {
 							push @out_best, "$multi_count\t$first_entry";
 							@seq_ID_entries = undef;
 						# If more than one read map, add them to @out_repeat and end the loop
 						}else {
 							push @out_repeat, "$multi_count\t$_" foreach @multi_mapped;
 							@seq_ID_entries = undef;
 						}
 					}
 				}
			}
			# Prepare the variables for the next comparison
			$line = $new_line;
			@seq_ID_entries = ($line);
			$sequence_id = $newsequence_ID;
		}
	}
	# Print each array to the appropriate files. The log data will be written in the sub blat_filter3
	open my $out_repeat, ">$filename_outRepeat" or die "Unable to open repeat out in BLAT filter 2: $!\n";
	open my $out_best, ">$filename_outBest" or die "Unable to open best out in BLAT filter 2: $!\n";
	print $out_best "$_\n" foreach @out_best;
	print $out_repeat "$_\n" foreach @out_repeat;

	
	close $data_in;
	close $out_best;
	close $out_repeat;
	
	return;
}


##########################################################################################

# This subroutine reads the filtered BLAT data and converts it into the preferred 
# results format. The read1 and read2 paired-end read coordinates go into their respective 
# hashes with the sequence ID as the key and are printed at the end to two files. One for further
# filtering and another for Uclust. Read1 sequences without corresponding read2's are 
# included in the output for Uclust.

sub blat_filter3 {
	
	my $folder = shift;
	open my $in_L, "<$folder/PEread_LTR_NonInternal_linker_long_UniquePair_F_blat.txt_filter2_bestMatch" or die "$folder/PEread_LTR_NonInternal_linker_long_UniquePair_F_blat.txt_filter2_bestMatch : $!\n";
	open my $in_R, "<$folder/PEread_LTR_NonInternal_linker_long_UniquePair_R_blat.txt_filter2_bestMatch" or die "$folder/PEread_LTR_NonInternal_linker_long_UniquePair_R_blat.txt_filter2_bestMatch : $!\n";
	
	my %ID2corL; # Hash to store the left reads.
	my %ID2corR; # Hash to store the right reads.
	my %id2count; # Hash to count how many times a set of coordinates is seen.
	
	# @parameters consists of two array references. Each one stores the left/right input
	# file handle and a reference to the read hash. 
	my @parameters = ([$in_L, \%ID2corL],[$in_R, \%ID2corR]);
	
		foreach my $p (@parameters) { # Perform this analysis for each array reference in @parameters.
			my $fh = $p->[0]; # $p->[0] is the file handle to read from.
			while (my $line = <$fh>) { # Iterate through the input file
				next if $line !~ m/\w/;
				chomp($line);
				my @elements = split(/\t/, $line); # Make an array of each tab-separated line element.

				my $id = $elements[10]; # The sequence ID with the unique pair count.
				my ($idwocount, $MID, $count) = (split /\#/,$id); # Separate the ID from the count.
				$id = $idwocount.'#'.$MID;
				$id2count{$id} = $count; # Store the count.
			
				my $chrcor; # The sequence coordinates.
				if($elements[9] eq "+") { # The line elements for the coordinates depend on the strand.
					$chrcor="$elements[14]"."\t+\t"."$elements[16]"; # Example: chr17	+	3618824
				}else{
					$chrcor="$elements[14]"."\t-\t"."$elements[17]";
				}
			
				$p->[1]{$id} = $chrcor; # Store the coordinates. 
			}	
		}
	close $in_L;
	close $in_R;
	
	my @bestmatches; # Create a list of the read 1 best match sequence IDs and coordinates.
	# Iterate through the read 1 hash
	while (my ($id, $cor)=each(%ID2corL)) {
		my $bestmatch = "$id\t$cor"; # Create the first portion of the output line 
		push (@bestmatches, $bestmatch); # Add the portion to @bestmatches
	}
	
	# Sort @bestmatches.
	my @sorted_bestmatches = sort{
	# Split the string into two parts
	my @a_fields = split /\t/, $a;
	my @b_fields = split /\t/, $b;
		
	$a_fields[1] cmp $b_fields[1]  # First sort by chromosome.
	||
	$a_fields[3] <=> $b_fields[3]  # Then by ascending nucleotide coordinate.
	} @bestmatches;
	
	# Print the coordinates for every left-hand read along with the corresponding
	# right-hand read if there is one.
	# Example 9	M02560:63:000000000-ACWWU:1:1101:17574:2571#9	chr17	+	3618824	chr17	-	3618893
	open my $matched_out, ">$folder/PEread_LTR_NonInternal_linker_long_UniquePairWct_F_R_blat.txt" or die "Cannot open >$folder/PEread_LTR_NonInternal_linker_long_UniquePairWct_F_R_blat.txt for printing.\n.$!\n";
	open my $Uclust_out, ">$folder/PEread_LTR_NonInternal_linker_long_UniquePairWct_F_R_for_Uclust.txt" or die "Cannot open >$folder/PEread_LTR_NonInternal_linker_long_UniquePairWct_F_R_blat.txt for printing.\n.$!\n";
	
	# Print the sorted data
	foreach my $RIS (@sorted_bestmatches) {
		my @elements=split /\t/, $RIS; # Split the line into its component parts
		my $id = $elements[0]; # Obtain the sequence ID
		# If this sequence ID has two coordinates, print it to both files. If it has only one,
		# only print it to the Uclust file
		if (exists $ID2corR{$id}) {
			print $matched_out "$id2count{$id}\t$RIS\t$ID2corR{$id}\n";
			print $Uclust_out "$id2count{$id}\t$RIS\t$ID2corR{$id}\n";
		}else{
			print $Uclust_out "$id2count{$id}\t$RIS\n";
		}
	}

	close $matched_out;
	close $Uclust_out;
	
	# Open repeat match to print multi-map meta data to the log file
	my $repeat_L =  "$folder/PEread_LTR_NonInternal_linker_long_UniquePair_F_blat.txt_filter2_repeatMatch" or die "$folder/PEread_LTR_NonInternal_linker_long_UniquePair_F_blat.txt_filter2_repeatMatch : $!\n";
	my $repeat_R =  "$folder/PEread_LTR_NonInternal_linker_long_UniquePair_R_blat.txt_filter2_repeatMatch" or die "$folder/PEread_LTR_NonInternal_linker_long_UniquePair_R_blat.txt_filter2_repeatMatch : $!\n";
	my %seen_repeats;
	# Open each file handle separately. Count the number of unique sequence IDs seen
	foreach ($repeat_L, $repeat_R) {
		open my $handle, "<", $_;
		while (my $line = <$handle>) {
			(my $id) = (split /\t/, $line)[10];
			$seen_repeats{$id} = "" unless exists $seen_repeats{$id};
		}
		close $handle;
	}
	my $repeat_count = keys %seen_repeats;
	# Return the number of unique sequence IDs for use in sub blat_filter4
	return($repeat_count);
}

##########################################################################################

# This subroutine removes integration sites from the data if the paired-reads aren't on the
# same chromosome, are on the same strand, or are farther than 10,000 bases away from
# each other.

sub blat_filter4 {

	my ($folder, $repeat_count) = @_;
	open my $data_in, "<$folder/PEread_LTR_NonInternal_linker_long_UniquePairWct_F_R_blat.txt" or die "$folder/PEread_LTR_NonInternal_linker_long_UniquePairWct_F_R_blat.txt";
	open my $data_out, ">$folder/PairedEndRIS_filter1.txt" or die "Can't open $folder/PairedEndRIS_filter1.txt for printing.\n$!\n";
	my $total_unique_pairs = my $total_raw_pairs = my $all_pairs = 0;
	
	# Iterate through the data
	while (my $line=<$data_in>) {
		next unless $line =~ m/\w/; # Skip any blank lines, like the last one.
		chomp($line);
		$all_pairs++; # Count ever line analyzed
		my @eles = split /\t/, $line; # Split the line into its component parts
		if( # If the two coordinates meet the criteria, print and count them. Otherwise,
			# exclude them from the next report.
			($eles[2] eq $eles[5]) and 
			($eles[3] ne $eles[6]) and
			(abs($eles[4]-$eles[7]) <10000)
		) {
			print $data_out "$line\n";
			# Count the number of pairs that passed and the number of raw reads
			# associated with each passing pair
			$total_unique_pairs += 1;
			$total_raw_pairs += $eles[0];
		}
			
			
	}

	close $data_in,;
	close $data_out;
	
	open my $log_fh, ">>$folder/blat_log.txt" or die "Unable to print to log file in blat_filter5 $!\n";
	print $log_fh "Filter unmatched read1 and read2 coordinates\nPairs analyzed:\t$all_pairs\n\tUnique pairs mapped to 1 location:\t$total_unique_pairs\n";
	print $log_fh "\tRaw pairs mapped to 1 location:\t$total_raw_pairs\n\tRaw pairs with at least one read\n\t\tmapping to multiple locations:\t$repeat_count\n\n";
	close $log_fh;
	
	return;
}

##########################################################################################

# This subroutine determines which reads should be checked for potential merging.
# Integration sites within 10 bases are sent to sub merged_lines for further review.

sub blat_filter5 {

	my $folder = shift;

	my $file = "$folder/PairedEndRIS_filter1.txt";
	open my $fh, "<", $file;

	# Read all the lines. The lines are sorted by integration site.
	my @data_lines = <$fh>;
	my $mapped_unique_pairs = @data_lines;
	
	my $line = shift @data_lines; # Define the first line to be compared
	chomp($line);
	 
	# Split the line into individual scalars. The integration site coordinates will be compared.
	my ($junction_chr, $junction_base) = (split "\t", $line)[2,4];
	
	my @lines_to_merge; # Seqs that may potentially be merged will be stored here.
	my @merged_data; # The final merged lines will be stored here.
	
	push @lines_to_merge, $line; # Always add the first instance of an integration site to @lines_to_merge.
	
	# Iterate through the data, sending lines that may be merged to subroutine merge_lines.
	foreach my $next_line (@data_lines) {
		chomp($next_line);
	
		# First split the next line into individual scalars
		my ($next_junction_chr, $next_junction_base) = (split "\t", $next_line)[2,4];
		
		# Sequences with integration sites on the same chromosome and within 10 bases of each
		# other may potentially be merged. Store them in @lines_to_merge for now.
		if ($junction_chr eq $next_junction_chr and abs($junction_base - $next_junction_base) <= 10) {
			push @lines_to_merge, $next_line;
			if ($next_line eq $data_lines[-1]) {
				push @merged_data, &merge_lines(\@lines_to_merge);
			}	
	 	}else { # Move on to the next integration site if the integration sites aren't near each other.
			# If there is only one instance of an integration site, add it to the final data.
			# If not, send @lines_to_merge to &merged_lines.
			@lines_to_merge == 1 ? push @merged_data, @lines_to_merge : push @merged_data, &merge_lines(\@lines_to_merge);
			
			if ($next_line eq $data_lines[-1]) {
				push @merged_data, $next_line;
			}
				
			# Reset the variables for the next comparisons
			($junction_chr, $junction_base) = (split "\t", $next_line)[2,4];
			$line = $next_line;
			@lines_to_merge = ($line);
		}	
	}
	
	# Separate the MID from the sequence ID and append the fragment length to the end 
	# of each line of data. The fragment length is the breakpoint coordinate minus the 
	# integration site coordinate. Print the results to a tab-delimited text file.
	# Tally the final number of read pairs.
	my $final_pair_tally = my $merged_unique_pairs = my $HIV_internal_pairs = my $host_pairs = 0;

	open my $out, ">", "$folder/PairedEndRIS_filter2.txt";
	# Iterate through the post-merging data
	foreach my $line (@merged_data) {
		$merged_unique_pairs++; # Count the number of entries after merging
	 	my @elements = split "\t", $line; # Split each line into its component parts
	 	$final_pair_tally += $elements[0]; # Count the number of raw pairs after merging. None should be lost
	 	$elements[2] =~ m/HIV/ ? $HIV_internal_pairs += $elements[0] : $host_pairs += $elements[0];
	 	# Remove the MID from the sequence ID for the output. It's now one of the output line elements
	 	my ($sequence_id, $mid) = ($elements[1] =~ m/(^.*?)\#([ATGCN]{10})/);
	 	# Append the fuzz designation back to the sequence ID. It will be handled later
	 	$sequence_id .= "#$1" if $elements[1] =~ m/(fuzz\?{0,1})/;
	 	my $line_out = "$elements[0]\t$sequence_id\t$mid\t".join("\t", @elements[2..7])."\t".abs($elements[7] - $elements[4])."\n";
		print $out $line_out;
	}
	close $out;
	
	open my $log_fh, ">>$folder/blat_log.txt" or die "Unable to print to log file in blat_filter6 merge sites. $!\n";
	print $log_fh "Merge mapped unique pairs\nUnique pairs before merging:\t$mapped_unique_pairs\n\tUnique pairs after merging:\t$merged_unique_pairs\n\tRaw pairs mapped:\t$final_pair_tally\n";
	print $log_fh "\t\tInternal HIV pairs:\t$HIV_internal_pairs\n\t\tHost integration site pairs:\t$host_pairs\n\n";
	close $log_fh;
}
	
##########################################################################################
# This subroutine merges integration sites determined to have come from the same read.	
# Sites are merged based if they meet all of the following criteria:
#	1.	The breakpoints are within 10 bases of each other.
#	2.	The MIDs differ by no more than two bases. (Hamming distance)
#	3.	The site with the larger count is greater than double the smaller count plus 1. (count score)
# Sites are then labelled as fuzz or potential fuzz if:
#	1.	If the breakpoints are within 5 base pairs of one another AND the smaller read count is 5% or 
#		less of the larger count, label it "fuzz"
#	2.	If the breakpoints are within 5 base pairs of one another AND the smaller read count is between
#		5% and 10% of the larger read count, label it "fuzz?"
#
	
sub merge_lines {
	my $lines_to_merge = shift;

	my %all_lines; # Unsorted lines wil go here
	my @sorted_lines; # Will hold the lines sorted by descending read count
	my @merged_lines; # Will hold the the merged lines from each pass
	my %already_merged; # Will indicate which lines have been merged on a previous pass.
	my @returned_lines; # The final array that will be returned
	
	# Each line will be stored as an array in %all_lines with $current_seq_id as the key
	foreach my $line (@{$lines_to_merge}) {
		my @split_line = split "\t", $line;
		$all_lines{$split_line[1]} = [@split_line];
	}

	# Sort the lines by descending read count.
	@sorted_lines = sort{
		$all_lines{$b}->[0] <=> $all_lines{$a}->[0]
	}keys %all_lines;


	# Merge the integration sites if all the following criteria are 1:
			#	1.	The breakpoints are within 10 bases of each other.
			#	2.	The MIDs differ by no more than two bases. (Hamming distance)
			#	3.	The site with the larger count is greater than double the smaller count plus 1. (count score)
			
	my $last_seq_ID = $sorted_lines[-1]; # Used to indicate the final array element.
		
	# Begin searching for matches.
	while (@sorted_lines != 0) {	
		
		# Compare the first (most reads) array element to the rest of the array.
		my $current_seq_id = shift @sorted_lines; # Remove the first element
		
		# Skip to the next loop if $current_seq_id has already been merged.
		next if exists $already_merged{$current_seq_id};
		
		# If $current_seq_id is the last element of the array AND it has already been
		# merged with a previous line, end the loop. If not, add it to the results.
		if ($current_seq_id eq $last_seq_ID) {
			exists $already_merged{$current_seq_id} ? last : push @merged_lines, $all_lines{$current_seq_id};
		}
	
		# Compare $current_seq_id to the rest of the sorted lines.
	 	foreach my $next_seq_id (@sorted_lines) {

			# First, check if the next array element has already been merged.
		 	if (exists $already_merged{$next_seq_id}) {
		 	
		 		# If the next element is the last element of the array and $current_seq_id
		 		# has already been merged in a previous loop, add current line to the results
		 		# array. If this is not the last element of the array, move on to the next loop.
		 		$next_seq_id eq $last_seq_ID ? (push @merged_lines, $all_lines{$current_seq_id} and next) : next;
		 	}
		 	
		 	# Begin the comparison.
	 		(my $mid) = ($all_lines{$current_seq_id}[1] =~ m/#([ATGCN]{10})$/); # Obtain the MID from the current sequence ID
	 		(my $next_mid) = ($all_lines{$next_seq_id}[1] =~ m/#([ATGCN]{10})$/); # Obtain the MID from the next sequence ID

	 		my $break_diff = abs($all_lines{$current_seq_id}[7] - $all_lines{$next_seq_id}[7]); # Calculate the distance between the two breakpoints.
	 		my $hamming = ($mid ^ $next_mid) =~ tr/\0//c; # Calculate the Hamming distance between the two MIDs.
	 		my $score = 2 * $all_lines{$next_seq_id}[0] + 1; # Calculate the count score of the comparison.
		
			# Merge and store the data_lines if all three criteria are met
			if ($break_diff == 0 or ($break_diff <= 10 and $hamming <= 2 and $all_lines{$current_seq_id}[0] > $score)) {
				
				$all_lines{$current_seq_id}[0] += $all_lines{$next_seq_id}[0]; # Add the smaller count to the larger count.
				$already_merged{$next_seq_id} = ""; # Indicate that the smaller count line has been merged.
				
				# This was the last element, add the current line to the results array and
				# indicate that it has been merged.
				if ($next_seq_id eq $last_seq_ID) {
					push @merged_lines, $all_lines{$current_seq_id};
					$already_merged{$current_seq_id} = "";
				}
			}elsif($next_seq_id eq $last_seq_ID) {
				push @merged_lines, $all_lines{$current_seq_id} unless exists $already_merged{$current_seq_id};
				$already_merged{$current_seq_id} = "";
			}
		}
	}
	# After all lines are merged that can be, check to see if any of the remaining reads are
	# considered fuzz.
	# Sites are then labelled as fuzz or potential fuzz if:
	#	1.	If the breakpoints are within 5 base pairs of one another AND the smaller read count is 5% or 
	#		less of the larger count, label it "fuzz"
	#	2.	If the breakpoints are within 5 base pairs of one another AND the smaller read count is between
	#		5% and 10% of the larger read count, label it "fuzz?"
	while (my $line = shift @merged_lines) {
		push @returned_lines, join("\t", @{$line});
		my @next_loop;
			while (my $next_line = shift @merged_lines) {
				if (abs($line->[7] - $next_line->[7]) < 5) {
				my $diff = $next_line->[0] / $line->[0];
					if ($next_line->[0] / $line->[0] <= 0.05) {
						$next_line->[1] .= "#fuzz";
						push @returned_lines, join("\t", @{$next_line});
					}elsif ($next_line->[0] / $line->[0] <= 0.1) {
						$next_line->[1] .= "#fuzz?";
						push @returned_lines, join("\t", @{$next_line});
					}else {
						push @returned_lines, join("\t", @{$next_line});
					}
				} else {
					push @next_loop, $next_line;
				}
			}
		@merged_lines = @next_loop;
	}	
	
 	
 	return @returned_lines;
}

1;
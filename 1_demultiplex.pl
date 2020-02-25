#!/usr/bin/perl

# 1_demultiplex.pl
# Dr. Xiaolin Wu & Daria W. Wells

# This script is used to filter Illumina paired-end ISA sequence data from the MiSeq, 
# HiSeq, or NextSeq. It will search for the user-specified barcodes, remove unspecified barcodes, and 
# and split the reads at the LTR/Linker junctions for downstream processing. It 
# requires knowledge of the bases immediately upstream to the junction. The sequence of the
# i7 barcode used for demultiplexing is not used, only the name. It will use the name 
# found in the settings file, so be sure to use the correct one.

use strict;

use IO::Zlib qw(:gzip_external 1);
use Carp;

use lib '/Users/wellsd2/Desktop/perl_scripts/ISA/current_pipeline'; # The path of the Illumina barcodes directory.
use IlluminaBarcodes; # The module that returns a hash of Illumina barcodes. Update as needed.

# Require settings file on the command line.
croak "No settings file specified on the command line." unless @ARGV;


# Import your settings from the settings file. Each value is a list of the settings under
# header, allowing you to look for multiple barcodes and analyze both LTRs or only one.
my (%settings, $header); # Initialize the settings hash and the header variable.
open my $settings_fh, "<$ARGV[0]" or die "Cannot read the settings file $ARGV[0]: $!"; # Open the settings file for reading.
while (my $settings_line = <$settings_fh>) { # 
	chomp($settings_line);
	if ($settings_line =~ /<(.*)>/) { # Text surrounded by < and > indicates a settings type header, e.g. <LTR>.
		$header = $1; # Assign the matched text to $header.
		next;
	}
	push @{$settings{$header}}, $settings_line; # If the line is not a header, add it to the list for the most recent header.
}
close $settings_fh; # Finished with settings file.

# Create hashes for to hold the barcodes and file names.
my %all_barcodes = IlluminaBarcodes::get_barcodes(); # Obtain the hash of all Illumina barcodes.
my %barcodes2filename; # Store the filenames in a hash with the barcode sequences as keys.

# Create file names to be used for data output
foreach my $LTR (@{$settings{'LTR'}}) {
	
	# Use the $LTR variable to create a file containing the primer and junction 
	# sequences for later use.
	mkdir $LTR;
	my $primer_file = ">$LTR/$LTR"."_primer.txt";
	my $primer_key = "$LTR"."_primer";
	open my $primer_fh, "$primer_file";
	print $primer_fh @{$settings{$primer_key}};
	my $junction_key = $LTR."_junction";
	print $primer_fh "\n@{$settings{$junction_key}}";
	
	# Cycle through the list of barcodes.
	foreach my $PE1_barcode (@{$settings{PE1barcodes}}) {
		# Make a file name for each barcode listed in the settings file.
		my $filename = "PE1_"."$PE1_barcode"."_PE2_"."@{$settings{PE2barcodes}}[0]".".txt";
		# Croak if the file already exists.
		my $test = "$LTR/$filename";
		croak "\nError: $LTR/$filename already exists\n" if -e $test;
		# Store the file names in a hash with the barcode sequence string as the key.
		# This hash will be used to identify sequences later
		# The file names wil be used to write data to files later			
		$barcodes2filename{$LTR}{$all_barcodes{$PE1_barcode}}=$filename;
	}
}


# Use subroutine choose_LTR_sub to return a subroutine reference to analyze either the 
# 3LTR only, 5LTR only, or both. Pass to choose_LTR_sub each LTR primer and junction 
# sequence, as well as the list of LTR's to be analyzed.
my $match_LTR_seq = choose_LTR_sub(@{$settings{"3LTR_primer"}},@{$settings{"3LTR_junction"}},@{$settings{"5LTR_primer"}},@{$settings{"5LTR_junction"}},@{$settings{"LTR"}});

# Create a list of the numbers from 1 to n in the format 001, ... 00n with 
# n being the number of Read 1 & Read 2 pairs for analysis.
my @files2read; # Initialize the array
for (my $i = 1; $i <= $settings{pairs}[0]; $i++) { # Perform this loop for all integers from 1 to n.
	my $pairs = sprintf("%03d", $i); # Format the number with leading zeros.
	push @files2read, $pairs; # Add the formatted number to @files2read.
}

# Create the file name template from the file name in the settings.
# The template is the text before R1 or R2, e.g. "DNAsample_S1_L001_"
my ($filenametemplate) = ($settings{file_name_template}[0] =~ /(.*)[R]\d_\d{3}\.fastq\.gz/) or croak "Incorrect fastq.gz file name format. Check your settings file.\n";

# Open a file handle for printing the log file.
my $log_out = ">>demultiplex_log.txt";
open my $log_fh, $log_out or die "Unable to open print to the log:$!\n";

# Print the file name template and the start time to the log file and/or terminal.
my $start_time = localtime;
print "\n$filenametemplate\n";
print $log_fh "Demultiplex script used: $0\n"; 
print $log_fh "Analysis started at $start_time\n\n";
print $log_fh "$filenametemplate\n";

# Add './' to the beginning of the file name template.
$filenametemplate="./$filenametemplate";

# Initialize the counter for the total number of paired end reads across all suquence files analyzed.
my $grand_total_reads = my $all_3LTR = my $all_5LTR = my $all_unidentified = 0;

# Send the required variables and references to the extract_and_print_data subroutine.
# This subroutine only returns the nubmer of reads per pass.
foreach my $pairs (@files2read) { # Perform this loop for each three digit Read1/Read2 pair number.
	(my $total_reads, my $_3LTR_count, my $_5LTR_count, my $unidentified) = extract_and_print_data(
		$pairs, $filenametemplate, \%barcodes2filename, $settings{'LTR'}, $match_LTR_seq, $log_out);
	# Count how many read pairs are present, how many were found for each LTR, and how many were rejected.	
	$grand_total_reads += $total_reads;
	$all_3LTR += $_3LTR_count;
	$all_5LTR += $_5LTR_count;
	$all_unidentified += $unidentified;	
}

my %LTRs; # This hash is only used to determine which LTRs should be reported in the log
$LTRs{$_} = undef foreach @{$settings{'LTR'}};
# Print the number of paired end reads to the terminal and log file.
open $log_fh, $log_out or die "Unable to print to the log:$!\n";

print "*" x 40 . "\n\n";
printf "Total 3LTR read pairs:\t%s\n", exists $LTRs{'3LTR'} ? $all_3LTR : "N/A";
printf "Total 5LTR read pairs:\t%s\n", exists $LTRs{'5LTR'} ? $all_5LTR : "N/A";
print "Total unidentified:\t$all_unidentified\n\n";
print "Grand total:\t$grand_total_reads\n\n\n";
	
print $log_fh "*" x 40 . "\n\n";
printf $log_fh "Total 3LTR read pairs:\t%s\n", exists $LTRs{'3LTR'} ? $all_3LTR : "N/A";
printf $log_fh "Total 5LTR read pairs:\t%s\n", exists $LTRs{'5LTR'} ? $all_5LTR : "N/A";
print $log_fh "Total unidentified pairs:\t$all_unidentified\n\n";
print $log_fh "Grand total:\t$grand_total_reads\n\n";
print "Compressing...\n\n\n";

# This script will create a file for every barcode combination, even if no data exists for it.
# This loop deletes the empty files.
foreach my $LTR (@{$settings{LTR}}) {
	my @files = <$LTR/PE1*PE2*.txt>;
	foreach (@files) {
		unlink unless -s;
		system("gzip $_") if -e;
	}
}

# Print the time that the analysis completed to the log.
my $end_time = localtime;
print $log_fh "Analysis completed at $end_time\n\n";
print $log_fh "#" x 90 . "\n\n\n";
close $log_fh;
exit;

# This subroutine extracts the data from the fastq.gz files, identifies the LTR region 
# based on the sequences in the settings file, and prints the data to a .txt file in a
# new directory named after the 3LTR or 5LTR. The data disappears after each 
# pass through the subroutine.
sub extract_and_print_data {
	# Initialize the passed-in variables.
	my($pairs, $filenametemplate, $barcodes2filename, $analyses, $match_LTR_seq, $log_out) = @_;
	
	# Open the log file handle for printing
	open my $log_fh, $log_out or die "Unable to print to the log:$!\n";
	
	# Extract the contents of the fastq.gz files to temporary locations for processing.

	# Create the names of the files to extract using the file name template and the current $pairs.
	my $R1="$filenametemplate"."R1_"."$pairs".".fastq.gz";
	die "\nCan't find file $R1\n" unless -e $R1;
	my $R2="$filenametemplate"."R2_"."$pairs".".fastq.gz";
	die "\nCan't find file $R2\n" unless -e $R2;
	
	# Print the file names to the monitor and the log.
	print "$R1\n$R2\n\n";
	print $log_fh "$R1\n$R2\n\n";

	# Open file handles for printing the sequence data.
	open(my $INread1, "gzip -dc $R1 |") or die "Unable to extract the read 1 sequences: $!";
	open(my $INread2, "gzip -dc $R2 |") or die "Unable to extract the read 2 sequences: $!";
	
	# Initialize to 0 the read counters for this pass.
	my $totalreads = my $passedreads = my $_3LTR_all_reads = my $_5LTR_all_reads = my $unidentified = 0;
	
	# Initialize a hash to store the data from this pass through the subroutine.
	# The data is printed and wiped from memory at the end.
	my $data_storage = build_data_hash($barcodes2filename, $analyses);
	
	
	# Begin unzipping and debarcoding the data
	SCAN: while (my $lineread1=<$INread1>) { # Begin with read 1
		chomp($lineread1);
	
		if($lineread1=~/^@(.*)( 1.*)/) { # Lines beginning with @ indicate sequence name 
			$totalreads++; 
			my $read1_id=$1; # Only want the part of the ID that's identical to read 2
			my $read1_seq=<$INread1>; # The next line is the sequence
			chomp($read1_seq);
			<$INread1>;  # Skip the + sign in between seq and quality
			my $read1_qual=<$INread1>; # Pull in the quality line
			chomp($read1_qual);
		
			# Repeat the same steps for read 2
			my $read2_id=<$INread2>;
			chomp($read2_id);
			$read2_id =~s/\@(.*) .*/$1/;
			my $read2_seq=<$INread2>;
			chomp($read2_seq);
			<$INread2>;
			my $read2_qual=<$INread2>;
			chomp($read2_qual);
			
			# Check to see if the sequence IDs match. Exit if they don't.
			if($read1_id ne $read2_id){
				print "Read 1 sequence ID does not match read 2:\n$read1_id\t$read2_id\n";
				print $log_fh "Read 1 sequence ID does not match read 2:\n$read1_id\t$read2_id\n";
				print "The read 1 and 2 sequence IDs don't match. There might be something wrong with your data.\n";
				print "All sequence ID must match to complete the analysis.";
				exit;
			}
			
			# Match the sequence to 3LTR, 5LTR, or both, depending on user input.
			my ($R1_junction, $LTR) = &$match_LTR_seq($read1_seq);
			
			# If neither nested primer was detected, skip to the next sequence.
			unless (defined $R1_junction and defined $LTR) {
				$unidentified++;
				next SCAN;
			}
			
			# Obtain the barcode sequence from read 1 to recall the file name
			my $read1_bc=substr($read1_seq,6,8); # Extract the read 1 barcode, skipping the dogtag.
			$unidentified++ and next SCAN unless exists $barcodes2filename->{$LTR}{$read1_bc};
			my $filename=$barcodes2filename->{$LTR}{$read1_bc}; # Recall file name from hash
			
			# Only compare the LTR and barcode sequences to the first 75 bases of the reads.
			my $read1_first75=substr $read1_seq,0,75; 
			my $read2_first75= substr $read2_seq,0,75;
			
			# Store the data for later printing if the junction and linker sequences.
			if	($read1_first75 =~ m/$R1_junction/ and $read2_first75 =~ m/TCCGCTTAGAGGACT/) {
				# Keep track of the number of 3LTR and 5LTR reads.
				$LTR eq "3LTR" ? $_3LTR_all_reads++ : $_5LTR_all_reads++;
				# Split read 1 based on the HIV/host junction sequence.
				my ($preLTRend, $postLTRend) = split (/$R1_junction/, $read1_seq);																						
				my $LTR_junction= "$preLTRend$R1_junction\t$postLTRend";
				# Insert a tab at the end of the linker sequence.
				my ($preLinker, $postLinker) = split (/TCCGCTTAGAGGACT/, $read2_seq);
				$read2_seq = "$preLinker"."TCCGCTTAGAGGACT\t$postLinker";
				# Get the molecular identifier (MID), the first 10 bases of read 2.
				my $MID = substr $read2_seq, 0, 10;
				# Add the passed sequences to the data storage hash.
				# Append a # then the MID to the end of each sequence ID.
				push @{$data_storage->{$LTR}{$filename}}, "$read1_id#$MID\t1\t$LTR_junction\t$read1_qual\t$read2_id#$MID\t2\t$read2_seq\t$read2_qual\t$MID";
			}else {$unidentified++;}
		}
	}		
	close $INread1;
	close $INread2; 
	
	# Print out the sequence data from this pass through the subroutine to files and print the 
	# read counts to the terminal and log file. 
	foreach my $LTR (sort keys %{$barcodes2filename}) { # Perform for each LTR's data.
		print "$LTR\n";
		print $log_fh "$LTR\n";
		
		# Use the %barcodes2filename keys to look up data from %data_storage
		foreach my $barcodes (sort keys %{$barcodes2filename{$LTR}}) { 
			my $filename = $barcodes2filename->{$LTR}{$barcodes};
			if (exists $data_storage->{$LTR}{$filename}) {
				my $counts = @{$data_storage->{$LTR}{$filename}};
				
				# Print the LTR counts to the terminal and the log file.
				$LTR eq "3LTR" ? { 
					print "\t$filename\n\t\tIdentified 3LTR read pairs:\t$_3LTR_all_reads\n\n"
					and print $log_fh "\t$filename\n\t\tIdentified 3LTR read pairs:\t$_3LTR_all_reads\n\n"
				}:	
				{
					print "\t$filename\n\t\tIdentified 5LTR read pairs:\t$_5LTR_all_reads\n\n"
					and print $log_fh "\t$filename\n\t\tIdentified 5LTR read pairs:\t$_5LTR_all_reads\n\n"
				};
			}
			
			# Print the data to a file.
			my $file2print="$LTR/$filename";
			open my $data_out, ">>$file2print" or die "Can't open $file2print for printing: $!\n";			
			foreach my $line (@{$data_storage->{$LTR}{$filename}}) {
				print $data_out "$line\n";
			}
			close $data_out;
		}
	}
	
	# Print the unidentified read counts to the terminal and the log file.
	print "Unidentified pairs:\t$unidentified\n\n";
	print $log_fh "Unidentified pairs:\t$unidentified\n\n";
	return ($totalreads, $_3LTR_all_reads, $_5LTR_all_reads, $unidentified);	
}

###################################################################################################


# This subroutine initializes the hash to store the data each pass through the main data subroutine. 
# The subroutine receives references for %barcodes2filename and $settings{LTR}.
# The sequences for output will be stored in a hash divided by LTR, then by filename,
# then as a list of sequences for each file, e.g. $data_storage{3LTR}{file_name.txt}[$seq_1, $seq_2, $seq_3].
# %data_storage disappears at the end of each pass through sub extract_and_print_data before
# parsing the next pair of reads, if applicable.
sub build_data_hash {
	my ($barcodes2filename, $analyses) = @_; # Obtain the input references.
	my %data_storage; # Initialize the data storage hash.
	foreach my $LTR (@{$analyses}) { # Set up the hash for one or both LTR's.
		foreach my $barcodes (keys %{$barcodes2filename}) { # Use %barcodes2filename to assign file names to the data storage hash.
			my $filename = $barcodes2filename->{$barcodes}; # Assign the current file name to $filename.
			$data_storage{$LTR}{$filename}=[]; # Initialize the anonymous array for storing sequence.
		}
	}
	return \%data_storage; # Return the a reference for the data storage hash.
}



###################################################################################################



# This subroutine returns a subroutine reference that will be used to analyze the LTR data
sub choose_LTR_sub {
	my %IUPAC = (  # This hash replaces the ambiguous bases in the primer with a regex
	"R"	=> "[AG]", # character class for downstream matching and substitutions.
	"Y"	=> "[CT]",
	"S"	=> "[GC]",
	"W"	=> "[AT]",
	"K"	=> "[GT]",
	"M"	=> "[AC]",
	"B"	=> "[CGT]",
	"D"	=> "[AGT]",
	"H"	=> "[ACT]",
	"V"	=> "[ACG]",
	"N" => "[ATGC]",
	);
	
	my($_3LTR_primer, $_3LTR_junction, $_5LTR_primer, $_5LTR_junction, @analyses)=@_; # Initialize the passed in variables.
	$_ = uc $_ foreach ($_3LTR_primer, $_3LTR_junction, $_5LTR_primer, $_5LTR_junction); # Make all base symbols uppercase if they're not already
	my $_3LTR_primer_length = length($_3LTR_primer); # Store the 3LTR primer length as a scalar variable.
	my $_5LTR_primer_length = length($_5LTR_primer); # Store the 5LTR primer length as a scalar variable.
	# Replace the ambiguous bases in the primers and junctions with character classes.
	foreach my $seq ($_3LTR_primer, $_3LTR_junction, $_5LTR_primer, $_5LTR_junction) {
		while ($seq =~ /([RYSWKMBDHVN])/g) { # Global match searches every base in the sequence.
			$seq =~ s/($1)/$IUPAC{$1}/; # Substitute the ambiguous base with a character class of its possbile bases.
		}
	}
	
	# Decide which analysis subroutine to use.  
	my $num_choices = @analyses; # Store the number of items in @analyses as a scalar variable.
	foreach (@analyses) {croak "Invalid LTR option selected. Check your settings file.\n" unless m/^[35]LTR$/;} # Croak if anything other than 3LTR or 5LTR is in %analyses.
	croak "More than 2 options were entered for LTR analysis. Check your settings file.\n" if $num_choices > 2;  # Croak if more than 2 items were entered for LTR analyses.
	croak "You entered $analyses[0] twice. Check your settings file." if $num_choices == 2 and $analyses[0] eq $analyses[1]; # Croak if the same LTR was entered twice.
	if ($num_choices == 2) { # If two valid LTR options were selected, return the subroutine defined from lines 318 to 332.
		return sub {
			my ($read1_seq) = @_; # Assign the read 1 sequence to a scalar. It will be the only value passed to this subroutine. 
			my $read1_nestedPrimer; # Initialize the primer variable.
			# Starting at base 15 of read1_seq, obtain a substring equal to the length of the 3LTR primer and assign it to $read1_nestedPrimer.
			# if the $read1_nestedPrimer is an exact match, return the LTR junction sequence and a string indicating that this sequence is from the 3LTR.
			if ($read1_nestedPrimer=substr($read1_seq,14,$_3LTR_primer_length) and $read1_nestedPrimer =~ /$_3LTR_primer/) {
				my ($R1_junction) = "$_3LTR_junction";
				return($R1_junction, "3LTR");
				
			# If the sequence doesn't match the 3LTR primer, repeat the steps to check if it matches the 5LTR.
			}elsif ($read1_nestedPrimer=substr($read1_seq,14,$_5LTR_primer_length) and $read1_nestedPrimer =~ /$_5LTR_primer/) {
				my ($R1_junction) = "$_5LTR_junction";
				return($R1_junction, "5LTR");
			
			# If the sequence matches neither LTR primer, return 0, instructing the script to skip this sequence.	
			}else {
				return 0;
			}
		};
	
	# Use the same steps above but only check for a match to the 3LTR primer.
	}elsif ($analyses[0] eq "3LTR") {
		return sub {
			my ($read1_seq) = @_;
			my $read1_nestedPrimer=substr($read1_seq,14,$_3LTR_primer_length); # Extract nested primer sequence
			return 0 unless $read1_nestedPrimer =~ /$_3LTR_primer/;
			my ($R1_junction) = "$_3LTR_junction";
			return($R1_junction, "3LTR");
		};
		
	# Use the same steps above but only check for a match to the 5LTR primer.	
	}elsif ($analyses[0] eq "5LTR") {
		return sub {
			my ($read1_seq) = @_;
			my $read1_nestedPrimer=substr($read1_seq,14,$_5LTR_primer_length); # Extract nested primer sequence
			return 0 unless $read1_nestedPrimer =~ /$_5LTR_primer/;
			my ($R1_junction) = "$_5LTR_junction";
			return($R1_junction, "5LTR");
		};
	}
}
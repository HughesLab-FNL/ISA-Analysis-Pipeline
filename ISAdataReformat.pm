#!/usr/bin/perl

# ISAdataReformat.pm
# Dr. Xiaolin Wu & Daria W. Wells

# This module contains the subroutines required to parse the integration site data 
# into different formats.

package ISAdataReformat;

use strict;
use warnings;
use Exporter qw(import);
use IO::Zlib qw(:gzip_external 1);



our @EXPORT_OK = qw(annotate_blat_file sort_fasta_by_length resort_Uclust extract_sequences);

##########################################################################################################################

# This subroutine adds information to the blat file. It adds the 30 bases upstrean and
# downstream of the LTR junction and checks if in the 30 upstream bases there is a stretch 
# of 10 bases where 7 or more bases match the last 10 bases of the LTR primer. It also
# adds the raw read 1 and read 2 sequences to the lines.

sub annotate_blat_file {
	my ($blat_file,$LTR,$name,$path2hg19,$folder) = @_;
	# Remove the "chr" from the end of the hg19 path for the getseq_fast subroutine.
	$path2hg19 =~ s/chr$//;
	# Open the primer sequence file and extract the primer information
	my $primer_file = "$LTR/$LTR"."_primer.txt";
	open my $primer_in, "$primer_file" or die "$primer_file: $!";
	my $primer = <$primer_in>;
	my $LTR_junction = <$primer_in>;
	close $primer_in;
	my $primer_end = substr($primer, -10, 10);
	my @primer_end_bases = split_primer($primer_end);
	$blat_file = "$folder/$blat_file";
	# Open file to read
	open my $blat_in,"<", $blat_file or die "$!: $blat_file\n";
	
	# Import integration site coordinates.  Pass the coordinates to subroutine getseq_fast
	# to extract the sequence and indicate whether it will likely cause mis-priming.
	# Store the data for printing after the loop;
	my @all_lines; # The output lines are stored in @output.
	my %indicator_results; # %indicator_results stores the result of the mis-prime test.
	# Test each line individually.
	while (<$blat_in>) {
		chomp;
		my @elements = split /\t/, $_;
		# Oobtain the genomic coordinate
		my $coordinate = "$elements[3]$elements[4]$elements[5]";
		# If this coordinate hasn't been seen already, skip it. If it hasn't, determine if the coordinate
		# is spike-in control, possibly the result of mispriming, and/or if the genomic DNA surrounding the
		# integration site is contains an exact match to the LTR junction.
		unless (exists $indicator_results{$coordinate}) {

			my ($indicator, $sequence) = getseq_fast($elements[1],$elements[3], $elements[4], $elements[5],$path2hg19,$LTR, $LTR_junction, \@primer_end_bases);
			my $first30bases = substr $sequence,0,30,"\t";
			$sequence = $first30bases.$sequence;
			if ($first30bases =~ m/$LTR_junction/) {
				if ($indicator =~ m/no_indicator/) {
					$indicator =~ s/no_indicator/Exact match to LTR junction/;
				}else {
					$indicator =~ s/\t(.*)/\tExact match to LTR junction \& $1/;
				}
			}
			$indicator_results{$coordinate} = "$indicator\t$sequence";
		}
		my $line_out = join("\t",@elements, $indicator_results{$coordinate});
		push @all_lines, $line_out;	
	}
	close $blat_in;

	$_ =~ s/#\d+$// foreach @all_lines;
	# Add the raw sequences to the file.
	my $dbc_file = "$LTR/$name.txt.gz";
	my %data_lines = extract_sequences(\@all_lines,$dbc_file);
	my @sorted_lines = sort {
		$data_lines{$a}[3] cmp $data_lines{$b}[3]
		||
		$data_lines{$a}[5] <=> $data_lines{$b}[5]
		||
		$data_lines{$a}[9] <=> $data_lines{$b}[9]
	} keys %data_lines;
	
	# Print the results
	my $out_file = my $out_fuzz = "$blat_file";
	$out_file =~ s/\.txt/_annotated\.txt/;
	$out_fuzz =~ s/\.txt/_annotated_with_fuzz\.txt/;
	open my $blat_out,">", "$out_file" or die "$out_file: $!";
	open my $fuzz_out,">", "$out_fuzz" or die "$out_fuzz: $!";
	print $blat_out "Count\tSeq ID\tMID\tJunction chr\tJunction strand\tJunction base\tBreak point chr\tBreak point strand\tBreak point base\tLength\tHIV LTR\tProvirus orientation to chromosome\tPrediction\tUpstream 30\tDownstream 30\tRead 1\tRead 2\n";
	print $fuzz_out "Count\tSeq ID\tMID\tJunction chr\tJunction strand\tJunction base\tBreak point chr\tBreak point strand\tBreak point base\tLength\tHIV LTR\tProvirus orientation to chromosome\tPrediction\tUpstream 30\tDownstream 30\tRead 1\tRead 2\n";

	foreach (@sorted_lines) {
		if (defined $data_lines{$_}[12]) {
 			$data_lines{$_}[12] =~ s/no_indicator//;
 			$data_lines{$_}[12] = $1 if $data_lines{$_}[1] =~ s/\#(fuzz\?{0,1})//; 
 		}
		my $line = join ("\t", @{$data_lines{$_}})."\n";
		print $blat_out $line unless $line =~ m/fuzz/;
		print $fuzz_out $line;
	}
	close $blat_out;
	close $fuzz_out;
	return;
}

##########################################################################################################################

# This subroutine is used for the mis-priming check. It extracts and returns a DNA sequence 
# from an hg19 chromosome fasta file based on the coordinates passed to it.
	
sub getseq_fast {
	
	# Initialize the sequence and coordinate variables.
	my $primer_bases = pop;
	my ($sequence_ID, $chr, $strand, $IS_base, $path2hg19, $LTR) = @_;
	# Get the virus orientation relative to the chromosome.
	my $virus_strand;
	if ($LTR eq "3LTR" and $strand eq "+") {
		$virus_strand = "$LTR(+)\t+";
	}elsif($LTR eq "3LTR" and $strand eq "-") {
		$virus_strand = "$LTR(+)\t-";
	}elsif($LTR eq "5LTR" and $strand eq "+") {
		$virus_strand = "$LTR(-)\t-";
	}elsif($LTR eq "5LTR" and $strand eq "-") {
		$virus_strand = "$LTR(-)\t+";
	}
	
	# Modify the name of chromosome unknown to match the fasta file.
	$chr = "chrUn" if $chr =~ /chrUn/;
	# Don't check for matches to the HIV genome. Just return an empty result.
	if ($chr =~ /HIV-1|HIV1/) {
		my $indicator = my $seq = "";
		return ($indicator, $seq);
	}
	
	# $file is the path to the chromosome fasta file.
	my $file = $path2hg19."$chr.fa";
	my $hg19_in; # Initialize the chromosome file handle.
	if(! open($hg19_in,$file)) { # End the program if a chromosome file can't be found.
	    print 'ERROR opening file $file in ISAdataReformat::getseq_fast.'."\n";
	    exit(1);
	}
	
	# Use the 30 bases upstream and downstream of the integration site.
	my ($start, $stop);
	$start = $IS_base-30;
	$stop = $IS_base+30;
	
	# The start coordinate must be smaller than the stop coordinate.
	if ($stop < $start) { die 'Coodinates wrong in ISAdataReformat::getseq_fast.'."\n";
	}
	# Verify that the file is in fasta format.
	my $title=<$hg19_in>;
	if($title!~/^>/) {
	    print 'Offset bogus or not fasta file in ISAdataReformat::getseq_fast.'."\n";
	    exit(1);
	}
	
	my $headeroffset=length($title);
	chomp($title);
	my $line=<$hg19_in>;
	my $basesperline=length($line)-1;
	
	# compute various offsets required...
	my $lfstostart=int(($start/$basesperline)); #how many lfeed bytes
	my $lfswithin=int($stop/$basesperline)-$lfstostart; #line feed bytes within sequence
	my $startbyte=$headeroffset+$lfstostart+$start; #fasta header + linefeeder + start
	my $bytestoread=$stop+$lfswithin-$start;
	
	# go get the sequence...
	seek($hg19_in,$startbyte,0);
	my $seq;
	my $basesread=read($hg19_in,$seq,$bytestoread);
	$seq=~s/\s//g;
	$seq = uc($seq);
		
	if($basesread != $bytestoread) {
	    print '$chr, $strand, $IS_base: ERROR in read in ISAdataReformat::getseq_fast'."\n";
	    close($hg19_in);
	    exit(1);
	}
	
	$sequence_ID =~ s/\#(fuzz\?{0,1})//;
	my $indicator;
#	No need to check chr17+3618824. Indicate that it's the control and return.
	if ($chr =~ m/chr17/i and $IS_base > 3618822 and $IS_base < 3618826) {
		$indicator = "$virus_strand\tSpike-in control";
	}else {
		# Get the reverse complement if the integration site is on the minus strand.
		$seq = rev_com($seq) if $strand eq "-";
		# Determine whether mis-priming is likely.
		$seq = uc($seq);
		$indicator = compare_to_primer($seq, $primer_bases,$virus_strand);
		$indicator .= "\t" unless $indicator =~ m/\w/;
 	}
	close($hg19_in);
	# Return the sequence and the prediction.

	return($indicator, $seq);
	
}

##########################################################################################################################

# This subroutine looks for matches between two DNA sequences as part of the mis-prime
# prediction step. 

sub compare_to_primer {

	# Initialize values.
	my($IS_bases, $primer_bases, $virus_strand) = @_;
	
	# Extract the 10 bases to compare.
	my @IS_bases = split //, $IS_bases;
	
	my $score = 0;
	my $indicator = "$virus_strand\tno_indicator";
	
	# Starting 30 bases upstream from the integration site, look for matches between the 
	# primer bases and the genomic sequence. If there are 7 or more matching base (70% match),
	# report that mis-priming is likely. If no match is found, move downstream 1 base and
	# search again for matches.  Repeat until a 70 % match is found or the last primer base
	# reaches the inegration site.
	
	SCAN: for (my $i = 0; $i <= 21; $i++) {
		foreach (0..9) {
			$score += 1 if $IS_bases[$_] =~ /$$primer_bases[$_]/;
		}
		if ($score >= 7) {
			$indicator = "$virus_strand\tpossible mis-prime ($score"."0% match)";
			last SCAN;
		}
		shift @IS_bases;
		$score = 0;
	}
	# Return the prediction.
	return $indicator;
}

##########################################################################################################################

# Subroutine to generate and return the reverse complement of the sequence it receives.

sub rev_com {
	my $seq = shift;
	my $revcom = reverse $seq;
	$revcom =~ tr/ATGC/TACG/;
	return $revcom;
}

##########################################################################################################################

# A subroutine that returns an array of the bases of the sequence it receives, swapping out
# IUPAC codes for expressions to use in a regex.

sub split_primer {
	my %IUPAC = (
	"R"	=> "[AG]",
	"Y"	=> "[CT]",
	"S"	=> "[GC]",
	"W"	=> "[AT]",
	"K"	=> "[GT]",
	"M"	=> "[AC]",
	"B"	=> "[CGT]",
	"D"	=> "[AGT]",
	"H"	=> "[ACT]",
	"V"	=> "[ACG]",
	);
	my @bases =  split //, $_[0];
	foreach (@bases) {
		$_ = $IUPAC{$_} if exists $IUPAC{$_};
	}
	return @bases;
}

##########################################################################################################################

# This subroutine uses sequence ID's to pull the raw sequences from the debarcoded data files.

sub extract_sequences {
	my ($blat_lines, $file) = @_;
	
	my (%data_lines);
	foreach my $line (@{$blat_lines}) {
		chomp($line);
		my @elements = split /\t/, $line;
		my $seq_id = $elements[1];
 		$seq_id =~ s/\#fuzz\?{0,1}//;
		$data_lines{$seq_id} = [@elements];
	}
	
	open(my $fh, "gzip -dc $file |") or die "Unable to find $file: $!\n";
	while (<$fh>) {
		chomp;
		my ($seq_id,$seq1pre,$seq1post,$seq2pre,$seq2post) = (split /\t/)[0,2,3,7,8];
		$seq_id =~ s/#[ATGCN]{10}//;
		# print join("\t", $seq_id, $seq1pre, $seq1post) . "\n";
		if (exists $data_lines{$seq_id}) {
			next if $data_lines{$seq_id}[3] =~ m/HIV-1|HIV1/;
			my $seq1 = $seq1pre.$seq1post;
			my $seq2 = $seq2pre.$seq2post;
			push @{$data_lines{$seq_id}},"$seq1\t$seq2";
		}
	}
	close $fh;
	return %data_lines;
}

# This subroutine sorts the Read1 fasta file by length for Uclust.

sub sort_fasta_by_length {
	
	my $folder = shift;
	# Initialize folder and file handle variables.
	open my $data_in, "<","$folder/PEread_LTR_nonInternal_linker_long_UniquePairWct_F.fa" or die "$folder/PEread_LTR_nonInternal_linker_long_UniquePairWct_F.fa: $!\n";
	open my $data_out, ">$folder/Read1_insert.fa" or die;
	
	my @sites; # @sites stores the data coming in.
	while (my $line=<$data_in>) { # Iterate through $data_in.
		chomp($line);
		if ($line=~/^>/) { # Sequence ID's begin with ">".
			$line =~ s/\#.*\#.*//;
			my $id=$line;
			my $seq=<$data_in>; # The next line is the sequence.
			chomp($seq);
			my $seqlen=length($seq); # Obtain the length of the sequence.
			# Store the sequence information as a hash reference.
			my $site = {
					"ID" => $id,
					"SEQ" => $seq,
					"LEN" => $seqlen,
			};
			# Only keep the sequences 20 bases or longer.
			if($seqlen>19) {
				push (@sites, $site);
			}
		}
	}
	# Sort the sites by descending length and then lexically by sequence.
	my @sortedsites = sort {
		$b->{LEN} <=> $a->{LEN}
				||
		$a->{SEQ} cmp $b->{SEQ} 
	} @sites;
	# Print the data.
	foreach my $site (@sortedsites) {
		print $data_out "$site->{ID}\n$site->{SEQ}\n";
	}
	
	close $data_in;
	close $data_out;
	
	return;
}

##########################################################################################################################

# This subroutine sorts the uClust output into different formats for easy reading.

sub resort_Uclust {

	my ($name,$folder) = @_;
	
	my $uc_file = "$folder/Read1_insert.uc";
	open my $in_uc, "<", $uc_file or die "$uc_file: $!\n";
	
	# Initialize the variables.
	my %seq_id2map; # %seq_id2map allows you to obtain the integration site using the sequence ID.
	my %seq_id2seq; # %seq_id2seq allows you to obtain the sequence using the sequence ID.
	my %cluster_id2breakpoints; # %cluster_id2breakpoints stores the lengths of the sequences associated with that cluster ID.
	my %cluster_id2total_reads; # %cluster_id2total_reads allows you to obtain the total number of sequences
								# in a cluster using the cluster ID.
	my %cluster_id2seq_id; # %cluster_id2seq_id makes a list of the sequence ID's for each sequence in the cluster.
						   # Since the uClust output is sorted by sequence length, so are the sequences in each list.
	my @all_cid; # @all_cid is a list of all the cluster ID's (cid is short for cluster ID).
	my %cid2seed; # %cid2seed matches the cluster ID to the one sequence that was used as the seed for that cluster.
	
	# Iterate through the uClust file and parse the important information. The first portion of the 
	# file contains all of the sequences sorted by length. The lines begin with either S or H and contain
	# information on how that sequence relates to its cluster. The second field of each line is the number 
	# of the cluster which is combines with "C" to make the cluster ID, e.g. "C1234".
	
	<$in_uc> foreach 1..8; # Skip the first 8 lines of the Read1_insert.uc.
	while (my $line=<$in_uc>) {
	
		chomp($line);
		my @eles=split(/\t/, $line); # Split the line by tab.
		my $cluster_id ="C"."$eles[1]"; # Add "C" to the ID because user-defined variables can't start with a number.
		my $seq_id = $eles[8]; # Obtain the sequence ID.
		$seq_id =~ s/\#.*//;

		if($line=~/^(S|H)/) { # For the lines starting with S or H.
			
			my $breakpoint=$eles[2]; # The length of the sequence.
			# Create a hash for each cluster ID where each value is a hash that uses the breakpoint as a key.
			# This is used to count the number of breakpoints in each cluster.
			$cluster_id2breakpoints{$cluster_id}{$breakpoint} = 1 unless exists $cluster_id2breakpoints{$cluster_id}{$breakpoint};
			# Initialize an array for $cluster_id2seq_id{$cluster_id} if this is the first time that this cluster ID has been seen.
			$cluster_id2seq_id{$cluster_id} = [] unless exists $cluster_id2seq_id{$cluster_id};
			# Store each sequence ID belonging to each cluster.
			push @{$cluster_id2seq_id{$cluster_id}}, $seq_id;
			
		# The second portion of the uClust file contains cluster-specific information. Each line begins with C.
		# The third field in each line contains the number of sequences in that cluster. The cluster ID's are
		# in numerical order.
		}elsif($line=~/^C/) { # For lines starting with C.
		
			$cluster_id2total_reads{$cluster_id} = $eles[2]; # Store the sequence total for each cluster.
			push @all_cid, $cluster_id; # Store the cluster ID's on a list retaining their numerical order.
			$cid2seed{$cluster_id} = $seq_id; # Store the sequence ID of the seed sequence for each cluster.
			
		}
	}		
	close $in_uc;
	
	# Pull in the sequences for each sequence ID.
	open my $in_txt, "<","$folder/Read1_insert.fa" or die "$folder/Read1_insert.fa: $!\n";
	while (my $line=<$in_txt>) {
		chomp($line);
		if ($line=~/^>(.*)/) { # Remove ">" from the sequence ID line.
			my $seq_id = $1;
			my $seq = <$in_txt>; # The sequence is next.
			chomp($seq);
			$seq_id2seq{$seq_id} = $seq; # Store the sequence.
		}
	}

	close $in_txt;
	
	# Pull in the integration site and breakpoint coordinates from the _blat file. Not 
	# every sequence ID from the uClust file has coordinates in the _blat file.
	my $blat_file = "$folder/PEread_LTR_NonInternal_linker_long_UniquePairWct_F_R_for_Uclust.txt";
	$blat_file =~ s/\.gz^//;
	open my $blat_in, "<", $blat_file or die "$blat_file: $!\n";
	while (my $line=<$blat_in>) { # Iterate through each line.
		chomp($line);
		my @eles=split(/\t/, $line);
		shift @eles;
		my $seq_id = shift @eles; # Obtain the sequence ID.
		$seq_id =~ s/\#.*//;
		next unless exists $seq_id2seq{$seq_id};
		my $map = join("\t", @eles); # Create one variable for the coordinates.
		$seq_id2map{$seq_id} = $map; # Store the coordinates.
	}
	close $blat_in;
	
	# Store the cluster ID and its number of breakpoints and reads to a hash reference
	# and push that reference onto an array.
 	my @all_clusters; # The array to store the cluster hash references.
 	
 	# Create a hash reference for each cluster.
 	foreach my $cluster_id (keys %cluster_id2total_reads) {
 		# Obtain the number of breakpoints by assigning the list of the breakpoint 
 		# lengths keys %{$cluster_id2breakpoints{$cluster_id}} to a scalar variable.
 		my $breakpoints = keys %{$cluster_id2breakpoints{$cluster_id}};
 		my $total_reads = $cluster_id2total_reads{$cluster_id}; # Obtain the total reads for the cluster.
 		my $cluster_count ={ # Populate the hash.
 				"ID" => "$cluster_id",
 				"CNT" => "$breakpoints",
 			"TotalCNT" => "$total_reads",
 	  	};
 		push(@all_clusters, $cluster_count); # Add the hash reference to the list @all_clusters.
 	}
 	
 	# Sort the clusters for a cleaner data file. The clusters will be presented by descending
 	# number of breakpoints in the cluster.
 	my @all_clusters_sorted = sort {
 		$b->{'CNT'} <=> $a->{'CNT'} # First sort by descending the number of breakpoints.
 				||
 		$b->{'TotalCNT'} <=> $a->{'TotalCNT'} # Then sort descending by the number of sequences.
 				||
 		$a->{'ID'} cmp $b->{'ID'}  # Then sort by cluster ID.
 	} @all_clusters;
 	
 	
 	# Print the sorted clusters to a file.
	open my $out_SeqMap, ">", "$folder/Uclust_results_reformatted_with_sequences.txt" or die "Unable to create $folder/Uclust_results_reformatted_with_sequences.txt: $!\n";
 	foreach my $cluster (@all_clusters_sorted) { # Perform for each cluster.
 		my $cluster_id = $cluster->{'ID'}; # Assign the cluster ID to a new scalar variable.
 		my $breakpoints = $cluster->{'CNT'}; # Assign the number of breakpoints to a new scalar variable.
 		my $total_reads = $cluster->{'TotalCNT'}; # Assign the total number of reads to a new scalar variable.
 		# Each sequence in the cluster will be on its own line. The sequences are in descending order by length.
  		foreach my $seq_id (@{$cluster_id2seq_id{$cluster_id}}) { 
  			# Print cluster information, the sequence, and the integration site coordinates if available.
  			my $line_out = "$cluster_id\t$breakpoints\t$total_reads\t$seq_id\t$seq_id2seq{$seq_id}";
  			$line_out .= "\t$seq_id2map{$seq_id}" if exists $seq_id2map{$seq_id};
			print $out_SeqMap "$line_out\n";
  		}
 	}
 	close $out_SeqMap;

 	# Print a file containing only the sequence ID of the seed for each cluster.
 	open my $out_ClusterBreak, ">","$folder/Uclust_list_of_cluster_seeds.txt" or die "Unable to create $folder/Uclust_list_of_cluster_seeds.txt: $!\n";
 	# Iterate through the list of cluster ID's in numercal order.
 	foreach my $cid (@all_cid) { 
 		# Obtain the total number of breakpoints in the current cluster.
 		my $count = keys %{$cluster_id2breakpoints{$cid}};
 		# Print the cluster ID, counts, and seed sequence ID.
 		print $out_ClusterBreak "$cid\t$count\t$cid2seed{$cid}\n";
 	}
 	close $out_ClusterBreak;

 	
 	return;	
}

1;

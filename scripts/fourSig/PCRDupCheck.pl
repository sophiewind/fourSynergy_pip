#! /usr/bin/perl

use strict;
use warnings;

# This script will go through the original paired end sequencing file and pull out the bar codes originally designed
# in the 4C primers. These are at a fixed length, set by the $barSize variable. It then goes through and finds all of
# the barcodes associated with each sequence read pair. A count of all bar codes is kept and printed.
#
# #### Be sure to check how fw and rv reads are identified in the ID. Sometimes it's something like 1:N:0 and 2:N:0, but
# this will be different based on other barcodes, etc. Make sure to change this on lines 67, 94, 185, and 214. ####
#
# If 4C primers were multiplexed, the program has the ability to scan for reads associated with individual 4C datasets.
# This is done by inputing a comma-separated list of anchor files and associated sequencing read files. Only the fw
# files are needed for this. If not multiplexed, just enter one anchor and its associated data. Since this is for
# allele-specific anchors, each anchor name will have two sequencing read files associated with it. The two files are
# separated by a colon. So the anchor list would look like this:
#
# Anchor1,Anchor2,etc...
#
# The anchor file list would look like this:
# 
# Anchor1_Strain1_seq_fw:Anchor1_Strain2_seq_fw,Anchor2_Strain1_seq_fw:Anchor2_Strain2_seq_fw, etc...
#
# Finally, you need a PATH/basename for your output files.
#
# Output (large sequencing files):
#
# Basename_largeFile_Counts.txt - A list with all of the barcodes and how many times they occurred in the large dataset.
# Basename_largeFile_Histogram.txt - A list with the number of duplicates and how many times that number of duplicates 
#	occurred. Useful for a histogram of how many times a potential bias could have occurred.
#
# Output (per anchor file):
#
# Basename_AnchorName_Strain_1_BarCodeIDs.txt - A list of Sequencing IDs associated with each barcode for strain 1.
# Basename_AnchorName_Strain_2_BarCodeIDs.txt - A list of Sequencing IDs associated with each barcode for strain 2.
# Basename_AnchorName_Strain_1_BarCodeCounts.txt - A list of counts for each barcode in strain 1.
# Basename_AnchorName_Strain_2_BarCodeCounts.txt - A list of counts for each barcode in strain 2.
# Basename_AnchorName_Strain_1_BarCodeHist.txt - A list with the number of duplicates and how many times that number of #  duplicates occurred in strain 1. It is identical to the one for the large file, but specific to the anchor/strain 
#  combination.
# Basename_AnchorName_Strain_2_BarCodeHist.txt - A list with the number of duplicates and how many times that number of #  duplicates occurred in strain 2. It is identical to the one for the large file, but specific to the anchor/strain 
#  combination.
#
# Example:
#
# ./PCRDupCheck.pl PATH/fw_fastq_file PATH/rv_fastq_file barCodeSize(in bp) anchorNameList(Comma_separated) listOfFiles(per anchor,comma-separated, include the PATH) PATH/OutputBase
#

my %largeData;
my %fourmerCounts;

my ($fwLargeFile, $rvLargeFile, $barSize, $anchorNameList, $anchorFileList, $bigNameBase) = @ARGV;

print STDERR "Opening the forward file.\n";

open (FWIN, "<".$fwLargeFile) or die "Cannot open ".$fwLargeFile."!\n";

while (my $fwID = <FWIN>){
	my $fwSeq = <FWIN>;
	my $fwSymbol = <FWIN>;
	my $fwPhred = <FWIN>;
	
	chomp($fwID);
	
	my $fwIDKey = $fwID;
	# Depending upon the version of sequencing software, you might have to change what the REGEX searches for.
	$fwIDKey =~ s/\s1:N:0://;
	
	my $fwBarCode = substr($fwSeq, 0, $barSize);
	
	if (!defined($largeData{$fwIDKey})){
		$largeData{$fwIDKey} = $fwBarCode;
	} else {
		print STDERR "The ID ".$fwIDKey." is duplicated!\n";
		exit;
	}
}

close(FWIN);

print STDERR "Finished with the forward file.\nOpening the reverse file.\n";

open (RVIN, "<".$rvLargeFile) or die "Cannot open ".$rvLargeFile."!\n";

while (my $rvID = <RVIN>){
	my $rvSeq = <RVIN>;
	my $rvSymbol = <RVIN>;
	my $rvPhred = <RVIN>;
	
	chomp($rvID);
	
	my $rvIDKey = $rvID;
	# Depending upon the version of sequencing software, you might have to change what the REGEX searches for.
	$rvIDKey =~ s/\s2:N:0://;
	
	my $rvBarCode = substr($rvSeq, 0, $barSize);
	
	if (!defined($largeData{$rvIDKey})){
		print STDERR "The ID ".$rvIDKey." is not present in the forward file!\n";
		exit;
	} else {
		if (!defined($fourmerCounts{$largeData{$rvIDKey}."_".$rvBarCode})){
			$fourmerCounts{$largeData{$rvIDKey}."_".$rvBarCode} = 1;
		} else {
			$fourmerCounts{$largeData{$rvIDKey}."_".$rvBarCode}++;
		}
		$largeData{$rvIDKey} .= "_".$rvBarCode;
	}
}

close(RVIN);

print STDERR "Finished with the reverse file.\nOpening and printing the count data.\n";

open (COUNTOUT, ">".$bigNameBase."_largeFile_Counts.txt") or die "Cannot open ".$bigNameBase."_largeFile_Counts.txt!\n";

print COUNTOUT "Barcode\tOccurrence\n";

my $lowestBigNumber = 100000000000;
my $highestBigNumber = 0;

my %bigCountHist;

foreach my $printCount (sort keys (%fourmerCounts)){
	print COUNTOUT $printCount."\t".$fourmerCounts{$printCount}."\n";
	if ($fourmerCounts{$printCount}<$lowestBigNumber){
		$lowestBigNumber = $fourmerCounts{$printCount};
	}
	if ($fourmerCounts{$printCount}>$highestBigNumber){
		$highestBigNumber = $fourmerCounts{$printCount};
	}
	
	if (!defined($bigCountHist{$fourmerCounts{$printCount}})){
		$bigCountHist{$fourmerCounts{$printCount}} = 1;
	} else {
		$bigCountHist{$fourmerCounts{$printCount}}++;
	}
}

close (COUNTOUT);

open (HISTOUT, ">".$bigNameBase."_largeFile_Histogram.txt") or die "Cannot open ".$bigNameBase."_largeFile_Histogram.txt!\n";

print HISTOUT "NumberDuplicates\tOccurrence\n";

foreach my $histCount (sort {$a<=>$b} keys (%bigCountHist)){
	print HISTOUT $histCount."\t".$bigCountHist{$histCount}."\n";
}

close (HISTOUT);

print STDERR "The largest number of duplicates for any given barcode was ".$highestBigNumber."\n".
			 "The lowest number of duplicates for any given barcode was ".$lowestBigNumber."\n\n";

my @anchorNames = split(/,/, $anchorNameList);
my @anchorFileSets = split(/:/, $anchorFileList);

my $anchorNameListPosition = 0;

foreach my $anchorFileSet (@anchorFileSets){
	
	my %anchorBBarCodes;
	my %anchorBBarCodeCounts;
	my %anchorCBarCodes;
	my %anchorCBarCodeCounts;
	
	my $analyzedAnchor = $anchorNames[$anchorNameListPosition];
	$anchorNameListPosition++;
	
	print STDERR "Analyzing files for ".$analyzedAnchor.".\n";
	
	my ($anchorBFwFile, $anchorCFwFile) = split(/,/, $anchorFileSet);
	
	open (AFWINB, "<".$anchorBFwFile) or die "Cannot open ".$anchorBFwFile."!\n";
	
	while (my $fwBID = <AFWINB>){
		my $fwBSeq = <AFWINB>;
		my $fwBSymbol = <AFWINB>;
		my $fwBPhred = <AFWINB>;
		
		chomp($fwBID);
		
		my $fwBcheckID = $fwBID;
		# Depending upon the version of sequencing software, you might have to change what the REGEX searches for.
		$fwBcheckID =~ s/\s1:N:0://;
		
		if (!defined($largeData{$fwBcheckID})){
			print STDERR "The ID ".$fwBcheckID." was not found in the large files!\n";
			exit;
		} else {
			if (!defined($anchorBBarCodeCounts{$largeData{$fwBcheckID}})){
				$anchorBBarCodeCounts{$largeData{$fwBcheckID}} = 1;
				$anchorBBarCodes{$largeData{$fwBcheckID}} = $fwBcheckID;
			} else {
				$anchorBBarCodeCounts{$largeData{$fwBcheckID}}++;
				$anchorBBarCodes{$largeData{$fwBcheckID}} .= ",".$fwBcheckID;
			}
		}
	}
	
	close (AFWINB);
	
	open (AFWINC, "<".$anchorCFwFile) or die "Cannot open ".$anchorCFwFile."!\n";
	
	while (my $fwCID = <AFWINC>){
		my $fwCSeq = <AFWINC>;
		my $fwCSymbol = <AFWINC>;
		my $fwCPhred = <AFWINC>;
		
		chomp($fwCID);
		
		my $fwCcheckID = $fwCID;
		# Depending upon the version of sequencing software, you might have to change what the REGEX searches for.
		$fwCcheckID =~ s/\s1:N:0://;
		
		if (!defined($largeData{$fwCcheckID})){
			print STDERR "The ID ".$fwCcheckID." was not found in the large files!\n";
			exit;
		} else {
			if (!defined($anchorCBarCodeCounts{$largeData{$fwCcheckID}})){
				$anchorCBarCodeCounts{$largeData{$fwCcheckID}} = 1;
				$anchorCBarCodes{$largeData{$fwCcheckID}} = $fwCcheckID;
			} else {
				$anchorCBarCodeCounts{$largeData{$fwCcheckID}}++;
				$anchorCBarCodes{$largeData{$fwCcheckID}} .= ",".$fwCcheckID;
			}
		}
	}
	
	close (AFWINC);
	
	my %BDuplicates;
	my %CDuplicates;
	
	open (AFOUTB, ">".$bigNameBase."_".$analyzedAnchor."_Strain_1_BarCodeIDs.txt") or die "Cannot open ".$bigNameBase."_".$analyzedAnchor."_Strain_1_BarCodeIDs.txt!\n";
	open (AFOUTCOUNTB, ">".$bigNameBase."_".$analyzedAnchor."_Strain_1_BarCodeCounts.txt") or die "Cannot open ".$bigNameBase."_".$analyzedAnchor."_Strain_1_BarCodeCounts.txt!\n";
	
	print AFOUTB "BarCode\tIDs\n";
	print AFOUTCOUNTB "BarCode\tCounts\n";
	
	my $smallestB = 10000000000;
	my $largestB = 0;
	
	foreach my $barCodeB (sort keys (%anchorBBarCodes)){
		print AFOUTB $barCodeB."\t".$anchorBBarCodes{$barCodeB}."\n";
		print AFOUTCOUNTB $barCodeB."\t".$anchorBBarCodeCounts{$barCodeB}."\n";
		if ($anchorBBarCodeCounts{$barCodeB} < $smallestB){
			$smallestB = $anchorBBarCodeCounts{$barCodeB};
		}
		if ($anchorBBarCodeCounts{$barCodeB} > $largestB){
			$largestB = $anchorBBarCodeCounts{$barCodeB};
		}
		if (!defined($BDuplicates{$anchorBBarCodeCounts{$barCodeB}})){
			$BDuplicates{$anchorBBarCodeCounts{$barCodeB}} = 1;
		} else {
			$BDuplicates{$anchorBBarCodeCounts{$barCodeB}}++;
		}
	}
	
	close (AFOUTB);
	close (AFOUTCOUNTB);
	
	open (AFOUTC, ">".$bigNameBase."_".$analyzedAnchor."_Strain_2_BarCodeIDs.txt") or die "Cannot open ".$bigNameBase."_".$analyzedAnchor."_Strain_2_BarCodeIDs.txt!\n";
	open (AFOUTCOUNTC, ">".$bigNameBase."_".$analyzedAnchor."_Strain_2_BarCodeCounts.txt") or die "Cannot open ".$bigNameBase."_".$analyzedAnchor."_Strain_2_BarCodeCounts.txt!\n";
	
	print AFOUTC "BarCode\tIDs\n";
	print AFOUTCOUNTC "BarCode\tCounts\n";
	
	my $smallestC = 10000000000;
	my $largestC = 0;
	
	foreach my $barCodeC (sort keys (%anchorCBarCodes)){
		print AFOUTC $barCodeC."\t".$anchorCBarCodes{$barCodeC}."\n";
		print AFOUTCOUNTC $barCodeC."\t".$anchorCBarCodeCounts{$barCodeC}."\n";
		if ($anchorCBarCodeCounts{$barCodeC} < $smallestC){
			$smallestC = $anchorCBarCodeCounts{$barCodeC};
		}
		if ($anchorCBarCodeCounts{$barCodeC} > $largestC){
			$largestC = $anchorCBarCodeCounts{$barCodeC};
		}
		if (!defined($CDuplicates{$anchorCBarCodeCounts{$barCodeC}})){
			$CDuplicates{$anchorCBarCodeCounts{$barCodeC}} = 1;
		} else {
			$CDuplicates{$anchorCBarCodeCounts{$barCodeC}}++;
		}
	}
	
	close (AFOUTC);
	close (AFOUTCOUNTC);
	
	open (OUTBHIST, ">".$bigNameBase."_".$analyzedAnchor."_Strain_1_BarCodeHist.txt") or die "Cannot open ".$bigNameBase."_".$analyzedAnchor."_Strain_1_BarCodeHist.txt!\n";
	
	print OUTBHIST "NumberDuplicates\tOccurrence\n";
	
	foreach my $printDuplicateB (sort {$a <=> $b} keys (%BDuplicates)){
		print OUTBHIST $printDuplicateB."\t".$BDuplicates{$printDuplicateB}."\n";
	}	
	
	close (OUTBHIST);
	
	open (OUTCHIST, ">".$bigNameBase."_".$analyzedAnchor."_Strain_2_BarCodeHist.txt") or die "Cannot open ".$bigNameBase."_".$analyzedAnchor."_Strain_2_BarCodeHist.txt!\n";
	
	print OUTCHIST "NumberDuplicates\tOccurrence\n";
	
	foreach my $printDuplicateC (sort {$a <=> $b} keys (%CDuplicates)){
		print OUTCHIST $printDuplicateC."\t".$CDuplicates{$printDuplicateC}."\n";
	}	
	
	close (OUTCHIST);
	
	print STDERR "The largest number of duplicates for any given fourmer pair in the Strain_1 data was ".
				 $largestB."\n".
			 	 "The lowest number of duplicates for any given fourmer pair in the Strain_1 data was ".
			 	 $smallestB."\n".
			 	 "The largest number of duplicates for any given fourmer pair in the Strain_2 data was ".
				 $largestC."\n".
			 	 "The lowest number of duplicates for any given fourmer pair in the Strain_2 data was ".
			 	 $smallestC."\n\n";

}	

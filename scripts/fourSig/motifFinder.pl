#!/usr/bin/perl

##########
#
# IF YOU ARE USING A GENOME OTHER THAN mouse OR human... scroll down until
# you see the words "EDIT ME!!!!" and follow the instructions there.
# 
##########

use strict;
use warnings;

use Getopt::Std;

use constant TRUE => 1;
use constant FALSE => 0;

my $sHelp = <<END_OF_HELP;

Usage: motifFinder.pl [options] MOTIF genomeCHRNUM.fa > OUTPUT.txt

Where:

MOTIF = the motif you are are looking for.  For example, if you want to find
  all of the HindIII sites, specify "AAGCTT"

genomeCHRNUM.fa = This a special way to specify the names of the FASTA files 
  that have the sequence for the genome you want to search for motifs in. You
  should have one file per chromosome.  The "genome" part should be replaced 
  with the path and the first part of the name of the files.  Leave "CHRNUM"
  in, as that part will be replaced by a chromosome number.
  For example, if the wanted to find motifs in the hg19 genome, we might
  specify "hg19/Sequence/Chromosomes/chrCHRNUM.fa", if, in
  "hg19/Sequence/Chromosomes/", we had files with the names, "chr1.fa", 
  "chr2.fa", "chr3.fa", etc.

Here is an example of how to identify HindIII motifs in the human genome:

./motifFinder.pl -H AAGCTT PATH/TO/hg19/chrCHRNUM.fa > hg19_hind3_sites.txt


Options

-h
    Print this help information.

-M
    The data is from the mouse genome.
    NOTE: If you are using something other than mouse or human, you can change 
    the code pretty easily to fit your needs.

-H
    The data is from the human genome.
    NOTE: If you are using something other than mouse or human, you can change 
    the code pretty easily to fit your needs.

-o  CHR_NAMES
    Omit the chromosomes listed in CHR_NAMES, which is a list,
    ie. "Y", or "1 2" etc.
    

END_OF_HELP

if (-t STDIN && !@ARGV) {
    print STDERR $sHelp;
    exit 1;
}


# process the command line arguments
my %opts; # a hash table to store file names passed in as agruments
getopts('hMHo:', \%opts);

if ($opts{'h'}) { # print help and exit
    print STDERR $sHelp;
    exit 1;
}

##########
#
# EDIT ME!!!!
#
# YOU NEED TO EDIT THE LIST OF CHROMOSOMES IN @chromosomes (default)
#
##########

## EDIT ME!!!!
my @chromosomes = ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y");

################
#
# END OF SECTION YOU NEED TO EDIT IF YOU WANT TO USE A GENOME OTHER THAN mouse
#
# OR human
#
################

if ($opts{'M'}) { # use mouse chromosomes (mm9)
    print STDERR "Using the mouse genome\n";
    @chromosomes = ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y");
} elsif ($opts{'H'}) {
    print STDERR "Using the human genome\n";    
    @chromosomes = ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20","21","22","X", "Y");

} else {
    print STDERR "Using default genome, initialy set to mouse, but you may have changed it.\n";
}

print STDERR "\nThese are the chromosomes that will be examined...\n";
foreach my $chr (@chromosomes) {
    print STDERR "\t".$chr."\n";
}
print STDERR "\n";


if (defined($opts{'o'})) {
    my @chrs = split(/\s+/, $opts{'o'});
    foreach my $chr (@chrs) {
	print STDERR "Omitting $chr from the analysis...\n";
	for(my $index = 0; $index < scalar(@chromosomes); $index++) {
	    if ($chromosomes[$index] eq $chr) {
		splice(@chromosomes, $index, 1);
	    }
	}
    }

    print STDERR "This is the new list of chromosomes...\n";
    foreach my $chr (@chromosomes) {
	print STDERR $chr."\n";
    }
}


my $isPalindrome = FALSE;


my $motif_seq = shift(@ARGV);
$motif_seq =~ s/N/\./g;
$motif_seq = lc($motif_seq);

my $rev_seq = reverse($motif_seq);
$rev_seq =~ tr/atcg/tagc/;

if ($motif_seq eq $rev_seq) {
    $isPalindrome = TRUE;
    print STDERR "The motif: ".$motif_seq." is a palindrome.\n";
} else {
    print STDERR "The motif: ".$motif_seq." is not a palindrome.\n";
    print STDERR "Looking for: ".$motif_seq." and ".$rev_seq."\n";
}

my $motif_length = length($motif_seq);

my $seqFile = shift(@ARGV);
chomp($seqFile);

#for (my $i=1; $i <= 20; $i++) {
foreach my $chr (@chromosomes) {

    my $chromosome_file = $seqFile;
    $chromosome_file =~ s/CHRNUM/$chr/g;

    print STDERR "looking for sites in: ".$chromosome_file."\n";

    open(CHRM_FILE, "<".$chromosome_file) 
	|| die("Could not open ".$chromosome_file."\n".$!."\n");

    my $comment = <CHRM_FILE>;

    my $seq = ""; #<CHRM_FILE>;
    my $position = 0;
    while(my $line = <CHRM_FILE>) {
	chomp($line);
	$seq .= lc($line);
    }

    my $seq_length = length($seq);

    for(my $j=0; $j<$seq_length-$motif_length; $j++) {
	my $genomeFragment = substr($seq, $j, $motif_length);
	if ($genomeFragment =~ /$motif_seq/) {
	    print $chr."\t".$j."\n";
	}
	if (!$isPalindrome) {
	    if ($genomeFragment =~ /$rev_seq/) {
		print $chr."\t".$j."\n";
	    }
	}
    }

}

#!/usr/bin/perl

##########
#
# IF YOU ARE USING A GENOME OTHER THAN mm9 OR hg19... scroll down until
# you see the words "EDIT ME!!!!" and follow the instructions there.
# 
##########


use strict;
use warnings;

use POSIX qw(floor);

use Getopt::Std;

use constant TRUE => 1;
use constant FALSE => 0;

use constant SAM_CHR => 2;
use constant SAM_POS => 3;

use constant MIN_READ_LENGTH => 20; 

use constant BIN_SIZE => 100000;

sub is_whole_number { 
    $_[0] =~ /^\d+$/
}


my $sHelp = <<END_OF_HELP;

Usage: samToReTab.pl [options] LENGTH_A MAX_SIZE MIN_SIZE SAM_File.sam 3C_RE_SITES.txt 4C_RE_SITES.txt LI_RE_SITES.txt > OUTPUT.tab

Where:

LENGTH_A = The length of the sequenced read that was in the anchor fragment. 
  This includes the length of the 4C primer that you used as well as any 
  additional basepairs that were sequenced prior to the 3C or 4C cut sites. In 
  this case, the value 18 was inferred from the fact that the reads were 76 bp 
  long, but only 58 bp were aligned.

MAX_SIZE = The maximum number of bp submitted for sequencing after the size 
  selection step in 4C.

MIN_SIZE = The minimum number of bp submitted for sequencing.

SAM_FILE.sam = The SAM file.

3C_RE_SITES.txt = A file with locations of the 3C restriction enzyme. Note, 
  before generating one of these files on your own, make sure there is not one 
  already in the re_sites directory. These were created using the mm9 mouse 
  sequence.

4C_RE_SITES.txt = A file with locations of the 4C restriction enzyme.

LI_RE_SITES.txt =  A file with locations of the linearization restriction 
  enzyme, or NONE, if there was no linearization step in the 4C protocol.

OUTPUT.tab = The name of the file you want the output to be redirected into.
  If you do not specify OUTPUT.tab, the output will be printed to STDOUT.

Options

-h
    Print this help information.

-M
    The data is from the mm9 mouse genome. THIS IS THE DEFAULT
    NOTE: If you are using something other than mm9 or hg19, you can change 
    the code pretty easily to fit your needs.

-H
    The data is from the hg19 human genome.
    NOTE: If you are using something other than mm9 or hg19, you can change 
    the code pretty easily to fit your needs.

-o  CHR_NAMES
    Omit the chromosomes listed in CHR_NAMES, which is a list,
    ie. "Y", or "1 2" etc.

-c  CHR_STRING
    By default, samToReTab assumes that each chromosome name is formated like
    this "chr1", where the chromosome number, in this case 1, is preceded with
    the string "chr". If your chromosome names have a different format, like
    "1GRCm38" (where the chromosome number, 1, is followed by "GRCm38"), you 
    can use this option to specify the alternative to the "chr" string.

END_OF_HELP

if (-t STDIN && !@ARGV) {
    print STDERR $sHelp;
    exit 1;
}


# process the command line arguments
my %opts; # a hash table to store file names passed in as agruments
getopts('hMHo:c:', \%opts);

if ($opts{'h'}) { # print help and exit
    print STDERR $sHelp;
    exit 1;
}

##########
#
# EDIT ME!!!!
#
# IF YOU ARE USING A GENOME OTHER THAN mm9 OR hg19, EDIT THE CHROMOSOME
# SIZES BELOW, ADDING/SUBTRACTING CHROMOSOMES AS NEEDED
#
# YOU ALSO NEED TO EDIT THE LIST OF CHROMOSOMES IN @chromosomes
#
##########

## EDIT ME!!!!
my %gChrSizes = ( 
    1 => 197195432, 	
    2 => 181748087, 
    3 => 159599783, 
    4 => 155630120, 
    5 => 152537259, 
    6 => 149517037, 
    7 => 152524553,
    8 => 131738871, 
    9 => 124076172, 
    10 => 129993255, 
    11 => 121843856, 
    12 => 121257530, 
    13 => 120284312, 
    14 => 125194864, 
    15 => 103494974, 
    16 => 98319150, 
    17 => 95272651, 
    18 => 90772031, 
    19 => 61342430, 
    X => 166650296, 
    Y => 15902555
);

## EDIT ME!!!!
my @chromosomes = ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y");

################
#
# END OF SECTION YOU NEED TO EDIT IF YOU WANT TO USE A GENOME OTHER THAN mm9 AND
# hg19
#
################

if ($opts{'M'}) { # use mouse chromosomes (mm9)
    print STDERR "Using mm9 mouse genome\n";
    my %gChrSizes = ( 
	1 => 197195432, 	
	2 => 181748087, 
	3 => 159599783, 
	4 => 155630120, 
	5 => 152537259, 
	6 => 149517037, 
	7 => 152524553,
	8 => 131738871, 
	9 => 124076172, 
	10 => 129993255, 
	11 => 121843856, 
	12 => 121257530, 
	13 => 120284312, 
	14 => 125194864, 
	15 => 103494974, 
	16 => 98319150, 
	17 => 95272651, 
	18 => 90772031, 
	19 => 61342430, 
	X => 166650296, 
	Y => 15902555
	);
    my @chromosomes = ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y");
} elsif ($opts{'H'}) {
    print STDERR "Using hg19 human genome\n";
    %gChrSizes = ( # these are for hg19
	1 => 249250621, 
	2 => 243199373, 
	3 => 198022430, 
	4 => 191154276, 
	5 => 180915260, 
	6 => 171115067, 
	7 => 159138663,
	8 => 146364022, 
	9 => 141213431, 
	10 => 135534747, 
	11 => 135006516, 
	12 => 133851895, 
	13 => 115169878, 
	14 => 107349540, 
	15 => 102531392, 
	16 => 90354753, 
	17 => 81195210, 
	18 => 78077248, 
	19 => 59128983, 
	20 => 63025520,
	21 => 48129895,
	22 => 51304566,
	X => 155270560, 
	Y => 59373566
	);
    
    @chromosomes = ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20","21","22","X", "Y");

} else {
    print STDERR "Using default genome, initialy set to mm9, but you may have changed it.\n";
}


if (defined($opts{'o'})) {
    my @chrs = split(/\s+/, $opts{'o'});
    foreach my $chr (@chrs) {
	print STDERR "Omitting $chr from the analysis...\n";
	delete($gChrSizes{$chr});
	for(my $index = 0; $index < scalar(@chromosomes); $index++) {
	    if ($chromosomes[$index] eq $chr) {
		splice(@chromosomes, $index, 1);
		#delete $chromosomes[$index];
	    }
	}
    }
}

my $chrString = "chr";
if (defined($opts{'c'})) {
    $chrString = $opts{'c'};
}
print STDERR "The chromosome name format string is: ".$chrString."\n";
print STDERR "   Use the -c option to change this.\n";

print STDERR "\nThese are the chromosomes that will be examined...\n";
foreach my $chr (@chromosomes) {
    print STDERR "\t".$chr."\n";
}
print STDERR "\n";

my %gThreeCCutterBins;
my %gThreeCCutterPositions;
my %gFourBaseCutterBins;
my %gSixBaseCutterBins;  # this is the linearization ensyzme
my %gMaxSixBaseBins;
my %gMaxFourBaseBins;

my $knownLength = shift(@ARGV);
if (!is_whole_number($knownLength)) {
    print STDERR "\nOops!  There's an error...\n";
    print STDERR "The first argument, the length, in base pairs, of the anchor fragment that was sequenced needs to be a whole number\n";
    exit();
}
my $maxFragmentLength = shift(@ARGV);
if (!is_whole_number($maxFragmentLength)) {
    print STDERR "\nOops!  There's an error...\n";
    print STDERR "The second argument, the maximum length, in base pairs, of the sequenced fragments, needs to be a whole number\n";
    exit();
}
my $minFragmentLength = shift(@ARGV);
if (!is_whole_number($minFragmentLength)) {
    print STDERR "\nOops!  There's an error...\n";
    print STDERR "The third argument, the minimum length, in base pairs, of the sequenced fragments needs to be a whole number\n";
    exit();
}


my $minUnknownLength = $minFragmentLength - $knownLength;
my $bowtieFile = shift(@ARGV);
my $threeCCutterFile = shift(@ARGV);
my $fourBaseCutterFile = shift(@ARGV);
my $sixBaseCutterFile = shift(@ARGV);

my $noSixBaseCutter = FALSE;
if ($sixBaseCutterFile eq "NONE") {
    print STDERR "Skipping the 6 bp cutter (linearization cutter) analysis \n";
    $noSixBaseCutter = TRUE;
}

# thing to do here is bin the genome and put the threeCCutter sites in bins.
# then figure out which bin a bowtie read fits in and just search that bin
# (and/or the ones +/- 1) for which threeCCutter fragment the read belongs to
print STDERR "Loading in Three C cutter info\n";
open(THREE_C_CUTTER_FILE, "<".$threeCCutterFile) 
    || die("Could not open $threeCCutterFile $!\n");
while(my $line = <THREE_C_CUTTER_FILE>) {
    chomp($line);
    if ($line =~ /\Amotif_seq/) {
	next;
    }
    if ($line =~ /\Arev_seq/) {
	next;
    }
    my ($chr, $pos, $genome) = split(/\t/, $line);
    # strip off 'chr'
    $chr =~ s/chr//;

    # save the position in a hash of arrays (one array per chr)
    push(@{$gThreeCCutterPositions{$chr}}, $pos);

    # figure out what bin to put the threeCCutter site into...
    if(!defined($pos)) {
	print STDERR "Problem reading in the 3C cutter file\n";
	print STDERR "This was an unexpected line...\n";
	print STDERR $line."\n";
	exit;
    }
    my $bin = floor($pos / BIN_SIZE);
    if (defined($gThreeCCutterBins{$chr."\t".$bin})) {
        $gThreeCCutterBins{$chr."\t".$bin} .= "\t".$pos;
    } else {
        $gThreeCCutterBins{$chr."\t".$bin} = $pos;
    }
}
close(THREE_C_CUTTER_FILE);

# This is the 4C enzyme
print STDERR "Loading in Four C cutter (4bp enzyme) info\n";
open(FOUR_BASE_CUTTER_FILE, "<".$fourBaseCutterFile) 
    || die("Could not open $fourBaseCutterFile $!\n");
while(my $line = <FOUR_BASE_CUTTER_FILE>) {
    chomp($line);
    if ($line =~ /\Amotif_seq/) {
	next;
    }
    if ($line =~ /\Arev_seq/) {
	next;
    }
    my ($chr, $pos, $genome) = split(/\t/, $line);
    # strip off "chr"
    $chr =~ s/chr//;

    # figure out what bin to put the fourBaseCutter site into...
    if(!defined($pos)) {
	print STDERR "Problem reading in the 4-base cutter file\n";
	print STDERR "This was an unexpected line...\n";
	print STDERR $line."\n";
	exit;
    }
    my $bin = floor($pos / BIN_SIZE);
    if (defined($gFourBaseCutterBins{$chr."\t".$bin})) {
        $gFourBaseCutterBins{$chr."\t".$bin} .= "\t".$pos;
    } else {
        $gFourBaseCutterBins{$chr."\t".$bin} = $pos;
    }

    my $currentMaxBin = $gMaxFourBaseBins{$chr};
    if (!defined($currentMaxBin)) {
	$gMaxFourBaseBins{$chr} = $bin;
    } elsif ($currentMaxBin < $bin) {
	$gMaxFourBaseBins{$chr} = $bin;
    }
}
close(FOUR_BASE_CUTTER_FILE);


## This is the enzyme that is used to linearize the 4C product.
if (!$noSixBaseCutter) {
    print STDERR "Loading in linearization enzyme (6bp enzyme) info\n";
    open(SIX_BASE_CUTTER_FILE, "<".$sixBaseCutterFile) 
	|| die("Could not open $sixBaseCutterFile $!\n");
    while(my $line = <SIX_BASE_CUTTER_FILE>) {
	chomp($line);
	if ($line =~ /\Amotif_seq/) {
	    next;
	}
	if ($line =~ /\Arev_seq/) {
	    next;
	}
	my ($chr, $pos, $genome) = split(/\t/, $line);
	# strip off "chr"
	$chr =~ s/chr//;
	
	# figure out what bin to put the sixBaseCutter site into...
	my $bin = floor($pos / BIN_SIZE);
	if(!defined($pos)) {
	    print STDERR "Problem reading in the 6-base cutter file\n";
	    print STDERR "This was an unexpected line...\n";
	    print STDERR $line."\n";
	    exit;
	}
	if (defined($gSixBaseCutterBins{$chr."\t".$bin})) {
	    $gSixBaseCutterBins{$chr."\t".$bin} .= "\t".$pos;
	} else {
	    $gSixBaseCutterBins{$chr."\t".$bin} = $pos;
	}
	
	my $currentMaxBin = $gMaxSixBaseBins{$chr};
	if (!defined($currentMaxBin)) {
	    $gMaxSixBaseBins{$chr} = $bin;
	} elsif ($currentMaxBin < $bin) {
	    $gMaxSixBaseBins{$chr} = $bin;
	}
    }
    close(SIX_BASE_CUTTER_FILE);
}

# now bin the reads into 3C fragments (each fragment is indexed by the 
# 3C cutter position on the left side of the fragment)
print STDERR "Figure out the three C fragments for the bowtie reads\n";
open(SAM_FILE, "<".$bowtieFile) || die("Could not open $bowtieFile $!\n");
my %bowtieBins;
while(my $line = <SAM_FILE>) {
    chomp($line);
    
    if ($line =~ /\A\@/) {
        # if the line starts with @, then it is a comment and we don't need it
        next;
    }
    
    my @data = split(/\t/, $line);
    my $bowtieChr = $data[SAM_CHR];
    # strip off $chrString
    $bowtieChr =~ s/$chrString//;

    ## strip off "chr"
    #$bowtieChr =~ s/chr//;

    my $bowtiePos = $data[SAM_POS];

    my $bowtieBin = floor($bowtiePos / BIN_SIZE);
    
    # NOTE: We going to let the 3C site on the left side be the one used to
    # keep track of where the reads are.  Not the 3C site on the right side.
    # So we start with 0 (the edge of the chromosome) as the first 3C site and
    # go from there.
    my $preThreeCCutterSite = 0;
    my $preBin = $bowtieBin - 1;
    my $found = FALSE;
    while ($preBin >= 0 && !$found) {
	# Here's what we are doing...
	# We're looking for the closest 3C cutter site that comes before
	# the current window (so that if none of the 3C cutter sites
	# in thie current window come before the bowtiePos, we know this
	# one must be the one that defines the "left" side of the fragment)
	if (defined($gThreeCCutterBins{$bowtieChr."\t".$preBin})) {
	    my $preThreeCCutterSites = $gThreeCCutterBins{$bowtieChr."\t".$preBin};
	    my @preThreeCCutterSites = split(/\t/, $preThreeCCutterSites);
	    # assuming the 3C sites are ordered from left to right, we are
	    # only interested in the right most site.
	    $preThreeCCutterSite = $preThreeCCutterSites[-1];
	    $found = TRUE; # all done!
	} else {
	    $preBin--; # there were no sites in the previous bin, keep looking!
	}
    }

    # now look in the current bowtiePos
    my @threeCCutterSites;
    if (defined($gThreeCCutterBins{$bowtieChr."\t".$bowtieBin})) {
	my $threeCCutterSites = $gThreeCCutterBins{$bowtieChr."\t".$bowtieBin};
	@threeCCutterSites = split(/\t/, $threeCCutterSites);
    }
    if (scalar(@threeCCutterSites) == 0) {
	# if there are no 3C sites in the bin for bowtiePos, then the read must
	# belong to the fragment started in an earlier bin.
	if (defined($bowtieBins{$bowtieChr."\t".$preThreeCCutterSite})) {
	    $bowtieBins{$bowtieChr."\t".$preThreeCCutterSite}++;
	} else {
	    $bowtieBins{$bowtieChr."\t".$preThreeCCutterSite} = 1;
	}
	next; # all done with this bowtie read!
    }

    # if we're still looking for the 3C fragment, tack the "pre" site onto
    # the list of sites found in the current bin.
    if (defined($preThreeCCutterSite)) {
	@threeCCutterSites = ($preThreeCCutterSite, @threeCCutterSites);
    }

    # look to see if the bowtiePos is between two 3C cut sites in the
    # current bin.
    $found = FALSE;
    for (my $i=1; $i<scalar(@threeCCutterSites); $i++) {
	my $threeCCutterSite = $threeCCutterSites[$i];
	my $prevSite = $threeCCutterSites[($i-1)];
	if ($threeCCutterSite > $bowtiePos) {
	    if (defined($bowtieBins{$bowtieChr."\t".$prevSite})) {
		$bowtieBins{$bowtieChr."\t".$prevSite}++;
	    } else {
		$bowtieBins{$bowtieChr."\t".$prevSite} = 1;
	    }
	    $found = TRUE;
	    last;
	}
    }
    # if the bowtiePos is not between two 3C cut sites, it must come after
    # the last 3C cut site in the bin.
    if (!$found) {
	my $threeCCutterSite = $threeCCutterSites[-1];
	if (defined($bowtieBins{$bowtieChr."\t".$threeCCutterSite})) {
	    $bowtieBins{$bowtieChr."\t".$threeCCutterSite}++;
	} else {
	    $bowtieBins{$bowtieChr."\t".$threeCCutterSite} = 1;
	}
    }
}

## now print the output in some sort of "sorted" way....

#print "chr\tpos\tnum reads\t6bp cutter\t4bp cutter\tfragment < 500bp\n";
#print "chr\tpos\tnum reads\tmappable\t4base cutter\tshort\t6base cutter inside 4C fragments\n";
print "chr\tpos\treads\tmap\t4bp\tshort\t6bp interferes\tunknown length\n";


foreach my $chr (@chromosomes) {
    print STDERR "Printing out stuff for chromosome ".$chr."\n";
    
    my $threeCposRef = \@{$gThreeCCutterPositions{$chr}};

    my $threeCleftPos = 0; # we start out at the edge of the chromosome
    for(my $threeCindex = 0; 
	$threeCindex <= scalar(@$threeCposRef); $threeCindex++) {

	my $threeCrightPos;
	if ($threeCindex == scalar(@$threeCposRef)) {
	    $threeCrightPos = $gChrSizes{$chr};
	} else {
	    $threeCrightPos = $$threeCposRef[$threeCindex];
	}

	# Make sure that the 3C fragment is mappable.
	#
	#  1: There is NO 4C cutter in the fragment, but the fragment is
	#     short (so that known + unknown < MAX_LENGTH)
	#    
	#    a - If there is no 6base linearize cut site in the unknown
	#           fragment, it is MAPPABLE.
	#    b - If there is a 6base linearize cut site in the unknown
	#           fragment, it is NOT MAPPABLE
	#
	#  2: There is NO 4C cutter in the fragment, and the fragment is
	#       long (so that known + unknown > MAX_LENGTH) - NOT MAPPABLE
	#
	#  3: There is a 4C cutter in the fragment.
	#     
	#     a - If the ligation (of either end of the unknown fragment) 
	#            can only produce something long (known + unknown end >
	#            MAX_LENGTH), it is NOT MAPPABLE
	#
	#     b - If at least one ligation will produce a short 
	#             fragment (known + unknown < MAX_LENGTH), and that
	#             short 
	#
	#         i: If at least one short fragment does not contain 
	#            a linearization cutter, it is MAPPABLE
	#
	#         ii: If all short fragments contain a linearization cutter
	#             it is NOT MAPPABLE.

	my $mappable = FALSE;
	my $fragmentIsSmall = FALSE;
	my $foundFourBase = FALSE;
	my $sixBaseInterferes = FALSE;
	my $sixBaseLinearizes = FALSE;
	my @fourBasePositions;
	
	my $cutterBin = floor($threeCleftPos / BIN_SIZE);
	my $fourCutterBin = $cutterBin;
	my $fourCutterMaxBin = $gMaxFourBaseBins{$chr};
	
	my $doneSearching = FALSE;
	if (!defined($fourCutterBin)) {
	    print STDERR "fourCutterBin was not defined!\n";
	    exit;
	}

	if (!defined($fourCutterMaxBin)) {
	    print STDERR "fourCutterMaxBin was not defined!\n";
	    print STDERR "\tThis error can occur when there are no restriction sites\n";
	    print STDERR "\tdefined in 4C_RE_SITES.txt for chromosome $chr\n";
	    exit;
	}

	while(!$doneSearching && ($fourCutterBin <= $fourCutterMaxBin)) {
	    my $fourBasePosString = $gFourBaseCutterBins{$chr."\t".$fourCutterBin};
	    if (!defined($fourBasePosString)) {
		$fourCutterBin++;
		next;
	    }
	    my @tempFourBasePositions = split(/\t/, $fourBasePosString);
	    for(my $i=0; $i<scalar(@tempFourBasePositions); $i++) {
		if (($tempFourBasePositions[$i] > $threeCleftPos) && 
		    ($tempFourBasePositions[$i] < $threeCrightPos)) {
		    ## We've found a 4 base cutter within the fragment!
		    $foundFourBase = TRUE;
		    $doneSearching = TRUE;
		    push(@fourBasePositions, $tempFourBasePositions[$i]);

#		    print STDERR "fourBasePositions: @fourBasePositions\n";
		} 
		
		if ($tempFourBasePositions[$i] >= $threeCrightPos) {
		    # in this case, there isn't a 4 base cutter within the
		    # current fragment.
		    $doneSearching = TRUE;
		    last;
		}
	    }
	    $fourCutterBin++;
	}

#	if ($foundFourBase) {
#	    print STDERR "threeCleftPos: ".$threeCleftPos."\n";
#	    print STDERR "threeCrightPos: ".$threeCrightPos."\n";
#	    print STDERR "fourBasePositions: @fourBasePositions\n";
#	    print STDERR "\n";
#	}

	my $unknownLength = "NA";

	if (!$foundFourBase) {
	    # now we need to make sure that the 3C fragment is short
	    # (that known + unknown < MAX_LENGTH) and that it does not
	    # contain a linearization enzyme site.
	    $unknownLength = $threeCrightPos - $threeCleftPos;
	    
	    if (($knownLength + $unknownLength) <= $maxFragmentLength) {
		$fragmentIsSmall = TRUE;
		
		# now make sure it does not contain a linearization enzyme.
		#my $foundSixBase = FALSE;
		if (!$noSixBaseCutter) {
		    my $doneSearching = FALSE;
		    my $cutterBin = floor($threeCleftPos / BIN_SIZE);
		    my $sixCutterBin = $cutterBin;
		    my $sixCutterMaxBin = $gMaxSixBaseBins{$chr};
		    while(!$doneSearching && ($sixCutterBin <= $sixCutterMaxBin)) {
			my $sixBasePosString = $gSixBaseCutterBins{$chr."\t".$sixCutterBin};
			if (!defined($sixBasePosString)) {
			    $sixCutterBin++;
			    next;
			}
			my @sixBasePositions = split(/\t/, $sixBasePosString);
			for(my $i=0; $i<scalar(@sixBasePositions); $i++) {
			    if (($sixBasePositions[$i] > $threeCleftPos) && 
				($sixBasePositions[$i] < $threeCrightPos)) {
				## We've found a 6 base cutter within the hind3 fragment!
				$sixBaseInterferes = TRUE;
#			    $foundSixBase = TRUE;
				$doneSearching = TRUE;
				last;
			    } 
			    
			    if ($sixBasePositions[$i] >= $threeCrightPos) {
				# in this case, there isn't a 6 base cutter within the
				# current hind3 fragment.
				$doneSearching = TRUE;
				last;
			    }
			}
			$sixCutterBin++;
		    }
		    #if (!$foundSixBase) {
		    if (!$sixBaseInterferes) {
			$mappable = TRUE;
		    }
		}
	    }
	} else {	    
	    # in this case, we have found 4C cut sites.

	    # so we need to make sure that at least one potential ligation 
	    # with the "known" fragment is <= MAX_LENGTH
	    my $leftFourBase = $fourBasePositions[0];
	    my $rightFourBase = $fourBasePositions[-1];

	    if (!defined($leftFourBase)) {
		print STDERR "leftFourBase not defined!\n";
		print STDERR "@fourBasePositions\n";
		exit;
	    }

	    if (!defined($rightFourBase)) {
		print STDERR "rightFourBase not defined!\n";
		print STDERR "@fourBasePositions\n";
		exit;
	    }


	    # |  | <- left 3C        Both the left and right fragments can
	    # |  |                   ligate with the "known" fragment.  However
	    # |**| <- left four base    we will not know which one before
	    # |**| <- right four base   we sequence.
	    # |  |
	    # |  | <- right 3C


	    my $leftUnknown = $leftFourBase - $threeCleftPos;
	    my $leftFragment = $knownLength + $leftUnknown;

	    my $rightUnknown = $threeCrightPos - $rightFourBase;
	    my $rightFragment = $knownLength + $rightUnknown;

	    if (($leftFragment <= $maxFragmentLength) || 
		($rightFragment <= $maxFragmentLength)) {

		if (($rightFragment <= $maxFragmentLength) && 
		    ($leftFragment <= $maxFragmentLength)) {
		    if ($leftUnknown < $rightUnknown) {
			$unknownLength = $rightUnknown;
		    } else {
			$unknownLength = $leftUnknown;
		    }
		} elsif ($leftFragment <= $maxFragmentLength) {
		    $unknownLength = $leftUnknown;
		} else {
		    $unknownLength = $rightUnknown;
		}
		
		# now we need to make sure that there is not a 6base
		# cut site in both the left and right fragments (if so
		# the site is unmappable.)
		
		my $foundSixBaseLeft = FALSE;
		my $foundSixBaseRight = FALSE;
		if (!$noSixBaseCutter) {
		    my $doneSearching = FALSE;
		    my $cutterBin = floor($threeCleftPos / BIN_SIZE);
		    my $sixCutterBin = $cutterBin;
		    my $sixCutterMaxBin = $gMaxSixBaseBins{$chr};
		    while(!$doneSearching && ($sixCutterBin <= $sixCutterMaxBin)) {
			my $sixBasePosString = $gSixBaseCutterBins{$chr."\t".$sixCutterBin};
			if (!defined($sixBasePosString)) {
			    $sixCutterBin++;
			    next;
			}
			
			my @sixBasePositions = split(/\t/, $sixBasePosString);
			for(my $i=0; $i<scalar(@sixBasePositions); $i++) {
			    
			    if (($sixBasePositions[$i] > $threeCleftPos) && 
				($sixBasePositions[$i] < $leftFourBase)) {
				$foundSixBaseLeft = TRUE;
			    } 
			    
			    if (($sixBasePositions[$i] > $rightFourBase) && 
				($sixBasePositions[$i] < $threeCrightPos)) {
				$foundSixBaseRight = TRUE;
			    } 
			    
			    if ($sixBasePositions[$i] >= $threeCrightPos) {
				# in this case, there isn't a 6 base cutter 
				$doneSearching = TRUE;
				last;
			    }
			    
			    if ($foundSixBaseLeft && $foundSixBaseRight) {
				#$foundSixBase = TRUE;
				$sixBaseInterferes = TRUE;
				$doneSearching = TRUE;
				last;
			    }
			}
			$sixCutterBin++;
		    }
		}
		if ((!$foundSixBaseLeft) || (!$foundSixBaseRight)) {
		    $mappable = TRUE;
		}
#		if ($foundFourBase && ($threeCleftPos == 3015360)) {
#		    print STDERR "mappable: ".$mappable."\n";
#		    print STDERR "foundSixBaseLeft: ".$foundSixBaseLeft."\n";
#		    print STDERR "foundSixBaseRight: ".$foundSixBaseRight."\n";
#		    exit;
#		}
	    }
	}
	
	if (($unknownLength ne "NA") && ($unknownLength <= MIN_READ_LENGTH)
	    && ($unknownLength <= $minUnknownLength)) {
	    $mappable = FALSE;
	}

	if (defined($bowtieBins{$chr."\t".$threeCleftPos})) {
	    print $chr."\t".$threeCleftPos."\t".$bowtieBins{$chr."\t".$threeCleftPos}."\t".$mappable."\t".$foundFourBase."\t".$fragmentIsSmall."\t".$sixBaseInterferes."\t".$unknownLength."\n";
	} else {
	    print $chr."\t".$threeCleftPos."\t0\t".$mappable."\t".$foundFourBase."\t".$fragmentIsSmall."\t".$sixBaseInterferes."\t".$unknownLength."\n";
	}

	$threeCleftPos = $threeCrightPos;
    }
}

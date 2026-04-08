#!/usr/bin/perl

use strict;
use warnings;

my $baitLength = shift(@ARGV);
my $fastqFile = shift(@ARGV);

if ($fastqFile =~ /\.gz\Z/) {
    open(INPUT, "gunzip -c $fastqFile |") 
	|| die("Could not open $fastqFile $!\n");
} else {
    open(INPUT, "<".$fastqFile) || die("Could not open $fastqFile $!\n");
}

while(my $comment1 = <INPUT>) {
    my $seq = <INPUT>;
    my $comment2 = <INPUT>;
    my $phred = <INPUT>;
    
    chomp($seq);
    chomp($phred);

    my $newSeq = substr($seq, $baitLength);
    my $newPhred = substr($phred, $baitLength);

    print $comment1;
    print $newSeq."\n";
    print $comment2;
    print $newPhred."\n";
}

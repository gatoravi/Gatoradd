#! /usr/bin/perl

# Script to convert a family-add taxa file to a random-add taxa file for inputs to gatorADD.
# Usage - perl randomize.pl family-add-taxaFile
# Avinash Ramu, University of Florida

use warnings;
use strict; 

my $ip = $ARGV[0];
my $op = "r".$ip;


open (IN, $ip )  or die"Unable to open input file !";
open (OUT, ">".$op ) or die"Unable to open output file !";

while (<IN>) {
  chomp;
  my @fields = split("\t", $_);
  print OUT"$fields[0]\tRANDOM\n"; #first field is the taxa
}

close(IN);
close(OUT);

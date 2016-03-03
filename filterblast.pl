#!/usr/bin/perl -w
use strict;

open(TEXTFILE, $ARGV[0]);
my$t=$ARGV[1];


while (my$line = <TEXTFILE>) {
    if ($line =~ /^Query/) {
        chomp ($line);
        my @line = split /\t/, $line;
        print $line[0], "\t", $line[1], "\t", $line[2], "\t", $line[3], "\t", $line[4], "\t", $line[7], "\n"; ### print headers
    }
    else {
        chomp ($line);
        my @line = split /\t/, $line;
        if ($line[0] eq $line[1]) {                        ### removes blast hit against itself, don't print
            next;
        }
        elsif (($line[4] < $t) || ($line[7] < $t)){        ### filter for length threshold defined by user, don't print
            next;
        }
        else {
            print $line[0], "\t", $line[1], "\t", $line[2], "\t", $line[3], "\t", $line[4], "\t", $line[7], "\n";  ### if alignment not filtered, print out
        }
    }
}
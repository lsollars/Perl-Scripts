#!/usr/bin/perl -w
use strict;
use List::MoreUtils qw(any);

### Written by Elizabeth Sollars, e.sollars@qmul.ac.uk
### This script takes the output of filterblast.pl to make paralog groups of genes.
### Usage: perl get_groups.pl <inputfile>
### Output: 1st column is identifier of group. Additional columns list IDs of genes within group. To calculate total number of groups run "wc -l <output>"
### Output should be piped into a new file

my $href;
open(TEXTFILE, $ARGV[0]) or die "Can't open file";

my@line2=();
my$line2='';
my $group=1;
my $current_group = '';
my$pair_number=1;
## initiate hashref
$href->{$group}->{$pair_number} = join (':', "000000.1", "00000.2");

### read in input file

OUTERLOOP: while (my$line = <TEXTFILE>) { 
    my$flag=1;
    if ($line =~ /^Query/) {
        next;
    } else {
        chomp $line;
        @line2 = split /\t/, $line;
        my $pair = join (':', $line2[0], $line2[1]);  ### join together identifiers of genes
        foreach my $group (keys %$href) {
            foreach my $pair_number (keys %{$href -> {$group}}) {
                if (($href->{$group}->{$pair_number} =~ $line2[0]) && ($href->{$group}->{$pair_number} =~ $line2[1])){    ##### check that pair isn't already in hash (pair will appear twice)###
                    $flag=2;                                                                                              ### if already in hash, give value flag=2 and go to next line of file, or keep looking through hash
                    next OUTERLOOP;
                } else {
                    next;
                }
            }
        }
        ### once determined pair is not already in hash, go through hash again, through each pair in each group. If one member of pair matches something in the hash, add the pair into the group
        foreach my $group (keys %$href) { 
            foreach my $pair_number (keys %{$href -> {$group}}) {   ### each pair gets an identifier
                if ((($href->{$group}->{$pair_number} =~ $line2[0]) || ($href->{$group}->{$pair_number} =~ $line2[1])) && ($flag ==1)){ ### if one member of pair matches something already in hash, give value flag=3
                    $flag=3;
                    my $count =  scalar (keys %{$href -> {$group}});    ### get number of members in group
                    $current_group = $group;                            ### current group is now this group being added to
                    my $new_count = $count + 1;                         ### add 1 to count
                    $href->{$group}->{$new_count} = $pair;              ### add new pair to group
                    #next OUTERLOOP;
                } elsif ((($href->{$group}->{$pair_number} =~ $line2[0]) || ($href->{$group}->{$pair_number} =~ $line2[1])) && ($flag ==3) && ($current_group ne $group)){ ###if one member of pair matches another group in hash
                    foreach my $pair_number (keys %{$href -> {$group}}) {                           ### for each of the pairs in matched group
                        my $count =  scalar (keys %{$href -> {$current_group}});                    ### get count of genes in current group (from the last if statement)
                        my $new_count = $count + 1;
                        $href->{$current_group}->{$new_count} = $href->{$group}->{$pair_number};    ### add pair to current group
                        $href->{$group}->{$pair_number} = 0;                                        ### give pair value 0 (done with all members of matched group as they now belong in the current group)
                    }
                } else {
                    next;
                }
            }
        }
        if ($flag == 1) {                                  ### if pair hasn't matched anything else in the hash, make a new group for this pair, with new group identifier
            my $group_count = scalar (keys %$href);
            my $new_group_count = $group_count + 1;
            $href->{$new_group_count}->{'1'} = $pair;      ### pair will be the first pair of this new group
        } 
    }
}

## print out set of all group members ##

delete $href->{1}; #delete the first key value of hashref, as it was just used to initiate
foreach my $group (sort keys %$href) { ###sort groups numerically on group identifier
    
    if ($href->{$group}->{1} eq 0){    ### if the pair has value 0, ignore it
        next;
    }
    elsif ((scalar (keys %{$href-> {$group}}) == 1) && ($href->{$group}->{$pair_number} ne 0)){  ### if group contains just one pair, split it on the colon, then print group identifier and both gene identifiers
        my @pair = split /:/, $href->{$group}->{'1'};
        print $group, " ", $pair[0], " ", $pair[1], "\n";
    } else {                                                                  ### if group contains more than one pair, make array of gene identifiers and print them
        my @members=();
        foreach my $pair_number (sort keys %{$href -> {$group}}) {
            my @pair = split /:/, $href->{$group}->{$pair_number};
            if (scalar (@members) == 0){
                push @members, $pair[0], $pair[1];
            }
            else {
            if (any {$pair[0] eq $_ } @members) {
            } else {
                push @members, $pair[0];
            }
            if (any {$pair[1] eq $_} @members) {
            } else {
                push @members, $pair[1];
            }   
        }
    }
    print $group, " ";
        foreach my $string (@members){
            print $string, " ";
        }
        print "\n";    
    }
}


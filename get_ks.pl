#!/usr/bin/perl -w
use strict;

### Written by Elizabeth Sollars, e.sollars@qmul.ac.uk
### This script takes the output of get_groups.pl to calculate a set of ks values for each pairwise combination of genes within a paralog group
### Usage: perl get_seq.pl <inputfile> <cds.fa> <peptide.fa>
### cds and peptide fasta files must have same gene identifiers / headers as those in the output of get_groups.pl, otherwise the script will not be able to find the gene
### output is automatically put in a file named "ks_values.txt", other small intermediate files are created and deleted along the way but will not remain once the script is finished
### Output: 1st column is group identifer, next two columns are gene identifiers, fourth column is the ks value for the two genes
### Requirements: clustalw2, Pal2Nal.pl (http://www.bork.embl.de/pal2nal/), PAML yn00 (http://abacus.gene.ucl.ac.uk/software/paml.html)
### NOTE: You will probably need to change the path to these programs in the 'system' commands on lines 78, 81, 89 and 129. 
### Some of the ks values in groups with more than two genes will be redundant (see Maere et al 2005), and should be corrected; a group with n genes should have n-1 duplications events, not n! events.

open(GROUPFILE, $ARGV[0]) or die "Can't open group file";
open(SEQFILE, $ARGV[1]) or die "Can't open seq file";
open(PEPFILE, $ARGV[2]) or die "Can't open pep file";
open(OUTFILE, ">>ks_values.txt") or die "Can't open output file";

my%seq_hash;
my%pep_hash;
my $start = 0;
my $seq_id = '';

## make hash of cds sequences, header is key and sequence is value ##

while (my $line = <SEQFILE>) {
    if ($line =~ /^>.*/) {
        chomp $line;
        $seq_id = substr $line, 1;
        $seq_id =~ s/\s//g;
        $start = 1;
    } elsif (($line !~ /^>/) && ($start == 1)){
        chomp $line;
        $seq_hash{$seq_id} = $line;
        $start = 2;
    } elsif (($line !~ /^>/) && ($start == 2)){
        chomp $line;
        $seq_hash{$seq_id} = $seq_hash{$seq_id}.$line;
    }
}

## make hash of peptide sequences, header is key and sequence is value ##

while (my $line = <PEPFILE>) {
    if ($line =~ /^>.*/) {
        chomp $line;
        $seq_id = substr $line, 1;
        $seq_id =~ s/\s//g;
        $start = 1;
    } elsif (($line !~ /^>/) && ($start == 1)){
        chomp $line;
        $pep_hash{$seq_id} = $line;
        $start = 2;
    } elsif (($line !~ /^>/) && ($start == 2)){
        chomp $line;
        $pep_hash{$seq_id} = $pep_hash{$seq_id}.$line;
    }
}

## get groups ##

while (my$line = <GROUPFILE>){
    chomp $line;
    my @group = split /\s/, $line;
    my $group_id = shift @group;                     ### get the group identifier 
    my $count = scalar(@group) - 1;                  ### get number of genes in group
    for (my$i=1; $i<=$count; $i++){
        my $first_gene = shift @group;               ### get first gene and remove it from array, calculate KS against all other genes in group, repeat with all other genes until last pair
        foreach my $string (@group){
            ## Get sequences from cds and peptide files ##
            open(TWOPEPFILE, ">two_peps.fa");
            open(TWOSEQFILE, ">two_seqs.fa");
            print TWOSEQFILE "\>", $first_gene, "\n", $seq_hash{$first_gene}, "\n\>", $string, "\n", $seq_hash{$string}, "\n";
            print "Seq file created\n";
            print TWOPEPFILE "\>", $first_gene, "\n", $pep_hash{$first_gene}, "\n\>", $string, "\n", $pep_hash{$string}, "\n";
            print "Pep file created\n";
            ## Get Clustal alignment ## 
            system ('clustalw2 -infile=two_peps.fa -align -type=protein -outfile=pep.aln -output=clustal -pwmatrix=BLOSUM');
            print "Clustal alignment created\n";
            ## Align using Pal2NAL ##
            system ('perl pal2nal.v14/pal2nal.pl pep.aln two_seqs.fa -output clustal -nogap > codon.aln');
            print "Pal2Nal alignment created\n";
            system ('rm pep.aln');
            close TWOPEPFILE;
            close TWOSEQFILE;
            system ('rm two_peps.fa');
            system ('rm two_seqs.fa');
            ## Convert to phylip format ##
            system ('clustalw2 -infile=codon.aln -convert -output=phylip -outfile=alignment.phy');
            print "Alignment converted to phylip\n";
            system ('rm codon.aln');
            ## Add I to phylip file ##
            my$filename="alignment.phy";
            if (-e $filename) {
            open(PHYLIPFILE, "alignment.phy") or die "Can't open phylip file";
            open(NEWPHYLIPFILE, ">>alignment2.phy") or die "Can't open phylip2 file";  ###create new phylip file with correct format for PAML
            my $p = 1;
            while (my $input = <PHYLIPFILE>) {
                if ($p == 1) {
                    chomp $input;
                    print NEWPHYLIPFILE $input, "\t I\n";
                    print $input, "\t I\n";
                    $p = 2;
                } elsif ($p == 2){
                    chomp $input;
                    my @sl = split /\s+/, $input;
                    shift @sl;
                    my $nl = join " ", @sl;
                    print NEWPHYLIPFILE "Seq1    ", $nl, "\n";
                    print "Seq1    ", $nl, "\n";
                    $p = 3;
                } elsif ($p == 3){
                    chomp $input;
                    my @sl = split /\s+/, $input;
                    shift @sl;
                    my $nl = join " ", @sl;
                    print NEWPHYLIPFILE "Seq2    ", $nl, "\n";
                    print "Seq2    ", $nl, "\n";
                    $p = 4;
                } else {
                    chomp $input;
                    print NEWPHYLIPFILE $input, "\n";
                    print $input, "\n";
                } 
            }
            close PHYLIPFILE;
            #system('rm alignment.phy');
            ## run yn00 ##
            system ('paml4.8/bin/yn00 paml4.8/yn00.ctl');
            print "yn00 completed\n";
            ## Read yn00 output file and get kS value ##
            open(YN00FILE, "2YN.dS") or die "Can't open yn00 output file";
            my $y=1;
            while (my $input = <YN00FILE>) {
                if ($y==3) {
                    chomp $input;
                    print "Input here: ", $input, "\n";
                    my @input = split /\s+/, $input;
                    print $input[0], "\n", $input[1], "\n";
                    print OUTFILE $group_id, "\t", $first_gene, "\t", $string, "\t", $input[1], "\n";  ###print results to the output file
                    $y++;
                } else {
                    $y++;
                    next;
                }
            }
            close NEWPHYLIPFILE;
            close YN00FILE;
            system ('rm 2YN.dS');
            system ('rm alignment2.phy');
            }
            system ('rm alignment.phy');
        }
    }
}


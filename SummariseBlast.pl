#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::SearchIO;
#use Test::More;

my ($help, $man)='0';
my $output = 'first';
my $outfile = 'AlignmentOutput.txt';
GetOptions ("outformat:s" => \$output, "outfile:s" => \$outfile, "help" => \$help, "man" => \$man);

pod2usage(-exitval => 0, -verbose => 1) if $help;
pod2usage(-exitval => 0, -verbose => 2) if $man;


if (!(defined $ARGV[0])) {
    print "\n***Please specify an input file***\n";
    pod2usage(-exitval => 0, -verbose => 1);
}

open(OUTFILE, ">>$outfile");

my$y=0;
my$p=0;
my($count, $previous_length, $previous_added)=0;
my$current='';

my $in = new Bio::SearchIO(-format => 'blast', 
                           -file   => <$ARGV[0]>);

## Write Header ##
print OUTFILE "Query Name\tHit Name\tHit Length\tAlignment Length\t% Hit Length\tNumber Conserved\t% Hit Conserved\t% Query Length\n";

## Read in each query result, initialise 'count' score ##
OUTERLOOP: while(my $result = $in->next_result) {
    $count=0;
##If query has no hits, skip it ##
    if ($result->num_hits == 0) {
        next;
    }
    
my@matrix=();
    
## Read in each hit, add 1 to count (count marks number of hit) ##
  while(my $hit = $result->next_hit) {
     $count++;
## When 'first' output mode is selected, skip anything that isn't the first hit ##
    if (($output eq 'first') && ($count > 1)) {
        next OUTERLOOP;
    }
    
## If hit has only 1 HSP, calculate percentage coverage of query and hit sequences, and conservation percentage of hit. Print results. ##
    if ($hit->num_hsps == 1){
      while(my $hsp = $hit->next_hsp) {
        my$hit_length_percent = sprintf '%.5f', (($hsp->length('hit') / $hit->length) * 100);
        my$hit_cons_percent = sprintf '%.5f',(($hsp->num_conserved / $hit->length) * 100);
        my$query_length_percent = sprintf '%.5f',(($hsp->length('query') / $result->query_length) * 100);
        
##if 'first' or 'all' output selected, print the results straight away. ##
        if (($output eq 'first') || ($output eq 'all')){
            #cmp_ok($hit_length_percent, '<=', 100);
            #cmp_ok($hit_cons_percent, '<=', 100);
            #cmp_ok($query_length_percent, '<=', 100);
            #cmp_ok($hit_cons_percent, '<=', $hit_length_percent);
            print OUTFILE $result->query_name, "\t",
            $hit->name, "\t",
            $hit->length, "\t",
            $hsp->length('hit'), "\t",
            $hit_length_percent, "\t",
            $hsp->num_conserved, "\t",
            $hit_cons_percent, "\t",
            $query_length_percent, "\n";
            }
        elsif ($output eq 'best'){
            
##for 'best' output, put all results for the query into 2-d array. $count distuingishes each hit from next##
            $matrix[$count][0] = $hit->name;
            $matrix[$count][1] = $hit->length;
            $matrix[$count][2] = $hsp->length('hit');
            $matrix[$count][3] = $hit_length_percent;
            $matrix[$count][4] = $hsp->num_conserved;
            $matrix[$count][5] = $hit_cons_percent;
            $matrix[$count][6] = $query_length_percent;
            
##Put lowest value of length percent (either of query sequence, or more rare, of protein sequence) in another column which will be sorted later##
            if ($hit_length_percent < $query_length_percent) {
                $matrix[$count][7] = $hit_length_percent;
            } else {
                $matrix[$count][7] = $query_length_percent;
            } 
        }
      }
    }
    
##If hit has more than 1 HSP, invoke a loop to test whether each position in query / hit sequence is covered by at least one HSP. ##
##This takes a bit longer but allows for overlaps between the HSPs ##
    elsif ($hit->num_hsps > 1) {
        my%hash;
        my@conserved_points=(); ##array for conserved positions ##
        my$conservation=0;
        my@hsp_query_points=(); ##array for query positions covered in HSP ##
        my@hsp_hit_points=(); ## array for hit positions covered in HSP ##
        my$query_length=0;
        my$hit_length=0;
        
##loop through all HSPs of hit, noting conserved positions, and start and end positions of hit and query sequences##
        while(my $hsp = $hit->next_hsp) { 
            push @conserved_points, $hsp->seq_inds('hit','conserved');
            push @hsp_query_points, $hsp->range('query');
            push @hsp_hit_points, $hsp-> range('hit');
        }
        
##Enter all conserved_points into hash, then we can test later which positions are 'defined' and calculate conservation score ##
        foreach my$value (@conserved_points){
            if (exists($hash{$value})) { ##if position is already in hash don't bother defining it again (will occur if HSPs overlap) ##
                next;
            } else {
            $hash{$value}='yes';
            }
        }
        
## loop through every position in hit, see if it is defined in conservation hash. If so, add 1 to the conservation score ##
        for (my$i=1; $i<= $hit->length; $i++){
            if (exists ($hash{$i})) {
                $conservation++;
            }
        }
        
## Loop through each position in query, and see if it is covered between a 'start' and 'end' position in hsp_query_points array. If so, add 1 to query length ##
        LOOP1: for ($y=1; $y<= ($result->query_length); $y++){
            for ($p=0; $p<= (scalar(@hsp_query_points)-1); $p=$p+2){
                if (($y >= $hsp_query_points[($p)]) && ($y <= $hsp_query_points[($p+1)])){
                    $query_length++;
                    next LOOP1;
                } else {
                    next;
                }
            }
        }
        
## Same as LOOP1, but with hit positions ##
        LOOP2: for ($y=1; $y<= ($hit->length); $y++){
            for ($p=0; $p<= (scalar(@hsp_hit_points)-1); $p=$p+2){
                if (($y >= $hsp_hit_points[($p)]) && ($y <= $hsp_hit_points[($p+1)])){
                    $hit_length++;
                    next LOOP2;
                } else {
                    next;
                }
            }
        }
        
## Calculated percentages from scores ##
        my$hit_cons_percent = sprintf '%.5f',(($conservation / $hit->length) * 100);
        my$query_length_percent = sprintf '%.5f',(($query_length / $result->query_length) * 100);
        my$hit_length_percent = sprintf '%.5f',(($hit_length / $hit->length) * 100);
        
##Print results for 'all' and 'first' output format##
        if (($output eq 'first') || ($output eq 'all')){
            #cmp_ok($hit_length_percent, '<=', 100);
            #cmp_ok($hit_cons_percent, '<=', 100);
            #cmp_ok($query_length_percent, '<=', 100);
            #cmp_ok($hit_cons_percent, '<=', $hit_length_percent);
            print OUTFILE $result->query_name, "\t",
            $hit->name, "\t",
            $hit->length, "\t",
            $hit_length, "\t",
            $hit_length_percent, "\t",
            $conservation, "\t",
            $hit_cons_percent, "\t",
            $query_length_percent, "\n";
        }
        elsif ($output eq 'best'){
            
##for 'best' output, put all results for the query into 2-d matrix. $count distuingishes each hit from next##
            $matrix[$count][0] = $hit->name;
            $matrix[$count][1] = $hit->length;
            $matrix[$count][2] = $hit_length;
            $matrix[$count][3] = $hit_length_percent;
            $matrix[$count][4] = $conservation;
            $matrix[$count][5] = $hit_cons_percent;
            $matrix[$count][6] = $query_length_percent;
            
##Put lowest value of length percent (either of query sequence, or more rare, of protein sequence) in another column which will be sorted later##
            if ($hit_length_percent < $query_length_percent) {
                $matrix[$count][7] = $hit_length_percent;
            } else {
                $matrix[$count][7] = $query_length_percent;
            } 
        } 
    }
}
  
    if ($output eq 'best') {

##Now calculate the 'best hit'##
    $previous_length = 0;
    $previous_added = 0;
##loop through last column of matrix to find the longest length##

    for (my$i=1; $i<=$count; $i++){
        if ($matrix[$i][7] > $previous_length) {
            $previous_length = $matrix[$i][7];
        } else {
            next;
        } 
    }
    
##find all lengths that are within 5% of the longest, then find the hit with largest conservation+length out of those results##

    for (my$i=1; $i<=$count; $i++){
        if ((($previous_length - $matrix[$i][7]) < 5) && (($matrix[$i][5] + $matrix[$i][7]) > $previous_added)) {
            $previous_added = ($matrix[$i][5] + $matrix[$i][7]);
            $current = $i;
        } else {
            next;
        }
    }

##print the best result##
    #cmp_ok($matrix[$count][3], '<=', 100);
    #cmp_ok($matrix[$count][5], '<=', 100);
    #cmp_ok($matrix[$count][6], '<=', 100);
    #cmp_ok($matrix[$count][5], '<=', $matrix[$count][3]);
    print OUTFILE $result -> query_name, "\t",
    $matrix[$current][0], "\t",
    $matrix[$current][1], "\t",
    $matrix[$current][2], "\t",
    $matrix[$current][3], "\t",
    $matrix[$current][4], "\t",
    $matrix[$current][5], "\t",
    $matrix[$current][6], "\n",
  }
}

#done_testing();

###Documentation###

__END__

=head1 NAME
 
SummariseBlast.pl

=head1 SYNOPSIS

perl SummariseBlast.pl <inputfile> [options]

=head1 OPTIONS

--outformat=first   (default and quickest) Prints only the first hit's result. 

           =all     Prints all BLAST alignment results.

           =best    Prints only the 'best' hit's result.

--outfile=<outputfile>  (By default, output is written to a file named "AlignmentOutput.txt", but an alternative can be specified). 

-help prints usage and options

-man prints full manual

=head1 DESCRIPTION

This Perl script takes input from a pairwise BLAST search result. It prints some summary information from the alignments, so that the user can get a quick overview of the quality and length of the alignment.
It is designed to be used on a BLAST search against a database of genes, transcripts, or proteins, and not on whole genomes. For example, it could be used to accompany a homology search, in order to
see how well genes from a transcript assembly match with mRNA from a related organism, or in functional annotation to see the matches for a set of unannotated proteins.
It takes into account multiple HSPs per hit, if present, and whether they overlap in either the query or hit sequence. 

The output is the following summary information:

Query Name

Hit Name (or Accession number)

Hit Length = total length of the hit sequence

Alignment Length = length of the alignment, including all HSPs.

% Hit Length = Percentage of the hit in the alignment. (Alignment length / Hit length) * 100.

Number conserved = number of conserved amino acids.

% Hit Conserved = percentage of amino acids in hit sequence that are conserved in the alignment. (Number conserved / Hit Length) * 100.

% Query Length = percentage of query sequence participating in alignment.


About the three --outformat options:

"all" prints the above information for every hit result. This means that if there are 500 hits per query sequence, then this script will write 500 rows of information for each query sequence.
If only one result per query is needed, consider the "best" and "first" options. "First" (the default and quickest option) only considers the very first hit from the BLAST output, which has the lowest e-value. It can be used to get a quick first overview of the Blast results. 
Whereas, "best" does some sorting to find the actual best hit (which is not necessarily the first hit, although often will be). It first finds the hit with the highest % Query Length score (or % hit score if this is lower),
and then gets all other hits that are within 5% of the top scorer (e.g. if top score is 45%, all results with 40% and higher will be included). These "top" hits, are then sorted based on the sum of: % conservation + % Query Length.
The highest of these is chosen as the 'best' hit, and its alignment information is printed. The rationale behind this option, is that a low e-value can come from a very short but specific alignment
(i.e. hitting only one domain of a protein) and may end up as the first hit. Whereas a longer but less specific alignment (i.e. hitting the whole length of a protein) may have a higher e-value and therefore be lower down on the list,
however many would actually consider this the better alignment, because a larger amount of the hit sequence is covered. 

The output is a tab-delimited file, which can be imported easily into other programs such as Excel and R. If a query sequence had no hits in the BLAST result, it is ignored and will be missing from the output file.

This script was written by Elizabeth Sollars, Queen Mary University of London. Any bugs or friendly suggestions, please email e.sollars@qmul.ac.uk

=cut

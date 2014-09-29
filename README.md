Perl-Scripts
============

SummariseBlast.pl
-----------------

Usage: perl SummariseBlast.pl <inputfile> [options]

This Perl script takes input from a pairwise BLAST search result. It prints some summary information from the alignments, so that the user can get a quick overview of the quality and length of the alignment.
It is designed to be used on a BLAST search against a database of genes, transcripts, or proteins; not on whole genomes. For example, it could be used to accompany a homology search, in order to see how well genes from a transcript assembly match with mRNA from a related organism, or in functional annotation to see the matches for a set of unannotated proteins. 
It takes into account multiple HSPs per hit, if present, and whether they overlap in either the query or hit sequence. 

Options:

           --outformat =first   (default and quickest) Prints only the first hit's result. 

                       =all     Prints all BLAST alignment results.

                       =best    Prints only the 'best' hit's result. (see further down for explanation).

           --outfile= <outputfile>  (By default, output is written to a file named "AlignmentOutput.txt", but an alternative can be                       specified). 
           -help prints usage and options

           -man prints full manual

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

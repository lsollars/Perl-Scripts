Perl-Scripts
============

Ks Pipeline
------------

These are a few scripts which together, form a pipeline for investigating historic Whole Genome Duplications (WGD) in a genome. They group genes with similar sequences into paralog groups, and calculate the synonymous substitutions (KS) between genes within each group.  The assumption behind the KS method is that these paralogs arose through duplications; either local or genome-wide. A higher KS value (more synonymous substitutions), suggests that the duplication event that created the two copies occurred longer ago and the two seqeunces have had more time to diverge. By calculating and plotting all KS scores, one should observe a peak in values where a WGD is likely to have occured, as many paralogs point to roughly the same KS value. The absence of a peak suggests that no WGD has occurred (or that it occured so long ago that the KS method cannot detect it), and the current paralogs arose through local duplications. The scripts work like this:

1). SummariseBlast.pl takes a pairwise BLAST result (i.e. an all-against-all BLAST of CDS or protein sequences from one genome) and prints some summary statistics such as alignment length.
2). FilterBlast.pl takes the summary statistics and filters them to remove: A) the hit where the genes matches itself and B) hits that are below a certain length threshold (default is 50% of the longest gene, but this can be changed).
3). get_groups.pl takes the filtered summary information and creates groups of genes that would be considered paralogs.
4). get_ks.pl calculates KS values between every pairwise combination of genes within each paralog group. To do this it calls the following programs: clustal, Pal2Nal.pl (http://www.bork.embl.de/pal2nal/), and PAML yn00 (http://abacus.gene.ucl.ac.uk/software/paml.html). 
5). reduce_redundant.py is a script written by my colleague Endymion Cooper at QMUL. It takes the complete set of KS values belonging to each paralog group and corrects them for redundant values. This step is important for the following reason: A group of n genes should only have n-1 single duplication events between them. Whereas by computing the KS values between all pairwise combinations of genes, we in fact obtain n(n-1)/2 values (Maere et al (2005) 10.1073/pnas.0501102102). This script corrects for these redundant values using a similar method as in Maere et al (2005). 
6). The Ks values can then be plotted as a histogram, e.g. in R. 

More information on each script is contained within the scripts themselves, so please have a read. As SummariseBlast.pl has a few options, they are explained below in more detail. 

SummariseBlast.pl
-----------------

Usage: perl SummariseBlast.pl inputfile [options]

This Perl script takes input from a pairwise BLAST search result. It prints some summary information from the alignments, so that the user can get a quick overview of the quality and length of the alignment.
It is designed to be used on a BLAST search against a database of genes, transcripts, or proteins; not on whole genomes. For example, it could be used to accompany a homology search, in order to see how well genes from a transcript assembly match with mRNA from a related organism, or in functional annotation to see the matches for a set of unannotated proteins. 
It takes into account multiple HSPs per hit, if present, and whether they overlap in either the query or hit sequence. 
The Bio::SearchIO module from BioPerl is used in this script (http://www.bioperl.org/wiki/Module:Bio::SearchIO), so please ensure you have BioPerl modules installed and accessible before using the script.

Options:

           --outformat =first   (default and quickest) Prints only the first hit's result. 

                       =all     Prints all BLAST alignment results.

                       =best    Prints only the 'best' hit's result. (see further down for explanation).

           --outfile= <outputfile>  (By default, output is written to a file named "AlignmentOutput.txt", 
                      but an alternative can be specified). 
           
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

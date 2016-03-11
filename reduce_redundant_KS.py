##########################################################################################
#                                                                                        #
#   Created by Endymion D. Cooper (endymion.dante.cooper@gmail.com)                      #
#   This is version 1.1 created on 2nd Sept 2015                                         #
#                                                                                        #
#   Modified from version 1.0 (25th August 2015) to incorporate a second method for      #
#   calculating the per duplication KS value (see details of modes below).               #
#                                                                                        #
#   reduce_redundant_KS.py takes a set of pairwise KS values for clustered sequences     #
#   and uses a UPGMA like clustering method (following Maere et al. PNAS 2005            #
#   DOI:10.1073/pnas.0501102102) to approximate gene duplications and reduce redundant   # 
#   KS values.                                                                           #
#                                                                                        #
#   Maere et al. (2005) Supporting methods:                                              #
#                                                                                        #
#       " Correction for redundant KS values.                                            #
#         A gene family of n members originates from n-1 retained single gene            #
#         duplications, whereas the number of possible pairwise comparisons (KS          #
#         measurements) within a gene family is n(n-1)/2. To correct for the             #
#         redundancy of KS values when building the age distribution for duplicated      #
#         genes, we constructed tentative phylogenetic trees for each gene family        #
#         with an average linkage clustering algorithm using KS as a distance            #
#         measure, similar to the approach adopted by Blanc and Wolfe (1). Starting      #
#         from each gene as a separate cluster, the two clusters with the lowest         #
#         mean inter-cluster KS value (i.e. the mean of all observed KS values           #
#         (edges) between two clusters) were iteratively merged. The splits in           #
#         the resulting average linkage tree represent the n-1 retained duplication      #
#         events. For each split, the m KS measurements between the two merged           #
#         gene clusters were added to the KS distribution with a weight 1/m. In          #
#         other words, all KS estimates for a particular duplication event were          #
#         added to the KS distribution, while the total weight of a single               #
#         duplication event sums up to one.  "                                           #
#                                                                                        #
#                                                                                        #
#   Dependancies: written for python 2.7, builtins and standard libraries.               #
#                                                                                        #
#   Usage: python reduce_redundant_KS.py <input_file> <output_prefix> <M1/M2>            #
#                                                                                        #
#   The input file lines should look like this (tab separated):                          #
#      <cluster_ID>   <sequence_ID>   <sequence_ID>   <KS_value>                         #
#                                                                                        #
#   Two outputs are generated:                                                           #
#      1. <output_prefix>_KS_by_cluster.txt                                              #
#      2. <output_prefix>_KS.txt                                                         #
#                                                                                        #
#   Two methods of calculating mean KS between sub-clusters are available.               #
#   Specify M1 (mode 1) or M2 (mode 2).                                                  #
#   The differences are explained below with an example.                                 #
#                                                                                        #
#   Example data:                                                                        #
#   <Cluster ID>  <Gene1 ID>  <Gene2 ID>  <KS>                                           #
#     cluster1        G1          G2      0.003                                          #
#     cluster1        G1          G3      0.326                                          #
#     cluster1        G1          G4      0.563                                          #
#     cluster1        G2          G3      0.245                                          #
#     cluster1        G2          G4      0.637                                          #
#     cluster1        G3          G4      0.476                                          #
#                                                                                        #
#   Denote:  - the KS score between two sequences/subclusters as e.g. G1:G2.             #
#            - a subcluster of sequences as e.g. G1/2                                    #
#                                                                                        #
#   First iteration:                                                                     #
#   The smallest KS value is between G1 & G2 therefore join as G1/2. Tree:               #
#        G1 __                                                                           #
#        G2 __|  <-- Duplication event 1. KS score G1:G2 = 0.003                         #
#        G3                                                                              #
#        G4                                                                              #
#   And recalculate the KS scores by averaging (note that mode 1 = mode 2 here):         #
#    -------------------------------------------------------------------------           #
#   | New Pairs |            Mode 1            |            Mode 2            |          #
#   |-----------|------------------------------|------------------------------|          #
#   | G1/2  G3  | ((G1:G3)+(G2:G3))/2 = 0.2855 | ((G1:G3)+(G2:G3))/2 = 0.2855 |          #
#   | G1/2  G4  | ((G1:G4)+(G2:G4))/2 = 0.6    | ((G1:G4)+(G2:G4))/2 = 0.6    |          #
#   | G3    G4  | (raw value)         = 0.476  | (raw value)         = 0.476  |          #
#    -------------------------------------------------------------------------           #
#   Second iteration:                                                                    #
#   The smallest KS value is between G1/2 and G3 therefore join as G1/2/3. Tree:         #
#        G1 __                                                                           #
#        G2 __|__                                                                        #
#        G3 _____|  <-- Duplication event 2. KS score G1/2:G3 = 0.2855                   #
#        G4                                                                              #
#   And recalculate the KS scores by averaging (mode 1 =/= mode 2):                      #
#    ----------------------------------------------------------------------------------  #
#   |  New Pairs |                Mode 1               |            Mode 2             | #
#   |------------|-------------------------------------|-------------------------------| #
#   | G1/2/3  G4 | ((G1:G4)+(G2:G4)+(G3:G4))/3 = 0.559 | ((G1/2:G4)+(G3:G4))/2 = 0.538 | #
#    ----------------------------------------------------------------------------------  #
#                                                                                        #
#   Third iteration 3:                                                                   #
#   The smallest KS value is between G1/2/3 and G4 therefore join as G1/2/3/4. Tree:     #
#        G1 __                                                                           #
#        G2 __|__                                                                        #
#        G3 _____|__                                                                     #
#        G4 ________|   <-- Duplication event 3. KS score G1/2/3:G4 = 0.559 (mode 1)     #
#                                                                   = 0.538 (mode 2)     #
#                                                                                        #
##########################################################################################

import sys

# Parse command line arguments.
input_file = sys.argv[1]
output_file = sys.argv[2]
mode = sys.argv[3]

# First open the file and read the lines, store the data in a dictionary indexed by 
# cluster name, with values as a list of lines for each cluster
file = open(input_file, 'r')
lines = file.readlines()
file.close()
data_dict = {}
for line in lines:
    if line.startswith('\n'):
        continue
    else:
        L1 = line.rstrip('\n')
        L2 = L1.split('\t')
        if L2[0] not in data_dict:
            data_dict[L2[0]]=[L2]
        else:
            data_dict[L2[0]].append(L2)

# Define the primary function to do the clustering and extract the 
# non-redundant KS values for each blast cluster.
def remove_redundant_KS(lines):
    # list to hold current clusters
    clusters = []
    # list to hold sequence names
    sequences = []
    # dictionary to hold cluster definitions
    cluster_defs = {}
    # dictionary to hold all KS scores
    all_KS = {}
    # dictionary to hold KS scores for current clusters
    current_KS = {}
    # list for final KS scores
    KS = []
    # a counter to append to subcluster IDs
    counter = 0

    # populate the lists and dictionaries with initial values.
    for line in lines:
        L2 = line
        L3 = frozenset([L2[1],L2[2]])
        # populate the KS value dictionary
        all_KS[L3] = L2[3]
        # populate the cluster and sequence lists
        if L2[1] not in clusters:
            clusters.append(L2[1])
            sequences.append(L2[1])
        if L2[2] not in clusters:
            clusters.append(L2[2])
            sequences.append(L2[2])
    # populate the cluster definitions
    for i in xrange(0,len(clusters)):
        cluster_defs[clusters[i]]=[clusters[i]]

    # define a function to get the KS scores for the current clusters
    # KS scores are held in a dictionary and sent to get_min_KS to extract the smallest KS value
    def get_current_KS():
        dict = {}
        temp = []
        for x in xrange(0,len(clusters)):
            for y in xrange(0,len(clusters)):
                if x != y:
                    # This is the simplest case.
                    if frozenset([clusters[x],clusters[y]]) in all_KS:
                        if frozenset([clusters[x],clusters[y]]) not in dict:
                            dict[frozenset([clusters[x],clusters[y]])] = all_KS[frozenset([clusters[x],clusters[y]])]
                        else:
                            continue
                    # Here for merged clusters. Put pairs in temporary list. 
                    else:
                        if frozenset([clusters[x],clusters[y]]) not in temp:
                            temp.append(frozenset([clusters[x],clusters[y]]))
                        else:
                            continue
                else:
                    continue
        # For merged clusters, calculate the mean of all KS between clusters.
        for key in temp:
            pairs = []
            values = []
            # Values in first cluster of pair
            for x in xrange(0,len(cluster_defs[list(key)[0]])):
                # Values in second cluster of pair
                for y in xrange(0,len(cluster_defs[list(key)[1]])):
                    combi = frozenset([list(cluster_defs[list(key)[0]])[x],list(cluster_defs[list(key)[1]])[y]])
                    if combi not in pairs:
                        pairs.append(combi) 
                        values.append(float(all_KS[combi]))
            dict[key]=str(sum(values)/len(values))
            # If mode is M2 we update all_KS to treat subclusters as if they were terminal sequences for next iteration.
            if mode == 'M2':
            	all_KS[key]=str(sum(values)/len(values))
        # If mode is M2 we update cluster definitions to treat subclusters as if they were terminal sequences for next iteration.
        if mode == 'M2':
        	for i in xrange(0,len(clusters)):
        		cluster_defs[clusters[i]]=[clusters[i]]
        return dict

    # get the min KS in current clustering
    def get_min_KS(ks_dict):
        min_KS = {}
        pair = min(ks_dict, key = ks_dict.get)
        min_KS[pair] = ks_dict[pair]
        # Append the KS score to the final list of KS scores.
        KS.append(ks_dict[pair])
        return min_KS

    # update the current clusters
    # having identified the smallest KS, those sequences/sub-clusters are merged and the cluster list modified.
    def update_clusters(cluster_list,min_KS):
        temp = []
        new_clusters = []
        for x in xrange(0,len(clusters)):
            for key in min_KS:
                # If KS value for a pair of sequences.
                if list(key)[0] in sequences and list(key)[1] in sequences:
                    cluster_defs["sub_clust_"+str(counter)]=list(key)
                    if clusters[x] in list(key):
                        temp.append(clusters[x])
                    else:
                        continue
                # If KS value for a sequence and a cluster of sequences.
                elif list(key)[0] in sequences or list(key)[1] in sequences:
                    seq_list = []
                    if list(key)[1] in sequences:
                        seq_list.append(str(list(key)[1]))
                        if len(cluster_defs[list(key)[0]]) > 1:
                            for y in xrange(0,len(cluster_defs[list(key)[0]])):
                                if cluster_defs[list(key)[0]][y] not in seq_list:
                                    seq_list.append(cluster_defs[list(key)[0]][y])
                    if list(key)[0] in sequences:
                        seq_list.append(str(list(key)[0]))
                        if len(cluster_defs[list(key)[1]]) > 1:
                            for y in xrange(0,len(cluster_defs[list(key)[1]])):
                                if cluster_defs[list(key)[1]][y] not in seq_list:
                                    seq_list.append(cluster_defs[list(key)[1]][y])
                    cluster_defs["sub_clust_"+str(counter)]=seq_list
                    if clusters[x] in list(key):
                        temp.append(clusters[x])
                    else:
                        continue
                # If KS value for two clusters of sequences.
                else:
                    seq_list = []
                    if len(cluster_defs[list(key)[0]]) > 1:
                        for y in xrange(0,len(cluster_defs[list(key)[0]])):
                            if cluster_defs[list(key)[0]][y] not in seq_list:
                                seq_list.append(cluster_defs[list(key)[0]][y])
                    if len(cluster_defs[list(key)[1]]) > 1:
                        for y in xrange(0,len(cluster_defs[list(key)[1]])):
                            if cluster_defs[list(key)[1]][y] not in seq_list:
                                seq_list.append(cluster_defs[list(key)[1]][y])
                    cluster_defs["sub_clust_"+str(counter)]=seq_list
                    if clusters[x] in list(key):
                        temp.append(clusters[x])
                    else:
                        continue
        for y in xrange(0,len(clusters)):
            if clusters[y] not in temp:
                new_clusters.append(clusters[y])
            else:
                if "sub_clust_"+str(counter) not in new_clusters:
                    new_clusters.append("sub_clust_"+str(counter))
        return new_clusters

    # An iterator to sequentially merge sequences/sub-clusters with shortest KS.
    while len(clusters) > 1:
        # if mode is M2 we update the sequence list so sub-clusters are treated as terminal sequences in next iteration.
        if mode == 'M2':
        	sequences = clusters
        counter = counter + 1
        current_KS = get_current_KS()
        min_KS = get_min_KS(current_KS)
        clusters = update_clusters(clusters,min_KS)
    return KS

# Open output files.
OUT1 = open(output_file+"_KS_by_cluster.txt", "w")
OUT2 = open(output_file+"_KS.txt", "w")

# Call main function and write results to outputs.
for key in data_dict:
    KS = remove_redundant_KS(data_dict[key])
    print >> OUT1,key,KS
    for i in xrange(0,len(KS)):
    	print >> OUT2,KS[i]
    
# Close output files.
OUT1.close()
OUT2.close()
    
    	

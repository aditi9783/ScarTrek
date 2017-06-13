#!/usr/bin/env python

import sys

# start of extractSeq ##########################
def extractSeq( fh ): # extract gene or protein seq from respective file handles
    seqdict = {}
    seq = ""
    for line in fh:
        if ">" in line: # header line
            if seq != "": # there is some sequence to save
                seqdict[gname] = seq.replace('\n','')
            content = line.split(' [')
            genestring = content[1].split('=')
            gname = genestring[1].rstrip(']')
            seq = ""
        else:
            seq += line
    return seqdict
# end of extractSeq ##########################

genes = [] # list of tuples of gene start, gene end, gene name such that start < end (can't identify complement genes)
gdict = {} # key: genename, value: gene start and end pos with start > end for complement genes
ntseq = {} # key: genename, value: nt seq from NCBI
aaseq = {} # key: genename, value: amino acid seq for the encoded protein (from NCBI)
refgeneloci = sys.argv[1] # paths where reference genome's gene info, gene seq, and protein seq are located
refgeneseq = sys.argv[2]
refprotseq = sys.argv[3]
genefh = open(refgeneloci, 'r') # H37Rv genes
ntgenesfh = open(refgeneseq, 'r') # nt seq for genes
aagenesfh = open(refprotseq, 'r') # aa seq for genes

for line in genefh: # genes are already sorted by position in this genefile
    line = line.rstrip("\n")
    content = line.split()
    start = int(content[1])
    end = int(content[2])
    # gene nt and aa seq files have gene names by default, and locus tag when gene name is absent
    # save gene start-end pos using this naming convention
    gnamelist = content[3].split('_')
    gname = ""
    if gnamelist[0] != '': # both the gene name and locus tag is present
        if gnamelist[0] == "PE": # PE_PGRS genes
            gname = gnamelist[0]+"_"+gnamelist[1]
        else:
            gname = gnamelist[0]
    else:
        gname = gnamelist[1] # only locus tag is present

    gdict[gname] = [start, end]

    if start < end:
        genes.append([ start, end, gname ])
    else:
        genes.append([ end, start, gname ])
genefh.close()
# extract nucleotide and protein sequences for H37Rv genes
ntseq = extractSeq( ntgenesfh )
aaseq = extractSeq( aagenesfh )
ntgenesfh.close()
aagenesfh.close()

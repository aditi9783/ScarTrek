#!/usr/bin/env python

import argparse
import sys
from scartrek.scars_functions import *
from scartrek.globalvars import extractSeq


# start of scars ###########################
def scars( slist, MAPRATE, COVTHRES, ntseq, gdict, genes, aaseq): # Input: tuple of [seqdir, sd], seqdir: full path to mapped reads for each strain, sd: strain name
    seqdir, sd = slist
    maprate = checkMappingRate( seqdir+"/mapped/" ) # Get the % of reads that mapped to the reference
    print((seqdir, maprate))
    if maprate > MAPRATE: # only proceed if 20% (default) or more reads mapped
        findIndels( seqdir+"/mapped/", sd, COVTHRES ) # find all indels from mpileup file, and write those in the mapped folder for each strain
        with open(seqdir+"/mapped/"+sd+".genewise.mutations2", 'w') as fh,\
                open(seqdir+"/mapped/"+sd+".stopcodon_causing_mut2", 'w') as fhstop,\
                open(seqdir+"/mapped/"+sd+".frameshifts2", 'w') as fhfs:
            genewise_indels = getStrainIndels(seqdir+"/mapped/", sd, genes, gdict.keys())
            for g in genewise_indels:
                # not all genes had their sequences in NCBI sequence files that were read into ntseq and aaseq
                # if a gene sequence is not available, skip
                if g not in ntseq.keys():
                    continue
                nummut = len(genewise_indels[g])
                if nummut > 0: # at least one mutation
                    # write genic mutations to file
                    fh.write(g+" : "+" : "+" ".join(genewise_indels[g])+"\n")
                if nummut > 1: # gene should have > 1 mutations/indels to check for scars or reverting stop codons
                    #print g, genewise_indels[g]
                    if g not in aaseq.keys(): # RNA genes like rrs won't have a protein sequence
                        continue # skip if no protein sequence
                    processedIndels, disruptivemut = checkProtein2( g, genewise_indels[g], ntseq, gdict ) # check the impact of multiple mutations on the protein
                    if len(disruptivemut) > 0: # write stop codon inducing mutations to file
                        fhstop.write(g+" : "+" ".join(disruptivemut)+"\n")
                    output = checkFrameshifts( g, processedIndels, ntseq, aaseq, sd )
                    fhfs.write(output)

    return
# end of scars ###########################


def main(argv = sys.argv):
    parser = argparse.ArgumentParser(description='Detect indel scars from mpileup files.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=True,
                        help='Input directory that has mpileup files for each sample. See tests/test1 for example')
    parser.add_argument('-m', '--maprate', default=20.0, type=float,
                        help='Minimum read mapping rate required to consider a sample')
    parser.add_argument('-c', '--covthres', default=20, type=int,
                        help='Minimum read coverage required at a position to detect an indel')
    parser.add_argument('-g', '--geneseq', default="../reference/H37Rv_genes.txt",
                        help='Gene sequences in the reference genome, default reference: M. tuberculosis')
    parser.add_argument('-p', '--protseq', default="../reference/H37Rv_proteins_from_genbank.txt",
                        help='Protein sequences for the reference organism, default: M. tuberculosis')

    args = parser.parse_args(argv)

    seqpath = args.input
    seqids = os.walk(seqpath).next()[1] # only get the subdir in the seqpath, not further down subdirs

    # seqlist = [ x for sd in seqids]
    # seqlist = [ (seqpath + sd, sd) for sd in seqids ]
    seqlist = []
    for sd in seqids:
        #if len(seqlist) == 2: # for testing
        #    break
        seqdir = seqpath+sd
        seqlist.append( [seqdir, sd] )

    genes = []  # list of tuples of gene start, gene end, gene name such that start < end (can't identify complement genes)
    gdict = {}  # key: genename, value: gene start and end pos with start > end for complement genes
    ntseq = {}  # key: genename, value: nt seq from NCBI
    aaseq = {}  # key: genename, value: amino acid seq for the encoded protein (from NCBI)

    # extract nucleotide and protein sequences for H37Rv genes
    # ntseq, gdict, genes = extractSeq("../reference/H37Rv_genes.txt", 0)
    # aaseq = extractSeq("../reference/H37Rv_proteins_from_genbank.txt", 1)
    ntseq, gdict, genes = extractSeq(args.geneseq, 0)
    aaseq = extractSeq(args.protseq, 1)
    for sq in seqlist:
        #print sq
        scars(sq, args.maprate, args.covthres, ntseq, gdict, genes, aaseq)


if __name__ == '__main__':
    main(sys.argv)

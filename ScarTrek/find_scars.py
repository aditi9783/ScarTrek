#!/usr/bin/env python

import os, subprocess, sys
from multiprocessing import Pool
from scars_functions import *
from globalvars import genes, gdict, ntseq, aaseq

# start of scars ###########################
def scars( slist ): # Input: tuple of [seqdir, sd], where seqdir is full path to mapped reads for each strain, and sd is the strain name
    seqdir, sd = slist
    maprate = checkMappingRate( seqdir ) # Get the % of reads that mapped to the reference
    print seqdir, maprate
    if maprate > 20.0: # only proceed if 20% or more reads mapped, so that we have sites with coverage >= 20 for downstream analyses
        findIndels( slist ) # find all indels from mpileup file, and write those in the mapped folder for each strain 
        fh = open(seqdir+"/mapped/"+sd+".genewise.mutations", 'w')
        fhstop = open(seqdir+"/mapped/"+sd+".stopcodon_causing_mut", 'w')
        fhfs = open(seqdir+"/mapped/"+sd+".frameshifts", 'w')
        genewise_indels = getStrainIndels(seqdir, sd, genes, gdict.keys())
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
                processedIndels, disruptivemut = checkProtein2( g, genewise_indels[g], ntseq, gdict ) # check the impact of multiple mutations on the protein
                if len(disruptivemut) > 0: # write stop codon inducing mutations to file
                    fhstop.write(g+" : "+" ".join(disruptivemut)+"\n")
                output = checkFrameshifts( g, processedIndels, ntseq, aaseq, sd )
                fhfs.write(output)
        fh.close()
        fhstop.close()
        fhfs.close()
    return 
# end of scars ###########################

if __name__ == '__main__':
    seqlist = []
    seqpath = sys.argv[1]
    seqids = os.walk(seqpath).next()[1] # only get the subdir in the seqpath, not further down subdirs
    for sd in seqids:
        #if len(seqlist) == 2: # for testing
        #    break   
        seqdir = seqpath+sd
        seqlist.append( [seqdir, sd] )

    for sq in seqlist:
        #print sq
        scars(sq)
    #print(seqlist)

    exit()
    pool = Pool(15) # parallel processors
    pool.map(scars, seqlist)
    pool.close()
    pool.join()

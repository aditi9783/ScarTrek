#!/usr/bin/env python

import os
from scars_functions import *
from globalvars import genes, gdict, ntseq, aaseq
import filecmp


# start of scars ###########################
def scars( slist ): # Input: tuple of [seqdir, sd], seqdir: full path to mapped reads for each strain, sd: strain name
    seqdir, sd = slist
    maprate = checkMappingRate( seqdir+"/mapped/" ) # Get the % of reads that mapped to the reference
    print seqdir, maprate
    if maprate > 20.0: # only proceed if 20% or more reads mapped
        findIndels( seqdir+"/mapped/", sd ) # find all indels from mpileup file, and write those in the mapped folder for each strain
        fh = open(seqdir+"/mapped/"+sd+".genewise.mutations2", 'w')
        fhstop = open(seqdir+"/mapped/"+sd+".stopcodon_causing_mut2", 'w')
        fhfs = open(seqdir+"/mapped/"+sd+".frameshifts2", 'w')
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
        fh.close()
        fhstop.close()
        fhfs.close()

    return maprate
# end of scars ###########################

if __name__ == '__main__':
    seqlist = []
    # seqpath = "/scratch/ag1349/Zhang2013NatureGen/NGS_analysis/"
    seqpath = "/Users/aditi/PycharmProjects/ScarTrek/scartrek/tests/test1/"

    # <inputDir>/<sampleNames>/mapped
    seqids = os.walk(seqpath).next()[1] # only get the subdir in the seqpath, not further down subdirs
    for sd in seqids:
        #if len(seqlist) == 2: # for testing
        #    break   
        seqdir = seqpath+sd
        onlyfiles = [f for f in os.listdir(seqdir+"/mapped") if os.path.isfile(seqdir+"/mapped/"+f)]
        seqlist.append( [seqdir, sd] )

    # print seqlist
    # maprate = {seqdir: scars((seqdir, sd)) for seqdir, sd in seqlist}
    # print maprate
    # assert maprate['/Users/aditi/PycharmProjects/ScarTrek/scartrek/tests/test1/sample1', 99.21]
    # for sq in seqlist:
    #     #print sq
    #     scars(sq)

    samples = [ "sample1", "sample2", "sample3" ]
    files_to_compare = [ "frameshifts2", "indels", "genewise.mutations2", "stopcodon_causing_mut2"]

    scar_results = {}
    for sample in samples:
        for filename in files_to_compare:
            sample_basedir = seqpath + "/" + sample + "/mapped"

            expected_file = sample_basedir + "/oldfiles/" + sample + "." + filename
            actual_file = sample_basedir + "/" + sample + "." + filename

            scar_results["{} {}".format(sample, filename)] = filecmp.cmp(expected_file, actual_file)

    failed_tests = [sample for sample, result in scar_results.items() if not result]
    if failed_tests:
        print "Failed "
        print [sample for sample in failed_tests]
    else:
        print "All passed!"

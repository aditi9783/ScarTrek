#!/usr/bin/env python

import os, subprocess, sys
from multiprocessing import Pool

def zipmpileup( seqdir ):
    os.system("gzip "+seqdir+"/aln.sorted.bam.mpileup") # compress mpileup output

if __name__ == '__main__':
    seqpath = sys.argv[1] 
    seqids = os.walk(seqpath).next()[1] # only get the subdir in the seqpath, not further down subdirs
    unfinished = []
    for sd in seqids:
        seqdir = seqpath+sd+"/mapped"
        if os.path.exists(seqdir): # if this directory exists
            allf = os.listdir(seqdir)
            for seqf in allf: # get all files
                if seqf == "aln.sorted.bam.mpileup": # uncompressed mpileup file present
                    for sn in allf: # check if incomplete compressed file is present as well
                        if "aln.sorted.bam.mpileup.gz" == sn:
                            os.system("rm "+seqdir+"/aln.sorted.bam.mpileup.gz") # delete incompletely compressed file
                            break
                    unfinished.append(seqdir)    
    print "Unfinished ", len(unfinished)

    pool = Pool(15) # parallel processors
    pool.map(zipmpileup, unfinished)
    pool.close()
    pool.join()

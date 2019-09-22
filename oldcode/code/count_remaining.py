#!/usr/bin/env python

import sys, os, subprocess

if __name__ == '__main__':
    seqlist = []
    seqpath = sys.argv[1] 
    seqids = os.walk(seqpath).next()[1] # only get the subdir in the seqpath, not further down subdirs
    unfinished = []
    for sd in seqids:
        #if len(seqlist) == 2: # for testing
        #    break   
        seqdir = seqpath+sd
        #print("\n",seqdir)
        allf = os.listdir(seqdir)
        for seqf in allf: # get paired end files
            if "_1.fastq.gz" in seqf: # raw fastq files were deleted after processing was finished. If these files are present => un processed
                unfinished.append( [seqdir, sd] )
    print "Number of unprocessed samples:", len(unfinished)

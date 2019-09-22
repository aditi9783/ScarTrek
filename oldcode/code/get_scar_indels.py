#!/usr/bin/python

import sys

fname = sys.argv[1] # all_frameshifts.out file
datadir = sys.argv[2]
foutname = fname+".res" # change it to more informative suffix

fh = open(fname, 'r') # all_frameshifts.out file
fhout = open(foutname, 'w')
genewise_fshifts = {} # key: gene name, value: frame shifts
for line in fh:
    contents = line.split('\t')
    gname = contents[1]
    isolatelist = contents[0].split(':')
    isolate = isolatelist[1]
    if gname in genewise_fshifts: # if gene key exists
        genewise_fshifts[gname].append( isolate )
    else: # create a new key for the gene name
        genewise_fshifts[gname] = [isolate]
fh.close()

for g in genewise_fshifts:
    n = len(genewise_fshifts[g])
    print g, "Number of fshifts:", n
    indels = {} # key: indel line, value: list of isolates with that indel line
    for isolate in genewise_fshifts[g]:
        fh2 = open(datadir+isolate+"/mapped/"+isolate+".frameshifts", 'r')
        gflag = 0
        for line in fh2:
            if g in line: 
                gflag = 1
            elif "Processed indels:" in line and gflag == 1: # print processed indels for the gene
                if line in indels: # if indel line already exists
                    indels[line].append( isolate )
                else:
                    indels[line] = [isolate]
                break
        fh2.close()
    for l in indels:
        print "\t", l, "\t", indels[l], "\n\n"
        for isolate in indels[l]: # get all isolates that have this indel
            fhout.write(datadir+isolate+"/mapped/"+"\t"+isolate+"\t"+g+"\t"+l)
fhout.close()
    

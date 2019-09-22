#!/usr/bin/env python

import os, subprocess
import sys
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams.update({'font.size': 10})

# start of findMutations #################
def findMutations( seqlist ):
    mutmatrix = []
    strainids = []
    #fout = open("distmat_strain_ids.out", 'w') # write the strain ids as key for the distance matrix
    counter = 0
    for sq in seqlist:
        seqdir, sd = sq
        strainids.append(sd)
    #    fout.write(str(counter)+"\t"+sd+"\n")
        mutlist = []
        fh = open(seqdir+"/mapped/"+sd+".vcf", 'r')
        for line in fh:
            if line.startswith("#"):
                    continue
            else:
                content = line.split("\t")
                if float(content[5]) > 100.0: # minimum quality score for calling a mutation
                    if len(content[3]) == 1 and len(content[4]) == 1: # both the ref base and mutated base are of len 1 => no indels, only SNPs 
                        mut = content[1]+content[4]
                        mutlist.append( mut )
        fh.close()
        mutmatrix.append( mutlist )
    #fout.close()
    return mutmatrix, strainids

# end of findMutations ###################

# start of calcDistances #################
def calcDistances( mutmatrix, strainids ):
    distmat = []
    for i in range(len(mutmatrix)):
        distmat.append([0 for x in range(len(mutmatrix))]) # initialize the distances
    
    for i in range(len(mutmatrix)):
        for j in range(i+1, len(mutmatrix)):
            intersection = list( set(mutmatrix[i]) & set(mutmatrix[j]) )
            union = list( set(mutmatrix[i]) | set(mutmatrix[j]) )
            num_shared_mut = len(intersection) # number of mutations common in both strains
            num_all_mut = len(union) # number of unique mutations present in both strains
            if num_all_mut == 0:
                print strainids[i], strainids[j]
                dist = 1.0
            else:
                dist = 1.0 - float(num_shared_mut)/float(num_all_mut)
            distmat[i][j] = round(dist, 3)
            #distmat[i][j] = num_shared_mut
            distmat[j][i] = distmat[i][j] # because the matrix is symmetric
    return distmat
    
# end of calcDistances ###################

# start of printPhylipDistMat ############
def printPhylipDistMat( strainids, distmat ):
    outfname = sys.argv[1]+"/distmat_phylip.txt"
    fout = open(outfname, 'w')
    numstrains = len(strainids)
    fout.write("\t"+str(numstrains)+"\n")

    for i in range(numstrains):
        if len(strainids[i]) == 9: # ENA acc numbers are either 9 or 10 in length. Phylip requires that strain names be of length 10
            fout.write(strainids[i]+" "+"\t"+"\t".join([str(v) for v in distmat[i]])+"\n")
        else: # strain name is of length 10, no space padding needed
            fout.write(strainids[i]+"\t"+"\t".join([str(v) for v in distmat[i]])+"\n")

    fout.close()
# end of printPhylipDistMat ############

# start of pairedDistances #############
def pairedDistances( distmat, strainids, pairedsamples ):
    pairedidx = [] # strainid indices that represent paired samples
    paireddist = [] # distance between paired samples
    unpaireddist = [] # distance between unpaired samples
    for pair in pairedsamples:
        pids = []
        for i in range(len(strainids)):
            if strainids[i] in pair: # if strainid belongs to the current pair
                pids.append( i )
                if len(pids) == 2: # both strain ids in the pair has been found
                    pid1, pid2 = pids
                    paireddist.append( distmat[pid1][pid2] )
                    print(strainids[pid1], "\t", strainids[pid2], "\t", distmat[pid1][pid2], "\n")
                    pid1dist = [item for i,item in enumerate(distmat[pid1]) if i > pid1 and i not in pids] # get distance of 1st strain in pair to other unpaired strains (upper triangle to avoid counting distances more than once)
                    pid2dist = [item for i,item in enumerate(distmat[pid2]) if i > pid2 and i not in pids] # get distance of 2nd strain in pair to every other unpaired strain
                    unpaireddist.extend( pid1dist )
                    unpaireddist.extend( pid2dist )
                    break # go to the next pair
    return paireddist, unpaireddist 
# end of pairedDistances #############

if __name__ == '__main__':
    seqlist = []
    seqids = []
    pairedsamples = []
    seqpath = sys.argv[1] 
    if len(sys.argv) > 1: # if a file is provided with isolate ids
        fh = open(sys.argv[1], 'r')
        for line in fh:
            line = line.rstrip("\n")
            content = line.split("\t")
            seqids.extend( content ) # add both pairs to the seqid list 
            pairedsamples.append( content ) 
    else: # if no file with isolate ids is provided, process all isolates present in the given path
        seqids = os.walk(seqpath).next()[1] # only get the subdir in the seqpath, not further down subdirs

    for sd in seqids:
        seqdir = seqpath+sd

        # since problematic cases were moved to other location to allow running of parallel computing, check if the path exists
        if os.path.isdir(seqdir):
            allf = os.listdir(seqdir+"/mapped")
            if len(allf) > 8: # Analyses was not done for isolates with poor coverage. If analyses was done, it generated > 8 files.
                if sd+".vcf.gz" in allf: # if the compressed .vcf file is present, uncompress it
                    os.system("gzip -d "+seqdir+"/mapped/"+sd+".vcf.gz")
                seqlist.append( [seqdir, sd] )
        else:
            print("Problem case:", sd, "\n")

    mutmatrix, strainids = findMutations( seqlist )
    distmat = calcDistances( mutmatrix, strainids )
    printPhylipDistMat( strainids, distmat ) # print in the distance matrix format read by Phylip

    exit() 
    # the following code applies only for paired samples
    pairdist, nonpairdist = pairedDistances( distmat, strainids, pairedsamples ) # return distances for paired samples, and distances between unpaired samples
    #print(pairedsamples)
    #print("pair dist ", pairdist)
    #print("nonpair dist ", nonpairdist)

    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].hist(pairdist, bins=100, color='gray')
    axarr[0].set_ylabel("Number of Linked Isolates", fontsize=9)
    axarr[1].hist(nonpairdist, bins=100, color='gray')
    axarr[1].set_ylabel("Number of Unlinked Isolates", fontsize=9)
    axarr[1].set_xlabel("Pairwise Distances", fontsize=9)
    f.savefig("isolate_distances.pdf")
    plt.close()

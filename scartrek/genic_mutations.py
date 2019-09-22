#!/usr/bin/env python

# read complete DST profile file and retain only those columns (mutations) that have non-synonymous genic mutations

from functions_genic_mutations import *

# open necessary files
dstfh = open("/home/ag1349/Hassan/Zhangetal2013/results/DST_profiles_substitutions.txt", 'r') # DST Profiles
genefh = open("/home/ag1349/Hassan/output/figures/genome_plots/data/genes_loci.txt", 'r') # H37Rv genes
ntgenesfh = open("/home/ag1349/Hassan/H37Rv/H37Rv_genes.txt", 'r') # nt seq for genes
aagenesfh = open("/home/ag1349/Hassan/H37Rv/H37Rv_proteins.txt", 'r') # aa seq for genes

# read genic regions
genes = [] # list of tuples of gene start, gene end, gene name
gdict = {} # key: genename, value: gene start and end pos
for line in genefh: # genes are already sorted by position in this genefile
    line = line.rstrip("\n")
    content = line.split()
    start = int(content[1])
    end = int(content[2])
    # gene nt and aa seq files have gene names by default, and locus tag when gene name is absent
    # save gene start-end pos using this naming convention
    genenamelist = content[3].split('_')
    gname = ""
    if genenamelist[0] != '': # both the gene name and locus tag is present
        gname = genenamelist[0]
    else:
        gname = genenamelist[1] # only locus tag is present

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

############ MAIN #########################

# do this test only once: check if all nt sequences make the correct amino acid sequences
#for tup in genes:
#    if tup[2] not in ntseq: # tRNA genes are not in ntseq and aaseq. So ignore those
#        continue
#    gene_nt = ntseq[tup[2]]
#    gene_aa = aaseq[tup[2]]
#    trans_aa = translateGene( gene_nt )
#    trans_aa = trans_aa.rstrip('*')
#    if trans_aa[1:] == gene_aa[1:]: # since first residue is sometimes not changed to Met
#        print tup[2], "- ok"
#    else:
#        print tup[2], "- Check sequence!"
# TEST PASSED. When comparing seq, ignore 1st residue.

#gid = "Rv0008c"
#print gdict[gid]
#print ntseq[gid]
#print aaseq[gid]
#exit()
mut = dstfh.readline().rstrip("\n") # read mutation header in DST file: has all mutations found in all strains
allmutations = mut.split()
num_mut = len(allmutations)-1
mutdecision = [1] # decision for whether a mutation is non-synonymous (1) or not (0). First entry is word "Strain", so that is given a 1 to preseve that column. Columns with value 0 will be removed. 

for i in range(1, num_mut+1):
#for i in range(1, 30): # read 1st five mutations, for testing
    mutdecision.append(checkMut( allmutations[i], genes, ntseq, aaseq, gdict ))

# write all genic mutations to a new file
dstoutfh = open("/home/ag1349/Hassan/Zhangetal2013/results/complete_DST_profiles_only_genicmut.txt", 'w')
nonsynmutations = []
for i in range(len(mutdecision)):
    if mutdecision[i] == 1:
        nonsynmutations.append(allmutations[i])
dstoutfh.write("\t".join(nonsynmutations)+"\n")

# get all the lines in the dst profile file, and save the nonsyn genic mutations only
mutmatrix = []
for line in dstfh:
    line = line.rstrip("\n")
    content = line.split()
    pruned = []
    for i in range(len(mutdecision)):
        if mutdecision[i] == 1:
            pruned.append( content[i] )
    dstoutfh.write("\t".join(pruned)+"\n") 
#    print content[0], len(pruned)
    mutmatrix.append(pruned)
dstoutfh.close()

strain_mutations = {} # key1: strain name, key2: "subs"/"indels", key3: gene name, value: list of "subs"/"indels" in that gene
genemutfreq = {} # key: gene name, value: number of times that gene was mutated in the dataset
for g in ntseq.keys(): # get gene names
    genemutfreq[g] = 0 # initialize the frequency hash
 
for s in mutmatrix:
    sname = s[0] # strain name
    fh = open("/home/ag1349/Hassan/Zhangetal2013/NGS_analysis/"+sname+"/mapped/"+sname+".genewise.mutations", 'w')
    print "\n=======", sname, "\n"
    genenames = ntseq.keys()
    genewise_subs, mutpos = getStrainSubs(genenames, [int(val) for val in s[1:]], genes, nonsynmutations[1:]) # pass genenames, and mutkey. 0th idx is "Strain", don't pass that
    s_cov = getCoverage( sname, mutpos, 20 ) # get coverage at mut sites. If mutation site is not in s_cov, coverage is less than the value supplied. 
    genewise_subs_pruned = prunedMutlist( genewise_subs, s_cov ) # get only those subs that have at least thres coverage
    # get indels for this strain
    genewise_indels = getStrainIndels( sname, genenames, genes )
    for g in genewise_subs_pruned:
        nummut = len(genewise_subs_pruned[g])+len(genewise_indels[g])
        if nummut > 0: # at least one mutation
            # write genic mutations to file
            fh.write(g+" : "+" ".join(genewise_subs_pruned[g])+" : "+" ".join(genewise_indels[g])+"\n")
        if nummut > 1: # gene should have > 1 mutations to check for scars: reverting stop codons
            #print g, genewise_subs_pruned[g], genewise_indels[g]
            processedIndels = checkProtein2( g, genewise_subs_pruned[g], genewise_indels[g], ntseq, gdict ) # check the impact of multiple mutations on the protein
        if len(genewise_indels[g]) > 1: # at least two indels for frameshift
        #    print g, genewise_subs_pruned[g], genewise_indels[g]
        #    processedIndels = checkProtein2( g, genewise_subs_pruned[g], genewise_indels[g], ntseq, gdict ) # check the impact of multiple mutations on the protein
            checkFrameshifts( g, ntseq[g], aaseq[g], processedIndels)
    fh.close()
    #if sname == "I-xz09061": # just for testing
    #    exit() 

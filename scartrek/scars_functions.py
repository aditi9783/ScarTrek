#!/usr/bin/python 

import os
from genetic_code import code
#from globalvars import genes, gdict, ntseq, aaseq

global COVTHRES
COVTHRES = 20 # min coverage for detecting indels

# define base complements
complement = {"A" : "T", "T" : "A", "C" : "G", "G" : "C"}

# start of checkMappingRate ####################
def checkMappingRate( fpath ):
    fh = open(fpath+"map.run", 'r')
    maprate = "0.0"
    for line in fh:
        if " overall alignment rate" in line:
            contents = line.split(" ")
            maprate = contents[0].rstrip("%")
    return float(maprate)
# end of checkMappingRate ####################

# start of findIndels ########################
def findIndels( fpath, sd ):

    mpileupfh = open(fpath+"aln.sorted.bam.mpileup", 'r')
    outfh = open(fpath+sd+".indels", 'w')
    n_idl = 0
    for line in mpileupfh:
        content = line.split()
        # if the mpileup file is messed up (less than 5 items in line), print the dir name and move on to the next file
        if len(content) > 3: # positions with no coverage has len(content) = 4, thus correct lines have len > 3
            try:
                int(content[3]) # and correct files have 4th column as the coverage- an integer
            except ValueError: # if the 4th column does not have an int -> raise exception
                print "Incorrect mpileup:", fpath, sd, line
                break
        else: # the line has less than 4 elements, incorrect line.
            print "Incorrect mpileup:", sd, fpath, line
            break

        cov = int(content[3])
        if cov > COVTHRES: # at least COVTHRES depth at this position
            readbases = content[4].upper() # contains bases/indels in the read. Conver to uppercase so that 1T and 1t are read as same indel
            if "+" in readbases or "-" in readbases: # indel 
                bases = list(readbases)
                indel = []
                for i in range(len(bases)):
                    if bases[i] == "+" or bases[i] == "-":
                        if i-1 >=0 and bases[i-1] == "^": # ^ is a char that can precede +/- in non-indel scenario
                            continue
                        indellen = int(bases[i+1])
                        digitlen = 1 # indel may run into double or higher digits, increase indel span based on that
                        if bases[i+2].isdigit(): # indel len runs in double digits
                            indellen = int("".join(bases[i+1:i+3]))
                            digitlen = 2
                        if i+3 < len(bases) and bases[i+3].isdigit(): # indel len runs in triple digits
                            indellen = int("".join(bases[i+1:i+4]))
                            digitlen = 3
                        indel.append( "".join(bases[i:i+1+digitlen+indellen]) ) # save from indelsign (+1), len (digitlen), and end of indel (indellen+1)
                indelcounts = {x:indel.count(x) for x in indel} # key: indel, value: freq of that indel
                for idl in indelcounts:
                    if indelcounts[idl] > cov/2: # at least half the reads should have this indel
                        n_idl += 1
                        #print line
                        #print idl, indelcounts[idl]
                        outfh.write(content[1]+"\t"+content[2]+"\t"+idl+"\n")
            else: # only interested in indels
                continue
    # print "Number of true indels:", n_idl
    mpileupfh.close()
    outfh.close()
# end of findIndels ########################

# start of checkRegion ##########################
def checkRegion( pos, genes ): # check if this position is in genic region
    glist = []
    for tup in genes:
        if pos >= tup[0] and pos <= tup[1]: # pos is within gene
            glist.append(tup[2])

    if len(glist) == 0: # no genes found that contained this mutation
        return 0, glist
    else:
        return 1, glist
# end of checkRegion ##########################

# start of translateGene ##########################
def translateGene( nt ):
    tseq = ["M"] # no matter the first codon, the amino acid seq always starts with M -> NOT TRUE! Few proteins in H37Rv don't start with M
    for i in range(3, len(nt), 3): # start from 2nd codon, and then read gene sequence in steps of 3
        codon = "".join(nt[i:i+3])
        if codon in code:
            tseq.append( code[codon] )
        else: # some indels change frame, so last codon might be incomplete
            tseq.append( "***" )
    trans_aa = "".join(tseq)
    trans_aa = trans_aa.rstrip('*') # remove any stop codons at the end
    return trans_aa
# end of translateGene ##########################

# start of findStopCodon ##########################
def findStopCodon( gseq ):
    trans_aa = translateGene( gseq ) # translate mutated protein
    if '*' in trans_aa: # if a stop codon is inserted due to mutation; or frameshift led to a codon < 3 in size => insert "***" in that case
        #print "STOP CODON "
        #print trans_aa
        return 1
    else: # no stop codon
        return 0
# end of findStopCodon ##########################

# start of addMutations ##########################
def addMutations( genename, gseq, indels ): # add all indels together, and return the mutated seq
    mutseq = gseq[:]
    # add all indels    
    sorted_idls = sorted(indels,key=lambda x: x[1]) # sort by 2nd column: mutation site
    shift = 0 # shift in reference position due to addition of indels
    for idl in sorted_idls: # idl is tuple: [index for indel start, wt bases, indel bases, index for end of wt bases that correspond to the indel]
        mutpos = idl[1] + shift
        if idl[0] == "+": # insertion
            mutseq = mutseq[0:mutpos+1] + idl[2] + mutseq[mutpos+1:] # 1st slice goes to mutsite+1 so that wt base is also included
            shift += len(idl[2]) # len of insert adds to the shifted index 
        elif idl[0] == "-": # deletion
            mutseq = mutseq[0:mutpos+1] + mutseq[mutpos+1+len(idl[2]):]
            shift -= len(idl[2]) # len of delete subtracts from the shifted index 
    return mutseq # return mutated seq
# end of addMutations ##########################

# start of checkProtein2 ##########################
def checkProtein2(gname, indellist, ntseq, gdict):
    gene_nt = list(ntseq[gname])
    gstart, gend = gdict[gname]
    stop_codon_flag = 0 # if any of the mutations lead to a stop codon in the gene, put all of the 
    indels = []
    disruptive_mut = [] # mutations that introduce stop codon

    # check if any indels introduce a stop codon
    for idl in indellist:
        m_wt, m_pos, m_mut = idl.split('_') 
        mutsite = int(m_pos) - gstart
        indel = list(m_mut)
        indelbases = [b for b in indel[1:] if b.isalpha()] # get only bases, not indel sign and length
        m_mut = indelbases[:] # copy by value. indel bases for normal (i.e. NOT complement) genes.
        mutated_nt = []
        if gstart > gend: # complement gene
            comp_m_mut = [complement[base] for base in indelbases]
            revcomp_m_mut = comp_m_mut[::-1] # reverse the complement string
            m_wt, m_mut, mutsite = [complement[m_wt], revcomp_m_mut, -mutsite] 
        #print "subs: refbase:", gene_nt[mutsite], "\t mutation:", m_wt, m_mut, idl, mutsite
        #print "refregion:", gene_nt[max(mutsite-5,0):mutsite+len(m_mut)+5]
        if gene_nt[mutsite] == m_wt: # reference base matches
            if indel[0] == "+": # an insertion
                indels.append(["+", mutsite, m_mut])
                mutated_nt = gene_nt[0:mutsite+1] + m_mut + gene_nt[mutsite+1:] # 1st slice goes to mutsite+1 so that wt base is also included
                #print "1st slice:", gene_nt[mutsite-10:mutsite+1], "\n", "insert seq:", m_mut, "\n2nd slice:", gene_nt[mutsite+1:mutsite+1+len(m_mut)+5]
            elif indel[0] == "-": # deletion
                indels.append(["-", mutsite, m_mut])
                deletionseq = gene_nt[mutsite+1:mutsite+1+len(m_mut)]
                mutated_nt = gene_nt[0:mutsite+1] + gene_nt[mutsite+1+len(m_mut):]
                #print "1st slice:", gene_nt[mutsite-10:mutsite+1], "\n", "deletion seq:", deletionseq, "\n2nd slice:", gene_nt[mutsite+1+len(m_mut):mutsite+1+len(m_mut)+5]
        #print "mutated nt region:", mutated_nt[max(mutsite-5,0):mutsite+len(m_mut)+5]
        stop_codon_flag = findStopCodon( mutated_nt ) # translate mutated protein and look for stop codons
        if stop_codon_flag == 1:
            disruptive_mut.append(idl)

    # print stop codon inducing indels to file
    #if len(disruptive_mut) > 0:
    #    print gname, "\tStop-codon introducing mutations: ", disruptive_mut

    return indels, disruptive_mut

# end of checkProtein2 #########################

# start of getStrainIndels #####################
def getStrainIndels( fpath, sd, genes, gnames ):
    genemut = {gn : [] for gn in gnames} # key: gene, val: list of mutations in that gene
    fh = open(fpath+sd+".indels", 'r')
    for line in fh:
        line = line.rstrip('\n')
        content = line.split()
        mpos = int(content[0])
        mflag, genes_mut = checkRegion( mpos, genes ) # find genes in which this indel resides
        for g in genes_mut:
            if g in genemut: # only consider coding genes as per H37Rv gene and protein files from NCBI
                genemut[g].append( content[1]+"_"+content[0]+"_"+content[2] ) # refbase_pos_indel 
    fh.close()
    return genemut 
# end of getStrainIndels #####################

# start of checkFrameshifts ######################
def checkFrameshifts( gname, pindels, ntseq, aaseq, sd ): # these indels are processed indels returned by checkProtein2 function
    # if indel len is in multiples of 3 => not frameshift indels
    ilenflag = 0
    output = "\n"+gname+":\nProcessed indels: "
    for idl in pindels:
        output = output+" ".join([str(x) for x in idl])
        output = output+"\t"
    output = output+"\n"
    #print("Processed indels:", gname, pindels)
    for idl in pindels: # each indel is a tuple of 3 elements, with last element being the indel seq
        if len(idl[2]) % 3 != 0: # indel size is not a multiple of 3 => frameshift
            ilenflag = 1
            break # exit this for loop
    if ilenflag == 0: # none of the indels had lengths other than in multiples of 3 => not a frameshift
        output  = output+"Indels multiple of 3 => no frameshift. Exiting.\n\n"
        return output
    gene_nt = ntseq[gname]
    gene_aa = aaseq[gname]
    #if gname.endswith("c"):
    #    print gname, gene_aa
    #print gname
    #print pindels
    mutatedseq = addMutations( gname, list(gene_nt), pindels )
    if findStopCodon(mutatedseq) == 0: # no stop codon
        output  = output+"\n===No more stop codon after adding all indels!\n\n"
        output  = output+gname+"\t checking for frameshifts.\nOriginal gene seq:\n"+gene_nt
        output  = output+"\nmutated seq:\n"+"".join(mutatedseq)
        translated = translateGene(mutatedseq)
        output  = output+"\nOriginal\n"+gene_aa+"\nTranslated\n"+translated+"\n"
        first_mismatch = 0 # first residue from the beginning of the protein that is mismatched
        last_mismatch_end = 0 # last residue (from end) that is mismatched between mutated and wt proteins
        # compare wt and mutated protein sequences from both ends to see if the frame was restored by the second mutation
        if translated == gene_aa:
            output  = output+"Mutated seq exactly same as wt. Restored.\n"
        else:
            orig = list(gene_aa)
            mut = list(translated)
            #print "orig:", orig, "\nmut:", mut
            for i in range(1,min(len(orig),len(mut))): # ignore 1st residue, as methionine is sometimes posttranslationally modified
            # deletions may reduce the length of the protein, thus min of orig and mut proteins
                if orig[i] == mut[i]:
                    continue
                else:
                    first_mismatch = i
                    output  = output+"\n1st mismatch from L->R: "+str(i)+"\n"
                    break 
            # check for 1st mismatch from the reverse side:
            rev_orig = orig[::-1]
            rev_mut = mut[::-1]
            for i in range(max(len(orig), len(mut))):
                if rev_orig[i] == rev_mut[i]:
                    continue
                else:
                    last_mismatch_end = i
                    output  = output+"1st mismatch from R->L (from end): "+str(i)+"\n"
                    break
        if first_mismatch > 0 and last_mismatch_end > 0:
            output  = output+sd+"\t"+gname+"\t==Frameshift!==\n\n"
    else:
        output  = output+"Stop codon present after adding all indels"
    return output
# end of checkFrameshifts ######################

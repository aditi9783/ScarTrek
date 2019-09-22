#!/usr/bin/env python

from genetic_code import code

# define base complements
complement = {"A" : "T", "T" : "A", "C" : "G", "G" : "C"}

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
def addMutations( genename, gseq, subs, indels ): # add all point mutations and indels together, and return the mutated seq
    mutseq = gseq[:]
    for s in subs: # add all point mutations first. s is a tuple:[mutsite, mut_base]
        mutseq[s[0]] = s[1]
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
    stopflag = findStopCodon( mutseq )
    if stopflag == 0:
        print(genename, "".join(gseq), subs, indels)
        print("\n===No more stop codon after adding all subs+indels!\n\n")
    return mutseq # return mutated seq
# end of addMutations ##########################

# start of checkProtein2 ##########################
def checkProtein2( gname, subslist, indellist, ntseq, gdict ):
    gene_nt = list(ntseq[gname])
    gstart, gend = gdict[gname]
    stop_codon_flag = 0 # if any of the mutations lead to a stop codon in the gene, put all of the 
    subs = []
    indels = []
    disruptive_mut = [] # mutations that introduce stop codon

    # check if any of the substitutions introduce a stop codon
    for m in subslist:
        m_wt, m_pos, m_mut = m.split('_')
        mutsite = int(m_pos) - gstart
        if gstart > gend: # complement gene
            m_wt, m_mut, mutsite = [complement[m_wt], complement[m_mut], -mutsite] 
        #print "subs: refbase:", gene_nt[mutsite], "\t mutation:", m_wt, m_mut, mutsite
        if gene_nt[mutsite] == m_wt: # reference base matches
            gene_nt[mutsite] = m_mut # insert mutation at this position
            subs.append([mutsite, m_mut])
            stop_codon_flag = findStopCodon( gene_nt ) # translate mutated protein and look for stop codons
            if stop_codon_flag == 1:
                disruptive_mut.append(m)
            gene_nt[mutsite] = m_wt # revert the mutation back to wt after detecting the impact of it 

    # check if any indels introduce a stop codon
    for idl in indellist:
        m_wt, m_pos, m_mut = idl.split('_') 
        mutsite = int(m_pos) - gstart
        indel = list(m_mut)
        indelbases = [b for b in indel[1:] if b.isalpha()] # get only bases, not indel sign and length
        m_mut = indelbases[:] # copy by value. indel bases for normal genes.
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

    if stop_codon_flag == 1: # add all mutations in this gene together to see if the stop codon is reversed
        print(gname, "\tStop-codon introducing mutations: ", disruptive_mut)
        mutseq = addMutations( gname, gene_nt, subs, indels ) # add all point mutations and indels together, and return the mutated seq

    return indels

# end of checkProtein2 #########################

# start of checkProtein ##########################
def checkProtein( gname, subslist, indellist, ntseq, gdict ):
    gene_nt = list(ntseq[gname])
    gstart, gend = gdict[gname]
    indels_to_add = [] # list of indels to add. Add those at the end, just so the sequence indices match the reference indices for base substitutions
    substitutions = []
    stop_codon_flag = 0 # if any of the mutations lead to a stop codon in the gene, put all of the 

    # for each mutation, check if the mutation introduces a stop codon 
    for m in mutlist:
        m_wt, m_pos, m_mut = m.split('_')
        mutsite = int(m_pos) - gstart
        if len(m_wt) > 1 or len(m_mut) > 1: # indel
            continue # don't look at indels for now
            if gstart < gend: # normal gene
                wt_bases = gene_nt[mutsite:mutsite+len(m_wt)]
                #print wt_bases, m_wt
                #print "1st slice:", gene_nt[0:mutsite], "\nmutation:", list(m_mut), "\n2nd slice", gene_nt[mutsite+len(m_wt):]                    
                mutated_nt = gene_nt[0:mutsite] + list(m_mut) + gene_nt[mutsite+len(m_wt):]
                stop_codon_flag = findStopCodon( mutated_nt ) # translate mutated protein and look for stop codons
                # the m_wt might have multiple bases, all of which need to be replaced by the m_mut, thus pass both
                indels_to_add.append([mutsite, wt_bases, list(m_mut), mutsite+len(m_wt)]) # start till index that makes 1st slice, mutation, index for start of 2nd slice                    
            else: # complement gene
                m_wt_comp = [complement[base] for base in list(m_wt)]
                m_wt_revcomp = m_wt_comp[::-1]
                wt_bases = gene_nt[-mutsite-len(m_wt)+1:-mutsite+1] # +1 because reverse indices start from -1, not -0
                #print wt_bases, m_wt_revcomp
                comp_m_mut = [complement[base] for base in list(m_mut)]
                revcomp_m_mut = comp_m_mut[::-1]
                #print "1st slice:", gene_nt[0:-(mutsite+len(m_wt)-1)], "\nmutation:", list(revcomp_m_mut), "\n2nd slice", gene_nt[-mutsite+1:]                    
                mutated_nt = gene_nt[0:-(mutsite+len(m_wt)-1)] + list(revcomp_m_mut) + gene_nt[-mutsite+1:]
                stop_codon_flag = findStopCodon( mutated_nt ) # translate mutated protein and look for stop codons
                indels_to_add.append([-(mutsite+len(m_wt)-1), wt_bases, revcomp_m_mut, -mutsite+1])

        else: # single base substitution
            if gstart < gend: # normal gene
                if gene_nt[mutsite] == m_wt: # reference base matches
                    gene_nt[mutsite] = m_mut # insert mutation at this position
                    stop_codon_flag = findStopCodon( gene_nt ) # translate mutated protein and look for stop codons
                    gene_nt[mutsite] = m_wt # revert the mutation back to wt after detecting the impact of it 
                    substitutions.append([mutsite, m_wt, m_mut])
                else:
                    print "Error! Reference base in gene does not match the one in mutation!\n\n"
            else: # complement gene. gene_nt is already reverse complement of genome seq.
                if gene_nt[-mutsite] == complement[m_wt]: # reverse complement match for ref nt
                    gene_nt[-mutsite] = complement[m_mut] # insert complement of mutation at this site
                    stop_codon_flag = findStopCodon( gene_nt ) # translate mutated protein and look for stop codons
                    gene_nt[-mutsite] = complement[m_wt] # revert the mutation back to wt after detecting the impact of it 
                    substitutions.append([-mutsite, complement[m_wt], complement[m_mut]])
                else:
                    print "Error! Reference base in gene does not match the one in mutation!\n\n"
        #if stop_codon_flag:
        #    print "\n", gname, m, "\n"
    # add indels. DId not add them before becaues they would mess up the reference position numbering, and we won't get correct base substitutions

    if stop_codon_flag == 1: # if a stop codon gets introduced in the middle of the protein, print all mutations for that protein
        print "\nmutlist for ", gname, ":", mutlist
        gene_allmut = addMutations( gname, gene_nt, substitutions, indels_to_add )
    return 1

# end of checkProtein ##########################

# start of checkImpact ##########################
def checkImpact( pos, gnames, m_wt, m_mut, ntseq, aaseq, gdict ):
    for g in gnames:
        if g not in ntseq: # if this gene name does not has sequence
            continue
        gene_nt = list(ntseq[g])
        gene_aa = aaseq[g]
        gstart, gend = gdict[g]
        mutsite = pos-gstart
        # print mutsite, mutsite/3
        if gstart < gend: # normal gene
            if gene_nt[mutsite] == m_wt: # reference base matches
                gene_nt[mutsite] = m_mut # insert mutation at this position
            else:
                print "Error! Reference base in gene does not match the one in mutation!\n\n"
        else: # complement gene. gene_nt is already reverse complement of genome seq.
            if gene_nt[-mutsite] == complement[m_wt]: # reverse complement match for ref nt
                gene_nt[-mutsite] = complement[m_mut] # insert complement of mutation at this site
            else:
                print "Error! Reference base in gene does not match the one in mutation!\n\n"

        trans_aa = translateGene( gene_nt )
        trans_aa = trans_aa.rstrip('*') # remove any stop codons at the end if there are any
        if trans_aa[1:] == gene_aa[1:]: # since 1st residue can be Met or something else, ignore that site
            return 0
        else: # protein has changed due to mutation: non-syn
            #print trans_aa, "\n", gene_aa
            return 1
# end of checkImpact ##########################

# start of checkMut ##########################
def checkMut( m, genes, ntseq, aaseq, gdict ): # check is a mutation is synonymous or non-syn
    m_wt, mpos, m_mut = m.split('_')
    # check if the mutation is in a coding region. If yes, check for syn/non-syn. If no, don't save the mut
    ipos = int(mpos)
    gflag = sflag = 0
    gname = []
    gflag, gname = checkRegion( ipos, genes )
    if gflag == 0: # mutation is in inter-genic region, do not save
        return 0
    else: # mutation is in genic region.
        sflag = checkImpact( ipos, gname, m_wt, m_mut, ntseq, aaseq, gdict )
        if sflag == 0: # synonymous mut, do not save
            return 1 # return syn mutations as well for now, just because with an indel they might not be syn anymore
        else: # non-synonymous, save
            return 1
# end of checkMut ##########################

# start of getCoverage ##########################
def getCoverage( strain, mutposlist, cov_thres ): # read mpileup file for this strain to get read depth at every position
    fh = open("/home/ag1349/Hassan/Zhangetal2013/NGS_analysis/"+strain+"/mapped/"+strain+".cov", 'r')
    coverage = {} # key is mutation pos, value is coverage at that pos
    for line in fh:
        line = line.rstrip('\n')
        content = line.split()
        gpos = int(content[0])
        poscov = int(content[1])
        # 1st index: genome pos, 3rd index: coverage
        if poscov >= cov_thres:
            if gpos in mutposlist: # genome position has a mutation
                coverage[gpos] = poscov
    fh.close()
    return coverage
# end of getCoverage ##########################    

# start of pruneMutlist ##########################
def prunedMutlist( genewise_mut, s_cov ):
    pruned = {} # filtered list of only those mutations that have coverage higher than a threshold
    for g in genewise_mut: # for each gene
        pruned[g] = []
        for m in genewise_mut[g]: # for each mutation in the gene
            m_wt, mpos, m_mut = m.split('_')
            ipos = int(mpos)
            if ipos in s_cov: # if the mutatios is in hash with high coverage sites, save it
                pruned[g].append( m )
    return pruned
# end of pruneMutlist ##########################

# start of getStrainSubs ######################
def getStrainSubs( gnames, mutkey, genes, mutations ):
    genemut = {gn : [] for gn in gnames} # key: gene, val: list of mutations in that gene
    mutpos_strain = [] # get all positions in this strain that have a mutation. This list is to get coverage for these pos.
    for i in range(len(mutkey)): # 0th index in mutkey was for "Strain", but that is not passed to this function
        if mutkey[i] == 1: # mutation is present in this strain
            m_wt, mpos, m_mut = mutations[i].split('_')
            mutpos_strain.append( int(mpos) )
            # get the genes in which this mutation exits
            mflag, genes_mut = checkRegion( int(mpos), genes )
            for g in genes_mut: # since this is genic mutation, there will be at least one gene with this mutation 
                if g in genemut: # only consider coding genes as per H37Rv gene and protein files from NCBI
                    genemut[g].append( mutations[i] )
    return [genemut, mutpos_strain] 
# end of getStrainSubs ######################

# start of getStrainIndels #####################
def getStrainIndels( strain, gnames, genes ):
    genemut = {gn : [] for gn in gnames} # key: gene, val: list of mutations in that gene
    fh = open("/home/ag1349/Hassan/Zhangetal2013/NGS_analysis/"+strain+"/mapped/"+strain+".indels", 'r')
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
def checkFrameshifts( gname, gene_nt, gene_aa, pindels ): # these indels are processed indels returned by checkProtein2 function
    mutatedseq = addMutations( gname, list(gene_nt), [], pindels )
    if findStopCodon(mutatedseq) == 0: # no stop codon
        print gname, "\t checking for frameshifts", pindels, gene_nt
        translated = translateGene(mutatedseq)
        print "\nOriginal\n", gene_aa, "\nTranslated\n", translated
        first_mismatch = 0 # first residue from the beginning of the protein that is mismatched
        last_mismatch_end = 0 # last residue (from end) that is mismatched between mutated and wt proteins
        # compare wt and mutated protein sequences from both ends to see if the frame was restored by the second mutation
        if translated == gene_aa:
            print "Mutated seq exactly same as wt. Restored?\n"
        else:
            orig = list(gene_aa)
            mut = list(translated)
            for i in range(1,len(orig)): # ignore 1st residue, as methionine is sometimes posttranslationally modified
                if orig[i] == mut[i]:
                    continue
                else:
                    first_mismatch = i
                    print "1st mismatch from L->R:", i
                    break 
            # check for 1st mismatch from the reverse side:
            rev_orig = orig[::-1]
            rev_mut = mut[::-1]
            for i in range(max(len(orig), len(mut))):
                if rev_orig[i] == rev_mut[i]:
                    continue
                else:
                    last_mismatch_end = i
                    print "1st mismatch from R->L (from end):", i
                    break
        if first_mismatch > 0 and last_mismatch_end > 0:
            print "==Frameshift!==\n\n"

# end of checkFrameshifts ######################

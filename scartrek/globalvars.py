#!/usr/bin/env python

# start of extractSeq ##########################
def extractSeq( fname, typeflag ): # extract gene or protein seq from respective file handles. typeflag = 0 for ntseq, typeflag = 1 for protein seq
    fh = open(fname, 'r')
    seqdict = {}
    gposdict = {}
    genes = []
    seq = ""
    for line in fh:
        if ">" in line: # header line
            if seq != "": # there is some sequence to save
                seqdict[gname] = seq.replace('\n','')
            content = line.split()
            info = {'gene' : '', 'locus_tag' : '', 'location' : ''}
            for i in range(1, len(content)):
                key, val = content[i].split('=')
                key = key.lstrip('[')
                val = val.rstrip(']')
                if key in info: # H37Rv_proteins has a lot more fields, only get these three
                    info[key] = val
                else:
                    continue
            if info['gene'] != '': # there is a gene name
                gname = info['gene']
            else:
                gname = info['locus_tag']
            seq = "" # reinitialize gene seq
            if typeflag == 0: # only gene seq file has location info
                locstr = info['location']
                if locstr.startswith("complement"): # info[2] is location string [location=12468..13016] or [location=complement(13133..13558)]
                    locstr = locstr.lstrip("complement(")
                    locstr = locstr.rstrip(")]")
                    #gname = gname+"_c" # not all complement genes end with a 'c'. Thus add this distinguishing feature to the gene name
                elif locstr.startswith("order"): # only one gene has this. [locus_tag=Rv3216] [location=order(3593369..3593437,3593439..3593852)]
                    # enter this gene manually, encompassing the entire region covered by this gene
                    locstr = "3593369..3593852"
                start, end = locstr.split('..')
                start = start.lstrip('<') # some gene start positions have location mentioned as [location=<2550340..2551326]. Remove the "<" sign.
                end = end.lstrip('>') # some gene end positions have location mentioned as 817531..>817866. Remove the '>' sign.
                istart = int(start)
                iend = int(end)
                gposdict[gname] = [istart, iend]
                if istart < iend:
                    genes.append([ istart, iend, gname ])
                else:
                    genes.append([ iend, istart, gname ])

                #print(gname, start, end)
        else:
            seq += line
    fh.close()
    # add the last seq to the dict        # bugfix 1
    seqdict[gname] = seq.replace('\n','') # bugfix 1
    if typeflag == 0: # ntseq
        return seqdict, gposdict, genes         
    elif typeflag == 1:
        return seqdict
# end of extractSeq ##########################




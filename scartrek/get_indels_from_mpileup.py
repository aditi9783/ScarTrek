#!/usr/bin/python

# start of findIndels ########################
def findIndels( fpath, sd, COVTHRES ):

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
                print("Incorrect mpileup:", fpath, sd, line)
                break
        else: # the line has less than 4 elements, incorrect line.
            print("Incorrect mpileup:", sd, fpath, line)
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

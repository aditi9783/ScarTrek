#!/usr/bin/env python

import sys, os, subprocess
from multiprocessing import Pool

def qcmap( seqinfo ): # Input: tuple of [seqdir, sd, fq_1, fq_2]
    trim_path = # add trimmomatic path, or ability to use installed module 
    adaptor_path = # same as above for trimmomatic adaptors
    seqdir, sd, fq_1, fq_2 = seqinfo
    # create a new dir for trimmomatic output (qc)
    qcdir = seqdir+"/qc"
    os.makedirs(qcdir)
    os.chdir(qcdir)
    trim_cmd = "java -jar "+trim_path+" PE -phred33 "+fq_1+" "+fq_2+" s1_pe.fq.gz s1_se.fq.gz s2_pe.fq.gz s2_se.fq.gz ILLUMINACLIP:"+adaptor_path+":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:20 &> qc.run"
    trim_out= subprocess.Popen(trim_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    #print('\nTrimmomatic command:\n',trim_cmd,'\n')
    print('Trimmomatic Output:\n',trim_out.communicate(),'\n')
    # run FastQC to detect qual control
    qcfiles = ["s1_pe.fq.gz", "s1_se.fq.gz", "s2_pe.fq.gz", "s2_se.fq.gz"]
    print("\nRunning FastQC")
    for f in qcfiles:
        subprocess.call([fastqcpath, qcdir+"/"+f])
    print("Completed QC analysis. Files are in ",qcdir)

    # map the reads to the reference H37Rv using bowtie2 (reference has already been indexed)
    # assume that bowtie2 has been loaded in the environment
    mapdir = seqdir+"/mapped"
    os.makedirs(mapdir)
    os.chdir(mapdir)
    bowtie_cmd = "bowtie2 -p 1 -x "+refpath+" -1 "+qcdir+"/s1_pe.fq.gz -2 "+qcdir+"/s2_pe.fq.gz -U "+qcdir+"/s1_se.fq.gz,"+qcdir+"/s2_se.fq.gz -S "+mapdir+"/aln.sam &> map.run"
    bowtie_out = subprocess.Popen(bowtie_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print("Bowtie output:", bowtie_out.communicate())
            
    # use samtools (assume pre-loaded) to convert sam to bam file, and then to sort the bam file
    os.system("samtools view -bS aln.sam > aln.bam")
    os.system("samtools sort aln.bam aln.sorted")
    os.system("samtools mpileup -f "+refgenomepath+" aln.sorted.bam > aln.sorted.bam.mpileup")
    # call variants: assume bcftools is loaded in the environment

    # the following commands work for samtools/1.2 and bcftools/1.2, both of which should be loaded in env
    os.system("samtools mpileup -go "+sd+".bcf -f "+refgenomepath+" aln.sorted.bam")
    os.system("bcftools call -vmO z -o "+sd+".vcf.gz "+sd+".bcf")
    os.system("rm "+sd+".bcf") # delete bcf file as mpileup output is already stored

    # change sorted bam file name to include strain name (subdir name) as well
    os.system("mv aln.sorted.bam "+sd+".aln.sorted.bam")

    # create index for reading bam file in Tablet
    os.system("samtools index "+sd+".aln.sorted.bam")

    # delete sam file, unsorted bam file, raw fastq files, and quality controlled fastq files
    os.system("rm aln.sam")
    os.system("rm aln.bam")
    os.system("gzip aln.sorted.bam.mpileup")
    os.chdir(qcdir)
    os.system("rm *_*.fq.gz")
    os.system("rm "+fq_1) # fq_1 has full path to the file
    os.system("rm "+fq_2) # fq_2 has full path to the file

if __name__ == '__main__':
    seqlist = []
    seqpath = sys.argv[1]
    for typedir in os.walk(seqpath): # seqd is a tuple of [subdir, sub-subdir, and list of files]
        seqids = typedir[1] # list of sub-dir names after seq ids
        for sd in seqids:
            #if len(seqlist) == 2: # for testing
            #    break   
            fq_1 = ""
            fq_2 = ""
            seqdir = typedir[0]+sd
            print("\n",seqdir)
            allf = os.listdir(seqdir)
            for seqf in allf: # get paired end files
                if "_1.fastq.gz" in seqf:
                    fq_1 = seqdir+"/"+seqf
                elif "_2.fastq.gz" in seqf:
                    fq_2 = seqdir+"/"+seqf
            seqlist.append([seqdir, sd, fq_1, fq_2])
    #print(seqlist)
    pool = Pool(15) # parallel processors
    pool.map(qcmap, seqlist)
    pool.close()
    pool.join()

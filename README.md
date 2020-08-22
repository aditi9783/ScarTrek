# ScarTrek
Python package for identifying gene disrupting and restoring indels in whole-genome sequencing data of Mycobacterium tuberculosis samples


-----------------------------------------------------------------------------------------
### Description

"Insertion" or "deletions" (affectionly called "indels") of nucleotides are two types of genetic changes besides 
single nucleotide polymorphims (SNPs) that can be detected by analyzing sequencing data. ScarTrek is an application
written in python to detect indels in whole-genome sequencing data mapped to a reference genome, and if gene information
for the referece organism is provided, then determine the genes that have the indels and if the indels influence the 
translated product of the gene. This analyses is useful to detect gene inactivations due to indels.

### Installation

[//]: # (##### Using pip)
[//]: # (`$ pip install scartrek`)

##### From GitHub

ScarTrek project can be downloaded directly from GitHub and python scripts can be run directly from the command line as shown in usage.

### Usage

[//]: # (If installed using pip, ScarTrek can be run as:)

[//]: # (`$ find-scars [-h] -i INPUT [-m MAPRATE] [-c COVTHRES] [-g GENESEQ]
                     [-p PROTSEQ]`)

If downloaded from GitHub, ScarTrek can be run in the directory ScarTrek/scartrek/ as:

`$ python find_scars.py [-h] -i INPUT [-m MAPRATE] [-c COVTHRES] [-q MAPQ]
                     [-s SBALANCE] [-g GENESEQ] [-p PROTSEQ]`

where the arguments are:
```
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input directory that has mpileup files for each
                        sample. See tests/test1 for example (default: None)
  -m MAPRATE, --maprate MAPRATE
                        Minimum read mapping rate required to consider a
                        sample (default: 20.0)
  -c COVTHRES, --covthres COVTHRES
                        Minimum read coverage required at a position to detect
                        an indel (default: 20)
  -q MAPQ, --mapq MAPQ  Minimum average mapping quality required at a position
                        to detect an indel. This requires Samtools as well as
                        .bam and .bai files. Set this parameter to 0 if do not
                        want to use this filter. (default: 10)
  -s SBALANCE, --sbalance SBALANCE
                        Minimum forward/reverse strand balance required at a
                        position to detect an indel (default: 0.05)
  -g GENESEQ, --geneseq GENESEQ
                        Gene sequences in the reference genome, default
                        reference: M. tuberculosis (default:
                        ../reference/H37Rv_genes.txt)
  -p PROTSEQ, --protseq PROTSEQ
                        Protein sequences for the reference organism, default:
                        M. tuberculosis (default:
                        ../reference/H37Rv_proteins_from_genbank.txt)
```

##### Example usage:
`$ python find_scars.py -i ../example/ -q 0`

The above example runs in <30 seconds. The MAPQ parameter is set to 0 (`-q 0`) because .bam and .bai files for the example are not included due to large file size.

ScarTrek has been tested on CentOS Linux 7 (Core) and Mac OS 10.13.5.


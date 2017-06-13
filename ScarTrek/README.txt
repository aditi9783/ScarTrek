
ScarTrek - Python package for identifying genetic changes in whole-genome sequencing data
-----------------------------------------------------------------------------------------

"Insertion" or "deletions" (affectionly called "indels") of nucleotides are two types of genetic changes besides 
single nucleotide polymorphims (SNPs) that can be detected by analyzing sequencing data. ScarTrek is an application
written in python to detect indels in whole-genome sequencing data mapped to a reference genome, and if gene information
for the referece organism is provided, then determine the genes that have the indels and if the indels influence the 
translated product of the gene. This analyses is useful to detect gene inactivations due to indels.

ScarTrek can be run in parallel if there are a large number of samples.

Additional scripts are provided to do the following upstream analyses: quality control of raw reads in fastq format and mapping the filtered reads to the reference genome. Although, the user will need to install (or load in their environment) third party bioinformatics tools to process fastq files and map the reads (see DEPENDENCIES.txt). High Performance Computing Centers at major universities already have these resources installed.

For any questions, comments, or suggestions, please contact Aditi Gupta at ag1349@njms.rutgers.edu or aditi9783@gmail.com.

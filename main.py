#! /usr/bin/python3.10
import subprocess as sb
from Bio import Entrez
import multiprocessing
import os
from Bio import SeqIO


BWAI = "bwa index %s"
BWA = "bwa mem -t %s %s %s"
SMS = "samtools sort -@ %s -O BAM -o %s"

def download(reference):
    cwd = os.getcwd()
    Entrez.email = 'A.N.Other@example.com'
    #taking data
    net_handle = Entrez.efetch(db='nucleotide', id = reference, rettype = 'fasta', retmode = 'text')
    data = net_handle.read()
    net_handle.close()

    ref = os.path.join(cwd, "my_file.fasta")

    out_handle = open(ref, "w")
    out_handle.write(data)
    out_handle.close()



    return(ref)

def bwa(reference, sraFile):
    #BWAI = "bwa index %s"
    #BWA = "bwa mem -t %s %s %s"

    bi = BWAI % reference
    bwa_index = sb.Popen(bi, shell=True)
    bwa_index.communicate()
    for sra in sraFile:
        if sra[1]:
            sraFileName = sra[0] + "_1.fastq " + sra[0] + "_2.fastq"
        else:
            sraFileName = sra[0] + ".fastq "
        bwa_mem = BWA % ("10", reference, sraFileName)
        samtools_sort = SMS % (10, sra[0] + ".sorted.bam")
        bwa_mem_out = sb.Popen(bwa_mem, stdout=sb.PIPE , shell=True)
        samtools_sort_out = sb.Popen(samtools_sort, stdin=bwa_mem_out.stdout, shell=True)
        samtools_sort_out.communicate()
        print("done")

def sra(sra_numbers):
    cwd = os.getcwd()
    pe_sra = []
    # samples
    #sra_numbers = [species]
    # this will download the .sra files to ~/ncbi/public/sra/ (will create directory if not present)
    # for sra_id in sra_numbers:
    #     #print ("Currently downloading: " + sra_id)
    #     prefetch = 'prefetch ' + sra_id + "-O " + cwd
    #     #print ("The command used was: " + prefetch)
    #     sb.call(prefetch, shell=True)

    # this will extract the .sra files from above into a folder named 'fastq'
    for sra_id in sra_numbers:
        #print ("Generating fastq for: " + sra_id)
        fastq_dump = "fasterq-dump " + sra_id + " -O " + cwd
        #print ("The command used was: " + fastq_dump)
        sb.call(fastq_dump, shell=True)

        if os.path.exists(os.path.join(cwd, sra_id + "_1.fastq")) :
            pe = True
        else:
            pe = False
        pe_sra.append((sra_id, pe))

    return(pe_sra)

def main():
    pe_sra = sra(["SRR21936789","SRR21936788"])
    reference = download("NZ_CP007265.1")
    bwa(reference, pe_sra)

if __name__ == '__main__':
    main()
#! /usr/bin/python3.10
import subprocess as sb
from Bio import Entrez
from Bio import SeqIO


BWAI = "bwa index %s"
BWA = "bwa mem -t %s %s %s"

def download(reference):
    Entrez.email = 'A.N.Other@example.com'
    #taking data
    handle = Entrez.efetch(db='nucleotide', id = reference, rettype = 'fasta', retmode = 'text')
    record = SeqIO.read(handle, 'fasta')
    #test
    #taking results
    seq = record.seq
    seqq = str(seq)
    #length of the result
    length = len(seqq)
    str_length = str(length)

    #output in a file
    reference= open("/tmp/rpova.txt", "w")
    reference.write(seqq+'\n'+str_length)
    reference.close()
    return("rpova.txt")

def bwa(reference, sraFile, pe):
    #BWAI = "bwa index %s"
    #BWA = "bwa mem -t %s %s %s"

    bi = BWAI % reference
    bwa_index = sb.Popen(bi, shell=True)
    bwa_index.communicate()
    if pe:
        sraFile = sraFile + "_1.fastq " + sraFile + "_2.fastq"
    else:
        sraFile = sraFile + ".fastq "
    bwa_mem = BWA % ("10", reference, sraFile)
    bwa_mem_out = sb.Popen(bwa_mem, shell=True)
    bwa_mem_out.communicate()
    print("done")
    # a = subprocess.getstatusoutput('mkdir index')
    # b = subprocess.getstatusoutput('mv ' + reference + ' index/')
    # c = subprocess.getstatusoutput('bwa index index/'+ reference)
    # d = subprocess.getstatusoutput('bwa mem -t 8 index/GCF_000013925.1_ASM1392v2_genomic.fna.gz P7741_R1.fastq.gz P7741_R2.fastq.gz > output.sam')

def sra(species):
    # samples
    sra_numbers = [species]

    # this will download the .sra files to ~/ncbi/public/sra/ (will create directory if not present)
    for sra_id in sra_numbers:
        #print ("Currently downloading: " + sra_id)
        prefetch = 'prefetch ' + sra_id
        #print ("The command used was: " + prefetch)
        sb.call(prefetch, shell=True)

    # this will extract the .sra files from above into a folder named 'fastq'
    for sra_id in sra_numbers:
        #print ("Generating fastq for: " + sra_id)
        fastq_dump = "fasterq-dump " + sra_id
        #print ("The command used was: " + fastq_dump)
        sb.call(fastq_dump, shell=True)
    # if ....:
    pe = True
    # else:
    #     pe = False
    return(pe, species)

def main():
    pe, sraFile = sra("SRR21936789")
    reference = download("NZ_CP007265.1")
    bwa(reference, sraFile, pe)

if __name__ == '__main__':
    main()
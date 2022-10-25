#! /usr/bin/python3.10
import subprocess as sb
from Bio import Entrez
from Bio import SeqIO


BWAI = "bwa index %s"
BWA = "bwa mem -t %s %s %s"

def download():
    Entrez.email = 'A.N.Other@example.com'
    #taking data
    handle = Entrez.efetch(db='nucleotide', id = 'NZ_CP007265.1', rettype = 'fasta', retmode = 'text')
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

def bwa(reference):
    #BWAI = "bwa index %s"
    #BWA = "bwa mem -t %s %s %s"

    bi = BWAI % reference
    bwa_index = sb.Popen(bi, shell=True, cwd="/tmp")
    bwa_index.communicate()
    # a = subprocess.getstatusoutput('mkdir index')
    # b = subprocess.getstatusoutput('mv ' + reference + ' index/')
    # c = subprocess.getstatusoutput('bwa index index/'+ reference)
    # d = subprocess.getstatusoutput('bwa mem -t 8 index/GCF_000013925.1_ASM1392v2_genomic.fna.gz P7741_R1.fastq.gz P7741_R2.fastq.gz > output.sam')

def main():
    reference = download()
    bwa(reference)

if __name__ == '__main__':
    main()
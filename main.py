#! /usr/bin/python3.10
import subprocess
from Bio import Entrez
from Bio import SeqIO


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
    reference= open("rpova.txt", "w")
    reference.write(seqq+'\n'+str_length)
    reference.close()
    return("rpova.txt")

def bwa(reference):

    a = subprocess.getstatusoutput('mkdir index')
    b = subprocess.getstatusoutput('mv ' + reference + ' index/')
    c = subprocess.getstatusoutput('bwa index index/'+ reference)
    d = subprocess.getstatusoutput('bwa mem -t 8 index/GCF_000013925.1_ASM1392v2_genomic.fna.gz P7741_R1.fastq.gz P7741_R2.fastq.gz > output.sam')

def main():
    reference = download()
    bwa(reference)

if __name__ == '__main__':
    main()
#! /usr/bin/python3.10
import subprocess as sb
from Bio import Entrez
import os
from Bio import SeqIO
import argparse
import tempfile



BWAI = "bwa index %s"
BWA = "bwa mem -t %s %s %s"
SMS = "samtools sort -@ %s -O BAM -o %s"
FREEBAYES = "freebayes -f %s %s"
SAMTOOLSRG = 'samtools  addreplacerg -r "@RG\tID:%s\tSM:%s" -o %s %s'

def options():
    parser = argparse.ArgumentParser(prog='quantify and detect pathogens', usage='%(prog)s [options]')
    parser.add_argument('-t','--threads', type=int, default=1, choices=range(1, os.cpu_count()))
    parser.add_argument('-r', '--reference', type=str, required=True, help="NCBI ID for the genome to use as reference; example is NZ_CP007265.1")
    parser.add_argument('-em', '--email', type=str, required=True, help="email address to use in NCBI downloads")
    parser.add_argument('-s', '--species', type=str, help="species to use in NCBI downloads")
    args = parser.parse_args()
    return (args)

def vcfScan(vcf_file):

    nfields = 9
    dict_depth = {}
    index_sample = []
    dinfo = {"AD":[],"RO":[],"AO":[]}
    list_to_remove = []
    with open(vcf_file, "r") as infile:
        for line in infile:
            if not line.startswith("##"):
                el = line.strip().split("\t")
                if line.startswith("#"):
                    while nfields < len(el):
                        a = el[nfields]
                        dict_depth[a] = []
                        index_sample.append((a,nfields))
                        nfields += 1
                if not line.startswith("#"):
                    for sample in index_sample:
                        pos = el[0] + "_" + el[1]
                        info = el[sample[1]].split(":")
                        if len(info) == 8:
                            dict_pos = {}
                            dict_pos[pos] = [el[3], el[4], info[2], info[3], info[5]]
                            dict_depth[sample[0]] = dict_depth[sample[0]] + [dict_pos]
                        else:
                            list_to_remove.append(pos)

        print("done")

                    #ntotal = len(el) - nfix





def vcfFreebayes(genoma, bams):
    vcf = tempfile.NamedTemporaryFile(dir= os.getcwd(), suffix=".vcf", delete=False, mode="w")
    free_comm = FREEBAYES % (genoma, " ".join(bams))
    print(free_comm)
    free_comm_out = sb.Popen(free_comm, stdout=vcf, shell=True)
    free_comm_out.communicate()
    print("done")
    return (vcf.name)

def download(args):
    cwd = os.getcwd()
    Entrez.email = args.email
    #taking data
    net_handle = Entrez.efetch(db='nucleotide', id = args.reference, rettype = 'fasta', retmode = 'text')
    data = net_handle.read()
    net_handle.close()

    ref = os.path.join(cwd, "my_file.fasta")

    out_handle = open(ref, "w")
    out_handle.write(data)
    out_handle.close()
    return(ref)

def bwa(reference, sraFile, args):
    #BWAI = "bwa index %s"
    #BWA = "bwa mem -t %s %s %s"
    bi = BWAI % reference
    bwa_index = sb.Popen(bi, shell=True)
    bwa_index.communicate()
    bams = []
    for sra in sraFile:
        if sra[1]:
            sraFileName = sra[0] + "_1.fastq " + sra[0] + "_2.fastq"
        else:
            sraFileName = sra[0] + ".fastq "
        bwa_mem = BWA % (str(args.threads), reference, sraFileName)
        print(bwa_mem)
        samout = sra[0] + ".sorted.bam"
        samoutrg = sra[0] + ".rg.sorted.bam"
        samtools_sort = SMS % (str(args.threads), samout)
        print(samtools_sort)
        bwa_mem_out = sb.Popen(bwa_mem, stdout=sb.PIPE , shell=True)
        samtools_sort_out = sb.Popen(samtools_sort, stdin=bwa_mem_out.stdout, shell=True)
        samtools_sort_out.communicate()
        samtoosls_rg_com = SAMTOOLSRG % (sra[0], sra[0], samoutrg, samout)
        print(samtoosls_rg_com)
        samtoosls_rg_com_out= sb.Popen(samtoosls_rg_com, shell=True)
        samtoosls_rg_com_out.communicate()
        print("done")
        bams.append(samoutrg)
    return(bams)

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
    #args = options()
    #pe_sra = sra(["SRR21936789","SRR21936788"])
    #reference = download(args)
    #bams = bwa(reference, pe_sra, args)
    #vcf_file=vcfFreebayes(reference, bams)
    vcfScan("/home/lfaino/Bioinformatics-in-plant-pathology/tmpv6c3kg3x.vcf")

if __name__ == '__main__':
    main()
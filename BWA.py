import subprocess

a = subprocess.getstatusoutput('mkdir index')
b = subprocess.getstatusoutput('mv GCF_000013925.1_ASM1392v2_genomic.fna.gz index/')
c = subprocess.getstatusoutput('bwa index index/GCF_000013925.1_ASM1392v2_genomic.fna.gz')
d = subprocess.getstatusoutput('bwa mem -t 8 index/GCF_000013925.1_ASM1392v2_genomic.fna.gz P7741_R1.fastq.gz P7741_R2.fastq.gz > output.sam')
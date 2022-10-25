from Bio import Entrez
from Bio import SeqIO

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
f= open("rpova.txt", "w")
f.write(seqq+'\n'+str_length)
f.close()


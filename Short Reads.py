import subprocess

# samples
sra_numbers = ["SRR21936789"]

# this will download the .sra files to ~/ncbi/public/sra/ (will create directory if not present)
for sra_id in sra_numbers:
    #print ("Currently downloading: " + sra_id)
    prefetch = 'prefetch ' + sra_id
    #print ("The command used was: " + prefetch)
    subprocess.call(prefetch, shell=True)

# this will extract the .sra files from above into a folder named 'fastq'
for sra_id in sra_numbers:
    #print ("Generating fastq for: " + sra_id)
    fastq_dump = "fasterq-dump " + sra_id
    #print ("The command used was: " + fastq_dump)
    subprocess.call(fastq_dump, shell=True)

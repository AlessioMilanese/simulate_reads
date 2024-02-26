# code for a python script to simulate reads (fastq files) from genomes.
# The input:
# - repo with the genomes fasta files
# - file with the abundances (each line is a genome that we want to use. It can also be only one genome. It's a tab separated with two columns, the second column is the abundance like 8)
# - number of reads to create

# As a result it produces one fastq file. The headers are the genome fasta name and a counter

# load inputs
import sys
import random

# repo
repo = sys.argv[1] # like "/path/to/repo/with/genomes/"
if repo[-1] != "/":
    repo += "/"
# abundances
abundances_file = sys.argv[2]
# like "/path/to/abundances.txt"
# the file has two columns, the first is the genome name and the second is the abundance
# like:
# genome1.fastq 8
# genome2.fastq 4

# number of reads
n_reads = int(sys.argv[3]) # like 100000
# output
output = sys.argv[4]

# length of the reads
length_reads = 150
# number of Ns to separate the headers of the genomes when there are multiple headers
n_Ns = 100
# max number of errors per read
max_errors = 4

# -------------------------------------------------------------------------------------------------
# load the abundances
abundances = {}
total_abundance = 0
o = open(abundances_file)
for line in o:
    line = line.strip()
    genome, abundance = line.split("\t")
    abundances[genome] = int(abundance)
    total_abundance += int(abundance)
o.close()

# -------------------------------------------------------------------------------------------------
# for each genome in the abundances, we check if the file exists and we
# check the length of the genome
import os
import subprocess
genomes_length = {}
tot_genome_lengths = 0
for genome in abundances:
    # check that the file exists
    if not os.path.exists(repo + genome):
        sys.stderr.write("ERROR: file not found " + repo + genome)
        sys.exit(1)
    # check the length of the genome
    command = "grep -v '>' " + repo + genome + " | wc -c"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    output_t = process.communicate()[0]
    length = int(output_t.strip())
    # add to the dictionary
    genomes_length[genome] = length
    # and update the total length
    tot_genome_lengths += length


# -------------------------------------------------------------------------------------------------
# we calculate how many reads to generate per each genome
# note that we want to create a number of reads that is proportional to the abundance
# but taking into account the length of the genome

# number of reads per genome (before normalisation)
reads_per_genome = {}
for genome in abundances:
    reads_per_genome[genome] = abundances[genome] * n_reads / total_abundance

# we normalise by the length of the genome
reads_per_genome_normalised = {}
for genome in reads_per_genome:
    reads_per_genome_normalised[genome] = reads_per_genome[genome] * genomes_length[genome] / tot_genome_lengths
# TOT: 50 reads
# rel ab: 5, 5 and 15
# reads_per_genome: 10, 10 and 30
# genome lengths: 100, 50 and 50
# reads_per_genome_normalised 1: 10*100/200=5, 10*50/200=2.5 and 30*50/200=7.5

# and now we normalise again to make sure that the total number of reads is the same as n_reads
# we calculate the total number of reads
tot_reads = sum([reads_per_genome_normalised[genome] for genome in reads_per_genome_normalised])
# we calculate the factor to normalise
factor = n_reads / tot_reads
# we normalise
for genome in reads_per_genome_normalised:
    reads_per_genome_normalised[genome] = int(reads_per_genome_normalised[genome] * factor)

# -------------------------------------------------------------------------------------------------
# function to create the reads
def create_reads(genome, sequence, n_reads, length_reads, output):
    # open the output, we append or create the file if it doesnt exist
    o = open(output, "a")
    # create the reads
    for rrr in range(n_reads):
        start = random.randint(0, len(sequence) - length_reads)
        read = sequence[start:start + length_reads]
        # for the read we can insert some errors
        n_errors = random.randint(0, max_errors)
        positions = random.sample(range(length_reads), n_errors)
        for position in positions:
            read = read[:position] + random.choice("ACGT") + read[position + 1:]
        # for the quality we choose between 34 and 40
        values = ["C", "D", "E", "F", "G", "H", "I"]
        quality = "".join([random.choice(values) for i in range(length_reads)])
        # if there are Ns in the read, we set the quality of the position where
        # there is a N to 7
        quality = list(quality)
        for i in range(length_reads):
            if read[i] == "N":
                quality[i] = "("
        quality = "".join(quality)

        # we write the four lines of the fastq files
        o.write("@" + genome + "_" + str(rrr) + "\n")
        o.write(read + "\n")
        o.write("+\n")
        o.write(quality + "\n")
    o.close()


# -------------------------------------------------------------------------------------------------
# load the genomes and generate reads
# note that there might be multiple fasta headers for the same genome, hence we separate them
# with Ns
import os
for genome in abundances:
    # open the file
    o = open(repo + genome)
    # read the file
    sequence = ""
    # skip first line
    o.readline()
    for line in o:
        line = line.strip()
        if line[0] != ">":
            sequence += line
        else:
            sequence += "N" * n_Ns
    o.close()
    create_reads(genome, sequence, reads_per_genome_normalised[genome], length_reads, output)





#!/usr/bin/env python3

import argparse
import re

""" Given two fasta reference sequences and one fastq read file, the script will identify reads containing k-mers from both reference files (i.e. chimeric reads).
    If the reference sequences are viral genomes with low nucleotide similarity, the existence of chimeric reads may indicate that recombination has occurred. """




def get_args():
    """make an object to store the command line arguments, with appropriate --help info"""

    parser = argparse.ArgumentParser(
        description='Detect k-mers of two fasta sequences in each read of a fastq file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-fq', '--fastq', help='a fastq file contaning sequence reads', type=str)
    parser.add_argument('-fa1', '--fasta_1', help='a fasta file containing first query sequence', type=str)
    parser.add_argument('-fa2', '--fasta_2', help='a fasta file containing second query sequence', type=str)
    parser.add_argument('-k', '--kmer_length', help='non-zero integer specifying the length of k', type=int, default=22)

    args = parser.parse_args()
    return args

# ---------------------------------------------------------------------------------------------------------------------
# Create a namespace object containing the neccessary arguments and the user inputs relating to each of them
args = get_args()
# ---------------------------------------------------------------------------------------------------------------------




def check_fastq_list(x, y):
    """Check if two lists are the same length"""
    if len(x) == len(y):
        return True
    else:
        return False

def get_fq_list():
    """Make a list with length=2 from a fastq file. The first element is a list of fq headers and the second element is the fq_seqs list (fq_list[1])"""

    with open(args.fastq) as fq:
        fq_split = fq.read().splitlines()

    fq_heads = fq_split[0::4]
    fq_seqs = fq_split[1::4]
    fq_list = (fq_heads, fq_seqs)
    
    if check_fastq_list(fq_list[0], fq_list[1]):
        return (fq_list)
    else:
        raise Exception("Error: fastq file is corrupt. number of headers != number of sequences")  

# ---------------------------------------------------------------------------------------------------------------------
# Create two lists. One containing the fastq headers, and the other containing the fastq reads
fq_list = get_fq_list()
# ---------------------------------------------------------------------------------------------------------------------




def get_str(x):
    """Extract a sequence from the fasta file. The File can be text-wrapped or not, but must contain only one fasta sequence"""

    fa_heads = []
    with open(x) as fa:
        fa_list = fa.read().splitlines() 

        for i in fa_list:
            if i.startswith('>'):
                fa_heads.append(i)
                head_pos = fa_list.index(i)
        
        if len(fa_heads) > 1:
            raise Exception("Error: input is a multi-fasta. This script only supports a single sequence per fasta file")
        elif len(fa_heads) == 0:
            raise Exception("Error: input fasta file is not formatted correctly. No header")
        elif head_pos != 0:
            raise Exception("Error: input fasta file is not formatted correctly. First line should be a header")
        else:
            fa1_str = "".join(fa_list[1::1])
        
        return(fa1_str)
    
def reverse_complement(seq): 
    """return the reverse-complement of a DNA or RNA sequence"""
    
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'U':'A', 'R':'Y', 'Y':'R', 'S':'S', 'W':'W', 'K':'M', 'M':'K', 'B':'B', 'V':'V', 'H':'H', 'N':'N'} 
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases

# ---------------------------------------------------------------------------------------------------------------------
# Create one string from the fasta sequence, and another from its reverse complement 
fa1_str = get_str(args.fasta_1)
rev_fa1_str = reverse_complement(fa1_str)
fa2_str = get_str(args.fasta_2)
rev_fa2_str = reverse_complement(fa2_str)
# ---------------------------------------------------------------------------------------------------------------------




def make_combined_kmer_set(seq1, seq2, k):
    """Make a combined kmer set from two strings"""

    kmer_list = []
    seq1_window = len(seq1) - k + 1
    seq2_window = len(seq2) - k + 1

    for i in range(seq1_window):
        kmer_list.append(seq1[i:i+k])
    for i in range(seq2_window):
        kmer_list.append(seq2[i:i+k])    
    return set(kmer_list)

# ---------------------------------------------------------------------------------------------------------------------
# Create a k-mer set consisting of all k-mers present in the fasta sequence and its reverse-complement
fa1_kmers = make_combined_kmer_set(fa1_str, rev_fa1_str, args.kmer_length)
fa2_kmers = make_combined_kmer_set(fa2_str, rev_fa2_str, args.kmer_length)
# ---------------------------------------------------------------------------------------------------------------------




def make_kmer_set(seq, k):
    """Make a set of unique kmers from a string"""
    kmer_list = []
    window = len(seq) - k + 1

    for i in range(window):
        kmer_list.append(seq[i:i+k])
    return set(kmer_list)

def kmer_search(fq_list, fa_kmer_set, k):
    """For each string in a list of strings, make a set of k-mers and check if any k-mers in that set match any of the k-mers in another k-mer set"""

    matches = []

    for i in fq_list:
        found=set.intersection(make_kmer_set(i, k), fa_kmer_set)  
        matches.append(len(found))
    return matches    

# ---------------------------------------------------------------------------------------------------------------------
# Create a list of integers, with each integer representing one read in the fastq file, and indicating the number of fasta k-mers detected in that read
fa1_counts = kmer_search(fq_list[1], fa1_kmers, args.kmer_length)
fa2_counts = kmer_search(fq_list[1], fa2_kmers, args.kmer_length)
# ---------------------------------------------------------------------------------------------------------------------



# ---------------------------------------------------------------------------------------------------------------------
# write a comma-separated file containing one row for each chimeric read, including read identifier, sequence, fasta1 k-mer counts, and fasta2 k-mer counts.
if re.search('^.+\.fq$', args.fastq):
    outname = re.search('(^.+)\.fq$', args.fastq).group(1)
elif re.search('^.+\.fastq$', args.fastq):
    outname = re.search('(^.+)\.fastq$', args.fastq).group(1)
else:
    outname = args.fastq

output = open(outname + '_output.csv', 'w+')
output.write('Read_Identifier' + ',' + 'Read_Sequence' + ',' + 'k_mers_' + args.fasta_1 + ',' + 'k_mers_' + args.fasta_2 + '\n')

for i, value in enumerate(fa1_counts):
    if fa1_counts[i] > 0 and fa2_counts[i] > 0:
        output.write(fq_list[0][i] + ',' + fq_list[1][i] + ',' + str(value) + ',' + str(fa2_counts[i]) + '\n')

output.close()
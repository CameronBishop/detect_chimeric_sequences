# detect_chimeric_sequences
An alignment-free way to find chimeric reads in a fastq file.


usage: detect_chimeric_sequences.py [-h] [-fq FASTQ] [-fa1 FASTA_1] [-fa2 FASTA_2] [-k KMER_LENGTH]

Detect k-mers of two fasta sequences in each read of a fastq file

arguments:
  -h, --help                                     show this help message and exit
  -fq FASTQ, --fastq FASTQ                       a fastq file contaning sequence reads (default: None)
  -fa1 FASTA_1, --fasta_1 FASTA_1                a fasta file containing first query sequence (default: None)
  -fa2 FASTA_2, --fasta_2 FASTA_2                a fasta file containing second query sequence (default: None)
  -k KMER_LENGTH, --kmer_length KMER_LENGTH      non-zero integer specifying the length of k (default: 22)

#!/usr/bin/env python 

from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import itertools
import argparse
import datetime

def get_arguments():
    parser = argparse.ArgumentParser(description='Count the frecuency of a k-mer set  given its size (k) along a genome or a set of genomes',
			formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input', type=str,help='Path to input fasta file')
    parser.add_argument('output', type=str,help='Path to put file/folder output')
    parser.add_argument('-k', type=int, default=1, help='size of k-mer')
    
    args = parser.parse_args()

    return args

# Generates the k-mer sets given k (size of kmer)
def generate_kmers(k):
    bases = ['A', 'C', 'T', 'G']
    kmers = [''.join(p) for p in itertools.product(bases,repeat=k)]
    return(kmers)

# Populates a pandas dataframe with the kmer and its frequency in genomes parsered 
def kmer_genome_population(kmers, L):
    kmer_count = pd.DataFrame([kmer,L[i].seq.count(kmer)] for kmer in kmers for i in range(len(L))) # if len(L[i].seq) > 6*10**6
    kmer_count.columns=['kmer', 'kmer_value']
    return(kmer_count)

def main():
    args = get_arguments()  

    # Importing the genomes as fastas and stored in list record by record (i.e by the ">" symbol)
    file = args.input
    sequences = SeqIO.parse(file, 'fasta')
    L = []
    for record in sequences:
        L.append(record)

    now = datetime.datetime.now()
    print("Current time at initialization: ")
    print(str(now.strftime('%H:%M:%S')))

    kmers = generate_kmers(args.k)
    kmer_count = kmer_genome_population(kmers,L)

    now = datetime.datetime.now()
    print("Current time after execution: ")
    print(str(now.strftime('%H:%M:%S')))

    kmer_count.to_csv(args.output)

if __name__ == '__main__':
    main()

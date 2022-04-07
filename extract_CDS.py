#!/usr/bin/env python 

import argparse
import sys 
import os
import itertools as product
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description='Extract CDS from a set of genbank files given the path to the directory and the name of the CDS (preferentially if its the exact protein name) and it outputs a single file with the multiple fasta sequences', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
inputs = parser.add_mutually_exclusive_group(required=True)

inputs.add_argument('-i', '--input', type=str,
                        help='Wirte the path to genbank files directory')
parser.add_argument('-p', '--protein', type=str,
                        help='Give a protein name or at least words that may match part of its complete name (it also finds the abreviation)')


parser.add_argument('-o', '--output', default="CDS.faa", help="Output a multfasta file", action="store")
    
args = parser.parse_args()


datadir=args.input
protein=args.protein
genomes = os.listdir(datadir)
output = open(args.output, "w")

print("\n")

for molecule in genomes:
    filename = datadir + molecule
    print("Extracting all sequences matching the '%s' input name for the %s sequence" % (protein,molecule))

    for record in list(SeqIO.parse(filename, 'genbank')):
        for feat in record.features:
            if feat.type == "CDS": # source, CDS, gene, tRNA, rRNA, misc_RNA, operon, intron, exon, mRNA, regulatory
                if protein in feat.qualifiers['product'][0]: 
                    call = feat.qualifiers['product'][0]
                    start = feat.location.start.position + 1
                    end = feat.location.end.position
                    pos = [start, end] 
                    length = end - start
                    output.write(">" + molecule.split('.')[0] + " " + "%s" % (call) + "," + " Position %d-%d" % (start, end) + ", %sbp" % (length) + "\n")
                    output.write(str(feat.extract(record.seq)) + "\n")

output.close
print("\n")

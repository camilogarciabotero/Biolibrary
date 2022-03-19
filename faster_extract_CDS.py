#!/usr/bin/env python

import argparse
import sys
import gb_io
from rich import print

def get_arguments():

    parser = argparse.ArgumentParser(description='Count the frecuency of a k-mer set  given its size (k) along a genome or a set of genomes', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', type=str, help='Wirte the path to a genbank file')
    parser.add_argument('-o', '--output', default="input-cds.fasta", help="Output a multifasta file", action="store")
    parser.add_argument('-p', '--protein', type=str, help='Give a protein name or at least words that may match part of its complete name (it also finds the abreviation)')
    
    args = parser.parse_args()

    return(args)
    
def main():
    args = get_arguments()  

    input = args.input
    protein = args.protein
    output = open(args.output, "w")

    with open(output, "w") as out:
        for record in gb_io.iter(input):
            for feature in filter(lambda feat: feat.type == "CDS", record.features):
                qualifiers = feature.qualifiers.to_dict()
                if protein in qualifiers["locus_tag"][0]:
                    start = feature.location.start.position + 1
                    end = feature.location.end.position
                    pos = [start, end] 
                    out.write(">{}\n".format(qualifiers["locus_tag"][0]))
                    out.write("{}\n".format(qualifiers["translation"][0]))


if __name__ == '__main__':
    main()

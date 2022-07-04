#!/usr/bin/env python

import argparse
from Bio import SeqIO
from BCBio import GFF

def get_arguments():

    parser = argparse.ArgumentParser(description='Convert a GFF file into GBK', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--fasta', type=str, help='Path to the FASTA file')
    parser.add_argument('-g', '--gff', type=str, help='Path to a GFF file')
    parser.add_argument('-o', '--output', default="genome.gbk", help="Output GBK file", action="store")
    
    args = parser.parse_args()

    return(args)


def main():

    args = get_arguments()

    fasta_input = args.fasta
    gff_input = args.gff
    gbk_output = args.output

    with open(gbk_output, "wt") as gbk_handler:
        fasta_handler = SeqIO.to_dict(SeqIO.parse(fasta_input, "fasta"))
        for record in GFF.parse(gff_input, fasta_handler):
            record.annotations["molecule_type"] = "DNA"
            SeqIO.write(record, gbk_handler, "genbank")

if __name__ == "__main__":
    main()



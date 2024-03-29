#!/usr/bin/env python

### This script was based on the StackExchange discussion here:
### https://bioinformatics.stackexchange.com/questions/18119/updating-the-gff3-fasta-to-genebank-code 

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from BCBio import GFF

def get_arguments():

    parser = argparse.ArgumentParser(description='Convert a GFF and a FASTA into GBK', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--fasta', type=str, help='Path to the FASTA file')
    parser.add_argument('-g', '--gff', type=str, help='Path to a GFF file')
    parser.add_argument('-o', '--output', default="genome.gbk", help="Output GBK file", action="store")
    args = parser.parse_args()

    return(args)


def convert(fasta_input, gff_input, gbk_output):
    with open(gbk_output, "wt") as gbk_handler:
        fasta_handler = SeqIO.to_dict(SeqIO.parse(fasta_input, "fasta"))
        for record in GFF.parse(gff_input, fasta_handler):
            for feature in record.features:
                if feature.strand == 1:
                    feature.qualifiers.update({'translation': Seq.translate(record.seq[feature.location.start.position:feature.location.end.position], to_stop=True)})
                else:
                    feature.qualifiers.update({'translation': Seq.translate(record.seq[feature.location.start.position:feature.location.end.position].reverse_complement(), to_stop=True)})
            record.annotations["molecule_type"] = "DNA"
            SeqIO.write(record, gbk_handler, "genbank")

def main():

    args = get_arguments()

    fasta_input = args.fasta
    gff_input = args.gff
    gbk_output = args.output

    convert(fasta_input, gff_input, gbk_output)

if __name__ == "__main__":
    main()
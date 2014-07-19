#!/usr/bin/env python
'''Reads a list of genes from a tophit file and

outputs all sequences belong to the gene.
'''

import sys
from Bio import SeqIO

def parse_tophit(tophit_file):
    '''parses gene name from a tophit file.'''

    genes = set()
    fp = open(tophit_file)
    fp.readline()
    for line in fp:
        items = line.split('\t')
        if items[1] != 'NA':
            genes.add(items[0])

    return genes

def get_seq(fasta_file, tophit_file):
    '''get sequences from a FASTA file

    and outputs only sequences belong to genes in
    a gene set.
    '''

    genes = parse_tophit(tophit_file)
    for rec in SeqIO.parse(fasta_file, 'fasta'):
        geneid = rec.id.split('.')[0]
        if geneid in genes:
            SeqIO.write(rec, sys.stdout, 'fasta')


def main():
    '''Main function'''

    fasta_file = sys.argv[1]
    tophit_file = sys.argv[2]
    get_seq(fasta_file, tophit_file)


if __name__=='__main__':
    main()


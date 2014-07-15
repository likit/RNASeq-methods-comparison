#!/usr/bin/env python
'''Description'''
import sys
from Bio import SeqIO
from collections import defaultdict

def parse(tophits):
    ensbl_genes = defaultdict(list)

    for line in tophits:
        items = line.split('\t')
        if items[1] == "NA":
            continue
        bits = float(items[3])
        geneid = items[1]
        seqid = items[0]
        ensbl_genes[geneid].append((seqid, bits))

    for key, value in ensbl_genes.iteritems():
        ensbl_genes[key] = sorted(value,
                                    key=lambda x: x[1],
                                    reverse=True)
    return ensbl_genes


def get_seq(fasta_file, tophit_genes):
    '''get sequences from a FASTA file

    and outputs only sequences belong to genes in
    a gene set.
    '''
    for rec in SeqIO.parse(fasta_file, 'fasta'):
        geneid = rec.id.split('.')[0]
        if geneid in tophit_genes:
            SeqIO.write(rec, sys.stdout, 'fasta')


def main():
    '''Main function'''

    tophits = open(sys.argv[1])  # a tophit file
    fasta_file = open(sys.argv[2])
    tophits.readline()  # skip the header

    ensbl_genes = parse(tophits)
    tophit_genes = set()
    for value in ensbl_genes.itervalues():
        tophit_genes.add(value[0][0])
    get_seq(fasta_file, tophit_genes)


if __name__=='__main__':
    main()


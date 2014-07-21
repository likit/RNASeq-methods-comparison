#!/usr/bin/env python
'''Selects annotated genes that are differential-expressed.'''
import sys
from collections import defaultdict

def parse_degenes(infile):
    '''parses genes from DE genes'''

    fp = open(infile)
    fp.readline()  # skip the header
    genes = set()
    for line in fp:
        genes.add(line.split('\t')[0].replace('"', ''))

    return genes


def select(degenes, tophits):
    '''reads top hit from a file and print genes that
    are in degenes to stdout.

    '''

    fp = open(tophits)
    fp.readline()  # skip the header

    for line in fp:
        cols = line.split('\t')
        gene = cols[0]
        ensembl = cols[1]
        if ensembl == 'NA':
            continue
        if gene in degenes:
            print '%s\t%s' % (gene, ensembl)


def main():
    '''Main function'''

    defile = sys.argv[1]  # EBseq output file
    tophits = sys.argv[2]  # a file with top hits

    degenes = parse_degenes(defile)
    select(degenes, tophits)


if __name__=='__main__':
    main()


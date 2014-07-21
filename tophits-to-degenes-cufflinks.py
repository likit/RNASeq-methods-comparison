#!/usr/bin/env python
'''Selects annotated genes that are differential-expressed.'''
import sys
from collections import defaultdict


def parse_known_isoforms(infile):
    '''create a gene-isoforms database'''

    db = {}

    for line in open(infile):
        gene, isoform = line.strip().split('\t')
        db[isoform] = gene

    return db


def parse_degenes(infile):
    '''parses genes from DE genes'''

    fp = open(infile)
    fp.readline()  # skip the header
    genes = set()
    for line in fp:
        genes.add(line.split('\t')[0].replace('"', ''))

    return genes


def select(degenes, isoformdb, tophits):
    '''reads top hit from a file and print genes that
    are in degenes to stdout.

    '''

    selects = defaultdict(list)  # select degenes

    fp = open(tophits)
    fp.readline()  # skip the header

    for line in fp:
        cols = line.split('\t')
        isoform = cols[0]
        ensembl = cols[1]
        if ensembl == 'NA':
            continue
        bits = float(cols[3])
        selects[isoformdb[isoform]].append((ensembl, bits))

    for gene, hits in selects.iteritems():
        if gene in degenes:
            print '%s\t%s' % (gene, sorted(hits,
                                key=lambda x: x[1],
                                reverse=True)[0][0])


def main():
    '''Main function'''

    defile = sys.argv[1]  # EBseq output file
    tophits = sys.argv[2]  # a file with top hits
    knownfile = sys.argv[3]  # a file with known isoforms

    isoformdb = parse_known_isoforms(knownfile)
    degenes = parse_degenes(defile)
    select(degenes, isoformdb, tophits)


if __name__=='__main__':
    main()


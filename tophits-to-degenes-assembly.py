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

    selects = defaultdict(list)
    fp = open(tophits)
    fp.readline()  # skip the header

    for line in fp:
        cols = line.split('\t')
        gene, isoform = cols[0].split()
        ensembl = cols[1]
        if ensembl == 'NA':
            continue
        bits = float(cols[3])
        selects[gene].append((ensembl, bits))

    for gene, hits in selects.iteritems():
        if gene in degenes:
            print '%s\t%s' % (gene, sorted(hits,
                                key=lambda x: x[1],
                                reverse=True)[0][0])


def main():
    '''Main function'''

    defile = sys.argv[1]  # EBseq output file
    tophits = sys.argv[2]  # a file with top hits

    degenes = parse_degenes(defile)
    select(degenes, tophits)


if __name__=='__main__':
    main()


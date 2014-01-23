import sys
from Bio import SeqIO
from collections import defaultdict

'''Finds the longest isoform of a gene in a FASTA file and writes
it to stdout.

'''

def get_longest_gene(infile):
    '''Returns a dictionary of sequence record objects of

    the longest sequences of each gene.
    '''

    db = {}
    for n, rec in enumerate(SeqIO.parse(infile, 'fasta'), start=1):
        try:
            seq = db[rec.id]
        except KeyError:
            db[rec.id] = rec
        else:
            if len(seq) < len(rec):
                db[rec.id] = rec
        if n % 1000 == 0:
            print >> sys.stderr, '...', n
    return db


if __name__=='__main__':
    for infile in sys.argv[1:]:
        for rec in get_longest_gene(infile).itervalues():
            SeqIO.write(rec, sys.stdout, 'fasta')

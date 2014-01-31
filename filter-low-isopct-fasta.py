'''Reads output from RSEM isoform expression and filter out
sequences with low isoform percentage.

'''

import sys
from Bio import SeqIO

if len(sys.argv) < 3:
    print >> sys.stderr, \
        'Usage: filter-low-isopct.py cutoff fasta isoforms-results'
    sys.exit(1)

cutoff = float(sys.argv[1])
fasta = sys.argv[2]
db = set()
for infile in sys.argv[3:]:
    fp = open(infile)
    _ = fp.readline()
    for line in fp:
        items = line.split('\t')
        transid = items[0]
        isopct = float(items[-1])

        if isopct > cutoff:
            db.add(transid)

for rec in SeqIO.parse(fasta, 'fasta'):
    rec.id = rec.id.split('/')[0]
    print rec.id
    if rec.id in db:
        SeqIO.write(rec, sys.stdout, 'fasta')

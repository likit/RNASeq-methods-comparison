'''Selects sequences of DE genes from RSEM output.'''

import sys
from Bio import SeqIO

if len(sys.argv) < 3:
    print >> sys.stderr, 'Usage: %s rsem-output fasta-file' % sys.argv[0]
    sys.exit(1)

rsemout = open(sys.argv[1])
fasta_file = sys.argv[2]
degenes = set()

_ = rsemout.readline()
for line in rsemout:
    geneid = line.split('\t')[0].replace('"', '')
    degenes.add(geneid)

for rec in SeqIO.parse(fasta_file, 'fasta'):
    geneid = '_'.join(rec.id.split('_')[:2])
    if geneid in degenes:
        SeqIO.write(rec, sys.stdout, 'fasta')

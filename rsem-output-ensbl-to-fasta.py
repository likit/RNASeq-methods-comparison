'''Selects sequences of DE genes from RSEM output.

The script is tailored specifically for sequences from Ensembl.
Fasta file has to be translated proteins from Ensembl.

'''

import sys
from Bio import SeqIO

if len(sys.argv) < 3:
    print >> sys.stderr, 'Usage: %s rsem-output fasta-file known-isoforms' \
            % sys.argv[0]
    sys.exit(1)

rsemout = open(sys.argv[1])
fasta_file = sys.argv[2]
known_isoforms = sys.argv[3]
degenes = set()

_ = rsemout.readline()
for line in rsemout:
    geneid = line.split('\t')[0].replace('"', '')
    degenes.add(geneid)

knowns = {}

for line in open(known_isoforms):
    geneid, transid = line.strip().split('\t')
    knowns[transid] = geneid

for rec in SeqIO.parse(fasta_file, 'fasta'):
    transid = rec.id
    geneid = knowns[transid]
    rec.id = geneid
    rec.description = transid
    if geneid in degenes:
        SeqIO.write(rec, sys.stdout, 'fasta')

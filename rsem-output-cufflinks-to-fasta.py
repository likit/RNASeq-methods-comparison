'''Selects sequences of DE genes from RSEM output.'''

import sys
from Bio import SeqIO

if len(sys.argv) < 3:
    print >> sys.stderr, 'Usage: %s rsem-output fasta-file gtf_file' % sys.argv[0]
    sys.exit(1)

rsemout = open(sys.argv[1])
fasta_file = sys.argv[2]
gtf_file = sys.argv[3]
degenes = set()  # differential-expressed genes
transdb = {}  # transcripts database

_ = rsemout.readline()
for line in rsemout:
    geneid = line.split('\t')[0].replace('"', '')
    degenes.add(geneid)

for line in open(gtf_file):
    attrs = line.split('\t')[-1].split()
    geneid = attrs[1].strip(';').replace('"', '')
    transid = attrs[3].strip(';').replace('"','')
    transdb[transid] = geneid

for rec in SeqIO.parse(fasta_file, 'fasta'):
    transid = rec.id
    geneid = transdb[transid]
    if geneid in degenes:
        rec.id = '%s.%s' % (geneid, transid)
        SeqIO.write(rec, sys.stdout, 'fasta')

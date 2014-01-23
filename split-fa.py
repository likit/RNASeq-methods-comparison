import sys
from Bio import SeqIO

seqs = []
chunk = int(sys.argv[2])
prefix = sys.argv[3]
for n, seq in enumerate(SeqIO.parse(sys.argv[1], 'fasta')):
    if (n > 0) and (n % int(chunk) == 0):
        SeqIO.write(seqs, '%s_%d.fa' % (prefix, n), 'fasta')
        seqs = []
    seqs.append(seq)

if seqs:
    SeqIO.write(seqs, '%s_%d.fa' % (prefix, n), 'fasta')

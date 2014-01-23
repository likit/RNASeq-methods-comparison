'''Description'''


import sys
from Bio import SeqIO

for rec in SeqIO.parse(sys.argv[1], 'fasta'):
    print rec.id
    print rec.description.split()[3].split(':')
    break

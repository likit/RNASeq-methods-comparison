'''Output Gene IDs with a flag 1 and 0 indicating DE or not DE.'''

import sys

if len(sys.argv) < 2:
    print 'Usage: rsem-to-goseq.py <matrix_file> <degenes_file>'
    sys.exit(1)

degenes = [row.split()[0] for row in open(sys.argv[2])]

for row in open(sys.argv[1]):
    geneid = row.split()[0]
    if geneid not in degenes:
        print '%s\t0' % geneid
    else:
        print '%s\t1' % geneid

'''creates knownIsoforms.txt file for RSEM from GTF file'''

import sys

for line in sys.stdin:
    items = line.strip().split('\t')[-1].split(';')
    geneid, transid = items[0], items[1]
    geneid = geneid.lstrip().replace('"', '').split(' ')[1]
    transid = transid.lstrip().replace('"', '').split(' ')[1]
    print '%s\t%s' % (geneid, transid)

import sys

op1 = open(sys.argv[1], 'w')  # fasta output file
op2 = open(sys.argv[2], 'w')  # gene - transcript ID list output file

for line in sys.stdin:
    if line.startswith('>'):
        items = line.lstrip('>').split('_')
        geneid = '%s_%s' % (items[0], items[1])
        transid = '%s_%s' % (items[2], items[3].split('/')[0])
        print >> op1, '>%s_%s' % (geneid, transid)
        print >> op2, '%s\t%s_%s' % (geneid, geneid, transid)
    else:
        print >> op1, line.strip()

op1.close()
op2.close()

import sys

for line in sys.stdin:
    if line.startswith('>'):
        tranid = line.lstrip('>').strip()
        geneid = tranid.split('.')[0]
        print '%s\t%s' % (geneid, tranid)

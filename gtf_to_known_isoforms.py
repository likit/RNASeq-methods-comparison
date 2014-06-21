import sys
import csv

reader = csv.reader(open(sys.argv[1]), dialect='excel-tab')
for rec in reader:
    if rec[2] == 'exon':
        geneid, transid = rec[-1].split('; ')[:2]
        geneid = geneid.split()[1].replace('"','')
        transid = transid.split()[1].replace('"','')
        print '%s\t%s' % (geneid, transid)

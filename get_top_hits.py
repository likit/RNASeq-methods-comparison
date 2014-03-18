import sys
from Bio.Blast import NCBIXML

infile = sys.argv[1]
print 'geneID\tscore\tbits\tE-value'
for n, blast_rec in enumerate(NCBIXML.parse(open(infile))):
    try:
        alignment = blast_rec.alignments[0]
    except IndexError:
        print '%s\tNA\tNA\tNA\tNA' % (blast_rec.query)
    else:
        geneid = alignment.title.split()[3].split(':')[1]
        hsp = alignment.hsps[0]
        print '%s\t%s\t%.2f\t%.2f\t%e' % \
                                (blast_rec.query.split('.')[0],
                                    geneid,
                                    hsp.score,
                                    hsp.bits,
                                    hsp.expect,
                                )

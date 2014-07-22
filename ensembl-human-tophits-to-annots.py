#!/usr/bin/env python
'''Description'''
import sys
import csv


def convert(topfile):
    '''Prints out a tab-delimited output with two columns.

    The first column is chicken gene IDs and the second column is
    matched human gene IDs.
    '''

    reader = csv.DictReader(open(topfile),
            fieldnames=['seqid', 'geneid', 'score', 'bits', 'evalue'],
            dialect='excel-tab')
    reader.next() # skip the header
    for row in reader:
        try:
            seqid = row['seqid'].split()[3].lstrip('gene:')
        except IndexError:
            continue  # some transcripts are novel with no gene ID

        if row['geneid'] == 'NA':
            continue  # if no match, skip
        else:
            print '%s\t%s\t%s\t%s\t%s' % (seqid,
                                            row['geneid'],
                                            row['score'],
                                            row['bits'],
                                            row['evalue'],
                                            )

def main():
    '''Main function'''

    tophits_file = sys.argv[1]
    convert(tophits_file)


if __name__=='__main__':
    main()


#!/usr/bin/env python
'''Description'''
import sys
from Bio import SeqIO
from collections import defaultdict

def parse(tophits, transcript_db):
    ensbl_genes = defaultdict(set)
    bits_db = {}

    for line in tophits:
        items = line.split('\t')
        if items[1] == "NA":
            continue
        bits = float(items[3])
        ensbl_id = items[1]
        seqid = items[0]
        geneid = transcript_db[seqid]
        ensbl_genes[ensbl_id].add(geneid)

        if geneid not in bits_db:
            bits_db[geneid] = bits


    for key, value in ensbl_genes.iteritems():
        ensbl_genes[key] = sorted(list(value),
                key=lambda x: bits_db[x],
                reverse=True)

    # for key, value in ensbl_genes.iteritems():
    #     for gene in value:
    #         print gene, bits_db[gene],
    #     print
    #     if len(value) > 2:
    #         break

    return ensbl_genes


def get_seq(fasta_file, tophit_genes):
    '''get sequences from a FASTA file

    and outputs only sequences belong to genes in
    a gene set.
    '''
    for rec in SeqIO.parse(fasta_file, 'fasta'):
        if rec.id in tophit_genes:
            SeqIO.write(rec, sys.stdout, 'fasta')

def build_db(known_isoforms):
    '''returns a transcript db'''
    transcript_db = {}
    for line in known_isoforms:
        geneid, transid = line.strip().split('\t')
        transcript_db[transid] = geneid
    return transcript_db


def main():
    '''Main function'''

    tophits = open(sys.argv[1])  # a tophit file
    fasta_file = open(sys.argv[2])
    known_isoforms = open(sys.argv[3])

    transcript_db = build_db(known_isoforms)
    tophits.readline()  # skip the header

    ensbl_genes = parse(tophits, transcript_db)
    tophit_genes = set()
    for value in ensbl_genes.itervalues():
        tophit_genes.add(value[0])
    print list(tophit_genes)[:5]
    tophit_transcripts = set()
    for trnx, gene in transcript_db.iteritems():
        if gene in tophit_genes:
            tophit_transcripts.add(trnx)

    get_seq(fasta_file, tophit_transcripts)


if __name__=='__main__':
    main()


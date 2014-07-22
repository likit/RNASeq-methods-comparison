import sys
import csv
from collections import defaultdict


def parse_pathway(annots_file):
    annots = {}

    with open(annots_file) as afile:
        reader = csv.reader(afile, dialect='excel-tab')
        for line in reader:
            geneid, pathid = line[0], line[-1]  # gene ID and pathway ID
            if geneid not in annots:
                annots[geneid] = set()
                if pathid != 'NA':
                    annots[geneid].add(pathid)
            else:
                annots[geneid].add(pathid)
    return annots

def parse_gene_annot(annot_file):
    '''Annot file is a result from BLAST top hits.'''

    gene_annots = {}
    fp = open(annot_file)
    reader = csv.DictReader(fp,
            fieldnames=["geneID", "matchID", "score", "bits", "evalue"],
            dialect='excel-tab')
    reader.next()  # skip the header
    for row in reader:
        geneid = row['geneID'].split()[0]
        matchid = row['matchID']

        if matchid == 'NA':
            continue

        bits = float(row['bits'])

        '''Find the best matched (highest bits score).'''
        if geneid not in gene_annots:
            gene_annots[geneid] = (matchid, bits)
        else:
            if bits > gene_annots[geneid][1]:
                gene_annots[geneid] = (matchid, bits)

    return gene_annots


def main():

    gene_annots = parse_gene_annot(sys.argv[1])

    first_path_annots = parse_pathway(sys.argv[2])
    second_path_annots = parse_pathway(sys.argv[3])
    combined_paths = defaultdict(set)
    # print >> sys.stderr, len(first_path_annots), len(second_path_annots)
    # print >> sys.stderr, len([len(x) for x in
    #                             first_path_annots.itervalues() if x])

    combined_pathways = defaultdict(set)
    for gene in first_path_annots:
        combined_pathways[gene].update(first_path_annots[gene])
        try:
            matchid = gene_annots[gene][0]
        except KeyError:
            continue  # some genes have no human homologs
        else:
            combined_pathways[gene].update(second_path_annots[matchid])

    for items in combined_pathways.iteritems():
        if len(items[1]) > 0:
            for pathid in items[1]:
                print '%s\t%s' % (items[0], pathid)


if __name__=='__main__':
    if len(sys.argv) < 4:
        sys.stderr.write('Error: Not enough input files\n')
        sys.stderr.write(
            'Usage: create_kegg_annots.py <annotation>' + \
            '<pathways1> <pathways2>\n')
        raise SystemExit
    else:
        main()

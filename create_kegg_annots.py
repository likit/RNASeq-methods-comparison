import sys
import csv
from collections import defaultdict


def parse(annot_file, pathway_file):
    annots = {}
    pathways = defaultdict(set)

    with open(annot_file) as annot_file:
        reader = csv.reader(annot_file, dialect='excel-tab')
        for geneid, annot in reader:
            if annot == 'NA':
                continue
            else:
                annots[geneid] = annot

    with open(pathway_file) as pathway_file:
        reader = csv.reader(pathway_file, dialect='excel-tab')
        for line in reader:
            annot, pathid = line[0], line[-1]
            if pathid == 'NA':
                continue
            else:
                if annot not in pathways:  # ignore redundant annotations
                    pathways[annot].add(pathid)
    return annots, pathways


def main():
    '''Annots are annotations of genes such as Ensembl.
    Genes are gene models such as those from assembly or Cufflinks.
    
    '''

    first_annots, first_pathways = parse(sys.argv[1], sys.argv[2])
    second_annots, second_pathways = parse(sys.argv[3], sys.argv[4])
    combined_paths = defaultdict(set)

    for gene in first_annots.keys():
        first_annot = first_annots[gene]
        if (first_annot in first_pathways and
                first_annot not in combined_paths):
            combined_paths[first_annot].update(first_pathways[first_annot])
        try:
            second_annot = second_annots[gene]
        except KeyError:
            pass
        else:
            combined_paths[first_annot].update(second_pathways[second_annot])

    with open(sys.argv[5]) as fp:
        reader = csv.reader(fp, dialect='excel-tab')
        for item in reader:
            if item[-1] != 'NA':
                combined_paths[item[0]].add(item[-1])

    for item in combined_paths.iteritems():
        if len(item[1]) == 0:
            continue
        else:
            for path in item[1]:
                print '%s\t%s' % (item[0], path)


if __name__=='__main__':
    if len(sys.argv) < 5:
        sys.stderr.write('Error: Not enough input files\n')
        sys.stderr.write(
            'Usage: create_kegg_annots.py <annotation> <pathways> ...\n')
        raise SystemExit
    else:
        main()

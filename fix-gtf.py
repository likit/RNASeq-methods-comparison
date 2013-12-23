'''Replaces "." strand in GTF file to "+".
Input is read from stdin.

'''

import sys

for line in sys.stdin:
    features = line.rstrip('\n').split('\t')
    if features[6] == ".":
        features[6] = "+"
    print "\t".join(features)

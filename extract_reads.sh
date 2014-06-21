#!/bin/sh

function usage {
    printf "Usage: %s <BAM file> <chromosomes.txt>\n" $0
    printf "A BAM file is required.\n"
} >&2

if [ $# -lt 2 ]
then
    usage
    exit 1
fi

if [ ! -e "$1".bai ]
then
    printf "BAM index file is not found, creating it...\n" >&2
    samtools index "$1"
fi

outdir=$(dirname "$1")

for chr in $(cat "$2")
do
    printf "extracting reads for %s..." "$chr" >&2
    samtools view -b -o "$outdir"/"$chr".bam "$1" "$chr"
    printf "done.\n" >&2
done

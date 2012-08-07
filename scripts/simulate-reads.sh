#!/bin/bash

ART_BIN=""
FASTA_IN=""
OUT_PREFIX="fusim-reads"
HELP=0

while getopts "hn:a:f:o:" o ; do
    case $o in
        a) ART_BIN=$OPTARG;;
        f) FASTA_IN=$OPTARG;;
        o) OUT_PREFIX=$OPTARG;;
        h) HELP=1;;
    esac
done

if [[ $HELP -eq 1 || ! -e $FASTA_IN ]]; then
    echo "usage: simulate-reads.sh "
    echo "  -a  path to ART binary"
    echo "  -f  fasta file of simulated fusions"
    echo "  -o  output prefix"
    echo "  -h  help"
    exit 1
fi

if [[ -z "$ART_BIN" ]]; then
    ART_BIN=`which art_illumina 2> /dev/null`
    if [[ -z "$ART_BIN" ]]; then
        echo "Can't find an ART executable. Please provide a valid path to ART."
        exit 1
    fi
fi

if [[ ! -x $ART_BIN ]]; then
    echo "$ART_BIN is not executable. Please provide a valid path to ART."
    exit 1
fi


echo "Simulating reads with ART.."
$ART_BIN -i $FASTA_IN -o $OUT_PREFIX \
         -l 75 -f 10 -p -m 400 -s 10

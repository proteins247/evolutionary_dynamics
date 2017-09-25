#!/usr/bin/env python
"""
generate_translation_params.py

Takes translation time data from Ciryam et al. (2013) PNAS, multiplies
the times by a scaling factor, and outputs it as (codon, time) pairs
for use as a parameter file for folding_evolution.

The codons are represented as hexadecimal numbers.

@author - Victor Zhao

"""

import sys
import datetime
import argparse
import numpy as np

translation_times = {           # in ms
    "UUU": 99.0  ,
    "UUC": 156.4 ,
    "UUG": 42.1  ,
    "UUA": 168.0 ,
    "UCU": 34.3  ,
    "UCC": 184.5 ,
    "UCG": 81.5  ,
    "UCA": 65.8  ,
    "UGU": 59.5  ,
    "UGC": 97.0  ,
    "UGG": 117.3 ,
    "UGA": 39.9  ,
    "UAU": 38.1  ,
    "UAC": 63.4  ,
    "UAG": 184.9 ,
    "UAA": 35.0  ,
    "CUU": 178.5 ,
    "CUC": 136.1 ,
    "CUG": 28.9  ,
    "CUA": 254.5 ,
    "CCU": 116.3 ,
    "CCC": 61.9  ,
    "CCG": 156.0 ,
    "CCA": 194.2 ,
    "CGU": 20.4  ,
    "CGC": 28.3  ,
    "CGG": 483.5 ,
    "CGA": 26.1  ,
    "CAU": 190.3 ,
    "CAC": 138.0 ,
    "CAG": 134.8 ,
    "CAA": 125.0 ,
    "GUU": 18.3  ,
    "GUC": 175.2 ,
    "GUG": 30.8  ,
    "GUA": 66.4  ,
    "GCU": 23.9  ,
    "GCC": 362.3 ,
    "GCG": 27.7  ,
    "GCA": 61.9  ,
    "GGU": 25.7  ,
    "GGC": 40.9  ,
    "GGG": 65.6  ,
    "GGA": 321.1 ,
    "GAU": 53.5  ,
    "GAC": 86.8  ,
    "GAG": 28.0  ,
    "GAA": 48.4  ,
    "AUU": 66.4  ,
    "AUC": 88.8  ,
    "AUG": 213.8 ,
    "AUA": 90.7  ,
    "ACU": 32.7  ,
    "ACC": 134.6 ,
    "ACG": 92.2  ,
    "ACA": 96.9  ,
    "AGU": 71.3  ,
    "AGC": 137.0 ,
    "AGG": 416.0 ,
    "AGA": 189.7 ,
    "AAU": 66.6  ,
    "AAC": 114.6 ,
    "AAG": 94.7  ,
    "AAA": 59.1
    }

letter_to_num = {
    'U' : 0,
    'C' : 1,
    'A' : 2,
    'G' : 3
    }


def codon_triplet_to_int(triplet):
    n1 = letter_to_num[triplet[0]]
    n2 = letter_to_num[triplet[1]]
    n3 = letter_to_num[triplet[2]]
    return (n1 << 8 ) | (n2 << 4) | n3;

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scaling', type=float, default=1.)
    return parser.parse_args()

def main(args):
    print "# Codon translation times"
    print "# Scaling factor:", args.scaling
    print "# Generated on:", datetime.date.today()
    sys.stdout.write('# ')
    print letter_to_num
    times = []
    for i, (codon, time) in enumerate(sorted(
            translation_times.items(), key=lambda x: codon_triplet_to_int(x[0]))):
        if (i % 4) == 0:
            print
        codon_as_int = codon_triplet_to_int(codon)
        scaled_time = int(time) * args.scaling
        times.append(scaled_time)
        print "0x{:04x} {:.0f}".format(codon_as_int, scaled_time)
    print
    print "# mean time:", np.mean(times)
    print "# median time:", np.median(times)
    return

if __name__ == "__main__":
    main(parse_args())

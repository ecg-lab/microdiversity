#!/usr/bin/env python
"""Get nt frequencies for sequence files
"""

import sys
import imp
try:
    imp.find_module('Bio')
except ImportError:
    print "No Bio Module. Wrong Python verison?"
    sys.exit(1)

import numpy

from Bio import SeqIO
from optparse import OptionParser

def parse_fasta(file, format='fasta'):
    """Get count of nt's for a sequence file
    """
    if not file:
        return
    freq = { 'A': 0, 'C': 0, 'G': 0, 'T': 0 }
    with open(file, 'rU') as f:
        seqs = SeqIO.parse(f, format)
        if not seqs:
            return
        for seq in seqs:
            for nt in freq.keys():
                freq[nt] += seq.seq.upper().count(nt)
    
    return freq

def parse_files(files=[], format='fasta'):
    """Get nt frequencies for several files
    """
    if not files:
        return
    
    freq = { 'A': 0, 'C': 0, 'G': 0, 'T': 0 }
    for file in files:
        temp = parse_fasta(file, format)
        for nt in temp.keys():
            freq[nt] += temp[nt]
    
    total = float(sum(freq.values()))
    if total == 0:
        print "Error getting sum"
        return
    for nt in freq.keys():
        freq[nt] /= total
    
    return freq

def main():
    parser = OptionParser()
    parser.add_option('-f', '--format', dest='format', default='fasta',
                      help='Format of sequence files. Default fasta')
    (options, args) = parser.parse_args()
    
    freq = parse_files(args, options.format)
    A = "{:.4f}".format(freq['A'])
    C = "{:.4f}".format(freq['C'])
    G = "{:.4f}".format(freq['G'])
    T = "{:.4f}".format(freq['T'])
    if float(A) + float(C) + float(G) + float(T) == 1:
        print A, C, G, T
    else:
        print >>sys.stderr, "WARNING! Frequencies when rounded do not add up to 1."
        print freq['A'], freq['C'], freq['G'], freq['T']
    
    
if __name__ == "__main__":
    main()

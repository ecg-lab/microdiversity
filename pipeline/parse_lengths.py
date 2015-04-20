#!/usr/bin/env python

import sys
import imp
try:
    imp.find_module('Bio')
except ImportError:
    print "No Bio Module. Wrong Python verison?"
    sys.exit(1)

import re
import cPickle
import gzip
import numpy

from Bio import SeqIO
from optparse import OptionParser

def parse_fasta(file, ungap=False):
    with open(file, 'rU') as f:
        seqs = {}
        temp = SeqIO.parse(f, 'fasta')
        for record in temp:
            if record.id not in seqs:
                seqs[record.id] = record
            else:
                print "Found duplicate gene id: %s" % record.id
    lengths = {}
    if seqs:
        GI_MATCH = re.compile(r'[fg]id?\|(\d+)')
        for seq in seqs:
            match = GI_MATCH.match(seq)
            if match:
                key = match.group(1)
            else:
                key = seq
            if ungap:
                length = len(seqs[seq].seq.ungap('-'))
            else:
                length = len(seqs[seq].seq)
            if key not in lengths:
                lengths[key] = length
            elif key in lengths and lengths[key] != length:
                print "Warnings: Same GI, but different length"
    return lengths

def parse_families(files, ungap=False):
    lengths = {}
    for file in files:
        lengths[file] = numpy.mean(parse_fasta(file, ungap).values())
    return lengths
    
def parse_genomes(files, ungap=False):
    lengths = {}
    for file in files:
        lengths.update(parse_fasta(file, ungap))
    return lengths

def save_lengths(lengths, file, pickle=True):
    if pickle:
        with gzip.open(file, 'wb') as w:
            cPickle.dump(lengths, w, -1)
    else:
        with open(file, 'wb') as w:
            w.write('\n'.join(['%s\t%i' % (seq, lengths[seq]) for seq in sorted(lengths)]))
            
def load_lengths(file, pickle=True):
    if pickle:
        with gzip.open(file, 'rb') as f:
            lengths = cPickle.load(f)
    else:
        with open(file, 'rb') as f:
            lengths = {}
            for line in f:
                temp = line.split('\t')
                lengths[temp[0]] = int(temp[1])
    return lengths

def filter_mcl_file(file, output, lengths={}):
    with open(file, 'rb') as r:
        with open(output, 'wb') as w:
            for line in r:
                genes = line.split()
                cluster_lengths = [lengths[gene] for gene in genes if gene in lengths]
                mean = numpy.mean(cluster_lengths)
                std = numpy.std(cluster_lengths)
                for gene in genes:
                    if gene in lengths:
                        if abs(lengths[gene] - mean) <= 2 * std:
                            continue
                    genes.remove(gene)
                w.write('\t'.join(genes)+'\n')

def main():
    parser = OptionParser(usage="\n\t%prog [options] -b -l LIBRARY file[...]\n\n\t%prog [options] -l LIBRARY -o OUT file[...]")
    parser.add_option('-o', '--out', dest='out', default='filter',
                      help='Specify extension to append to output file')
    parser.add_option('-l', '--lengths', dest='lengths',
                      help='Load or save a built dictionary of lengths')
    parser.add_option('-t', '--text', dest='text', action='store_true', default=False,
                      help='Save/load built dictionary to/from text file')
    parser.add_option('-b', '--build', dest='build', action='store_true', default=False,
                      help='Build a dictionary of sequence lengths (specify genomes as arguments)')
    parser.add_option('-f', '--family', dest='family', action='store_true', default=False,
                      help='Get average lengths for a set of gene families (will save results as text)')
    parser.add_option('-c', '--convert', dest='convert', action='store_true', default=False,
                      help='Convert saved length files between text and pickled gzipped file')
    parser.add_option('-u', '--ungap', dest='ungap', action='store_true', default=False,
                      help='Ungap alignment before calculating length')
    (options, args) = parser.parse_args()
    
    if not options.lengths:
        print "No lengths dictionary file specified"
        return
    
    if options.build:
        save_lengths(parse_genomes(args, options.ungap), options.lengths, not options.text)
        print "Saved dictionary to %s" % options.lengths
    elif options.family:
        save_lengths(parse_families(args, options.ungap), options.lengths, pickle=False)
        print "Saved dictionary to %s" % options.lengths
    elif options.convert:
        lengths = load_lengths(options.lengths, not options.text)
        if options.text:
            new = options.lengths + '.pkl.gz'
        else:
            new = options.lengths + '.txt'
        save_lengths(lengths, new, options.text)
    else:
        lengths = load_lengths(options.lengths, not options.text)
        for file in args:
            filter_mcl_file(file, '%s.%s'%(file, options.out), lengths)
    
if __name__ == "__main__":
    main()

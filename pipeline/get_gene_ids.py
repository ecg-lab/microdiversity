#!/usr/bin/env python

import sys
import imp
try:
    imp.find_module('Bio')
except ImportError:
    print "No Bio Module. Wrong Python verison?"
    sys.exit(1)

import re

from Bio import SeqIO
from optparse import OptionParser

GI_MATCH = re.compile(r'[fg]id?\|(\d+)')
LOCUS_MATCH = re.compile(r'locus\|(.*?)\|')

def get_gene_ids(filename="", format="fasta", locus=False):
    assert filename
    assert format
    genes = SeqIO.parse(filename, format)
    
    ids = set()
    for gene in genes:
        if locus:
            match = LOCUS_MATCH.search(gene.id)
        else:
            match = GI_MATCH.match(gene.id)
        if match:
            temp = match.group(1)
            if temp:
                ids.add(temp)
    
    

    return list(ids)
    
def parse_files(files=[], format="fasta", locus=False):
    assert files
    assert format
    
    ids = {}
    for file in files:
        ids[file] = get_gene_ids(file, format, locus)
    
    return ids

def write_id_file(filename="", ids={}, sep="\t"):
    assert filename
    assert ids
    assert sep
    
    with open(filename, 'wb') as w:
        for name in ids.keys():
            w.write('%s%s'%(name, sep)+sep.join(ids[name])+'\n')

def main():
    parser = OptionParser(usage="Usage: %prog [options] -o OUTPUT file [...]")
    parser.add_option('-f', '--format', dest='format', default='fasta',
                      help="Format of the gene files (default is fasta)")
    parser.add_option('-l', '--locus', dest='locus', action="store_true", default=False,
                      help="Output locus IDs instead of gene IDs")
    parser.add_option('-o', '--out', dest='out',
                      help="Output filename")
    parser.add_option('-s', '--sep', dest='sep', default='\t',
                      help='Separator in output file (default tab)')
    (options, args) = parser.parse_args()
    
    if not args:
        print "No files specified"
        return
        
    if not options.out:
        print "No output file specified"
        return
    
    if not options.sep:
        options.sep = '\t'
    
    print "Getting IDs for %i files" % len(args)
    ids = parse_files(args, options.format, options.locus)
    print "Found %i gene IDs" % (sum([len(x) for x in ids.values()]))
    write_id_file(options.out, ids, options.sep)
    print "Wrote gene IDs to %s" % options.out
    
if __name__ == "__main__":
    main()

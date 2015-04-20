#!/usr/bin/env python

import re
import os
import sys
import numpy

from optparse import OptionParser

def cluster_stats(file):
    assert file
    
    with open(file) as f:
        lengths = [len(line.split()) for line in f]
    
    print "There were a total of %i clusters" % len(lengths)
    print "The mean cluster size was %g" % numpy.mean(lengths)
    print "The SD of cluster size was %g" % numpy.std(lengths)
    print "The median cluster size was %g" % numpy.median(lengths)
    
    print "The max cluster size was %g" % max(lengths)
    print "The min cluster size was %g" % min(lengths)
    print "The number of clusters of <= 2 size were %g" % len([x for x in lengths if x <= 2])
    
    split_size = max(lengths) / 10
    for i in range(1, 10):
        size = split_size*i
        print "# of clusters >= %i: %i" % (size, len([x for x in lengths if x >= size]))
    
def filter_clusters(file, minimum=0, maximum=-1, n=-1):
    assert file
    assert type(minimum) == int
    assert type(maximum) == int
    
    if maximum > 0:
        assert maximum >= minimum
    
    
    if n > 0:
        with open(file) as f:
            lengths = {}
            for i, line in enumerate(f):
                lengths[i] = len(line.split())
        
        if n >= len(lengths):
            print "Note: n was set higher than the number of clusters in the file"
            with open(file) as f:
                return f.readlines()
        
        keep = sorted(lengths.keys(), key=lengths.get, reverse=True)[0:n]
        
        with open(file) as f:
            filtered = []
            for i, line in enumerate(f):
                if i in keep:
                    filtered.append(line)
            
            return filtered
    
    if maximum <= 0:
        with open(file) as f:
            return [line for line in f if len(line.split()) >= minimum]
    else:
        with open(file) as f:
            filtered = []
            for line in f:
                length = len(line.split())
                if length >= minimum and length <= maximum:
                    filtered.append(line)
            
            return filtered

def filter_out_genes(lines, gene_file):
    assert lines
    assert gene_file
    
    with open(gene_file, 'rb') as f:
        genes = []
        for line in f:
            temp = [x.strip('>;,\'"') for x in line.split()]
            genes.extend(temp)
        
        genes = set(genes)
        
    for line in lines:
        temp = line.split()
        for gene in temp:
            if gene in genes:
                lines.remove(line)
                break
                
    return lines

def keep_genes(lines, gene_file):
    assert lines
    assert gene_file
    
    with open(gene_file, 'rb') as f:
        genes = []
        for line in f:
            temp = [x.strip('>;,\'"') for x in line.split()]
            genes.extend(temp)
        
        genes = set(genes)
    
    keep = []
    for line in lines:
        temp = genes.intersection(line.split())
        if temp: keep.append('\t'.join(temp)+'\n')
    
    return keep

def write_file(file, lines):
    assert file
    assert lines
    
    with open(file, 'w') as w:
        for line in lines:
            w.write(line)
        return True

def main():
    parser = OptionParser(usage="\n\t%prog [options] file [...]\n\n\t%prog [options] -o OUT file\n\n\t%prog -s file [...]")
    parser.add_option("-o", "--out", dest="out", default=None,
                      help="Specify output filename. Only works if input is one file")
    parser.add_option("-m", "--min", dest="min", type="int", default=0,
                      help="Specify minimum cluster size to keep")
    parser.add_option("-x", "--max", dest="max", type="int", default=-1,
                      help="Specify maximum cluster size to keep")
    parser.add_option("-n", "--num", dest="num", type="int", default=-1,
                      help="Keep the largest n clusters (incompatible with min and max)")
    parser.add_option("-f", "--filter", dest="filter",
                      help="Specify file that specifies genes to filter out")
    parser.add_option("-k", "--keep", dest="keep",
                      help="Specify file of gene IDs to keep in (opposite of filter)")
    parser.add_option("-s", "--stats", dest="stats", action="store_true", default=False,
                      help="Output only statistics about MCL file. Do not filter.")
    (options, args) = parser.parse_args()
    
    for file in args:
        if options.stats:
            cluster_stats(file)
            return
        
        try:
            (name, ext) = re.match(r"(.*)\.(.*)", file).groups()
        except:
            name = file
            ext = 'mcl'
        if len(args) == 1 and options.out:
            output = options.out
        elif options.num > 0:
            output = '%s.top%i.%s'%(name, options.num, ext)
        elif options.max > 0 and options.min > 0:
            output = '%s.min%i.max%i.%s'%(name, options.min, options.max, ext)
        elif options.max <= 0 and options.min > 0:
            output = '%s.min%i.%s'%(name, options.min, ext)
        elif options.max > 0 and options.min <= 0:
            output = '%s.max%i.%s'%(name, options.max, ext)
        else:
            output = '%s.out.%s' % (name, ext)
        lines = filter_clusters(file, options.min, options.max, options.num)
        
        if options.filter:
            print "Prior to filtereing, there were %i clusters" % len(lines)
            lines = filter_out_genes(lines, options.filter)
            print "After filtering out genes, there are now %i clusters" % len(lines)
            output = output.replace(ext, "%s.%s" % (options.filter, ext))
        
        if options.keep:
            print "There were originally %i clusters" % len(lines)
            lines = keep_genes(lines, options.keep)
            print "Keeping only selected genes, there are now %i clusters" % len(lines)
            output = output.replace(ext, "%s.%s" % (options.keep, ext))
        
        write_file(output, lines)
        print "Output filtered clusters to %s"%output
    
if __name__ == "__main__":
    main()

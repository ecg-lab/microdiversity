#!/usr/bin/env python
"""Perform some filtering and modifying on MCL clusters
"""

import numpy

import get_gene_ids
from optparse import OptionParser

def read_gene_ids(file='', sep='\t'):
    """Read gene IDs from a text file
    """
    ids = {}
    with open(file, 'rb') as f:
        for line in f:
            temp = line.split(sep)
            ids[temp[0]] = temp[1:]
    
    return ids

def read_mcl_file(file=''):
    """Read clusters from the output of MCL
    """
    with open(file, 'rb') as f:
        clusters = [line.split() for line in f]
    
    return clusters

def write_mcl_file(file='', clusters=[]):
    """Write out clusters to a text file similar to MCL output
    """
    with open(file, 'wb') as w:
        for cluster in clusters:
            w.write('\t'.join([str(x) for x in cluster])+'\n')
    
def check_clusters(clusters = [], genomes = {}, remove=False, exclude_paralogs=-1, threshold=0):
    """Iterate through clusters, optionally removing paralogs or small clusters
    """
    results = {}
    
    for i, cluster in enumerate(clusters):
        results[i] = {}
        temp_remove = []
        tempset = set(cluster)
        for genome in genomes:
            temp = tempset.intersection(genomes[genome])
            if temp:
                results[i][genome] = len(temp)
                if remove and len(temp) > 1:
                    temp_remove.extend(temp)
        if remove:
            clusters[i] = list(tempset.difference(temp_remove))
    
    if remove:
        clusters = filter(None, clusters)
    
    if threshold:
        remove_me = []
        for i in results.keys():
            if len(results[i].keys()) < threshold:
                remove_me.append(clusters[i])
                #del results[i]
                
        for cluster in remove_me:
            clusters.remove(cluster)
    
    if exclude_paralogs > 0:
        remove_me = []
        for i in results.keys():
            if numpy.mean(results[i].values()) > exclude_paralogs:
                remove_me.append(clusters[i])
        
        for cluster in remove_me:
            clusters.remove(cluster)
                
    return results, clusters
    
def genome_count(results = {}, file='', sep='\t'):
    """Determine the number of genomes represented in the clusters
    """
    info = {}
    for cluster in results:
        info[cluster] = {}
        info[cluster]['genome'] = len(results[cluster].keys())
        info[cluster]['mean'] = numpy.mean(results[cluster].values())

    if file:
        with open(file, 'wb') as w:
            w.write(sep.join(['Cluster', 'Genome', 'Mean'])+'\n')
            for cluster in sorted(info.keys()):
                w.write('%i%s'%(cluster, sep))
                w.write('%i%s'%(info[cluster]['genome'], sep))
                w.write('%f\n'%(info[cluster]['mean']))
    else:
        print sep.join(['Cluster','Genome','Mean'])
        for cluster in sorted(info.keys()):
            print '%i\t%f\t%f'%(cluster, info[cluster]['genome'], info[cluster]['mean'])

def filter_cluster(filter, cluster, out):
    """Filter clusters based on a text file
    
    I don't know what this function is for anymore.
    """
    with open(filter, 'rb') as f:
        filter = set([int(x) for x in f])
    
    with open(cluster, 'rb') as f:
        clusters = []
        for i, line in enumerate(f):
            if i in filter:
                clusters.append(line)
    
    if clusters:
        with open(out, 'wb') as w:
            w.write(''.join(clusters))
            
def main():
    parser = OptionParser(usage="\n\t%prog [options] -m MCL -o OUTPUT file [...]\n\n\t%prog [options] -m MCL -i INPUT")
    parser.add_option('-m', '--mcl', dest='mcl',
                      help='Specify mcl cluster file')
    parser.add_option('-i', '--in', dest='input',
                      help='Specify file to import gene IDs from')
    parser.add_option('-f', '--format', dest='format', default='fasta',
                      help="Format of the gene files (default is fasta)")
    parser.add_option('-l', '--locus', dest='locus', action="store_true", default=False,
                      help="Output locus IDs instead of gene IDs")
    parser.add_option('-o', '--out', dest='out',
                      help="Output filename for gene IDs")
    parser.add_option('-t', '--threshold', dest='threshold', type='int', default=None,
                      help='Specify threshold for minimum genome count')
    parser.add_option('-c', '--clusters', dest='clusters', default='',
                      help='Specify output filename for cluster genome count')
    parser.add_option('-e', '--exclude', dest='exclude', type='float', default=-1,
                      help='Exclude gene families with average gene count above this')
    parser.add_option('-r', '--remove', dest='remove', action='store_true', default=False,
                      help='Remove paralogs from clusters, then save clusters (to new file)')
    parser.add_option('-s', '--sep', dest='sep', default='\t',
                      help='Separator in input/output file (default tab)')
    (options, args) = parser.parse_args()
    
    if not options.mcl:
        print "No MCL file specified"
        return
    
    if not options.sep:
        options.sep = '\t'
    
    if options.input:
        ids = read_gene_ids(options.input, options.sep)
        print "Found %i genomes in file" % len(ids.keys())
    elif not args:
        print "No input file or fasta files specified"
        return
    else:
        ids = get_gene_ids.parse_files(args, options.format, options.locus)
        if ids and options.out:
            get_gene_ids.write_id_file(options.out, ids, options.sep)
        
    clusters = read_mcl_file(options.mcl)
    print "Found %i clusters in file" % len(clusters)
    results, clusters = check_clusters(clusters, ids, options.remove, options.exclude, options.threshold)
    
    try:
        genome_count(results, options.clusters)
    except Exception as e:
        print e
    
    if options.remove or options.threshold or options.exclude > 0:
        try:
            write_mcl_file('%s.filtered'%options.mcl, clusters)
        except Exception as e:
            print e
    
if __name__ == "__main__":
    main()
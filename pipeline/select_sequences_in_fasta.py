#!/usr/bin/env python

import re
import sys

from optparse import OptionParser

import imp
try:
    imp.find_module('Bio')
except ImportError:
    print "No Bio Module. Wrong Python verison?"
    sys.exit(1)

from Bio import SeqIO

GI_MATCH = re.compile(r'[fg]id?\|(\d+)')

def get_gi(string=""):
    """Tries to find GI in string. Returns GI if found. Otherwise returns string
    """
    if not(string and type(string) == str):
        print "String was blank"
        return ""
    search = GI_MATCH.search(string)
    if search and search.group(1):
        return search.group(1)
    else:
        return string

def get_selection(filename):
    selection = []
    with open(filename) as f:
        for line in f:
            selection.append(get_gi(line.strip()))
    return selection

def get_mcl_clusters(filename):
    assert filename
    clusters = []
    with open(filename) as f:
        for line in f:
            clusters.append(line.split())
    
    return clusters

def get_keys(database):
    return dict((get_gi(record), record) for record in database.keys())
    
def load_database(file):
    assert file
    with open(file, 'rU') as raw:
        #return dict((get_gi(str(record.id)), record) for record in SeqIO.parse(raw, "fasta"))
        database = {}        
        for record in SeqIO.parse(raw, 'fasta'):
            temp = get_gi(str(record.id))
            if get_gi(str(record.id)) not in database:
                database[temp] = record
        return database, get_keys(database)
    #database = SeqIO.index(file, 'fasta')
    #return database, get_keys(database)

def select_genes(database, keys, output, selection=[]):
    assert database
    assert keys
    assert output
    assert selection
    
    sequences = [database[keys[key]] for key in selection if key in keys and keys[key] in database]
    
    with open(output, 'w') as w:
        SeqIO.write(sequences, w, "fasta")
        
    print "Wrote %i sequences to %s" % (len(sequences), output)
    return True


def main():
    parser = OptionParser(usage='\n\t%prog [options] -i INPUT -o OUTPUT\n\n\t%prog [options] -n -i INPUT')
    parser.add_option("-i", "--input", dest="input",
                      help="Input fasta file")
    parser.add_option("-o", "--output", dest="output", 
                      help="Output fasta file")
    parser.add_option("-s", "--selection", dest="selection",
                      help="Selection file")
    parser.add_option("-m", "--mcl", dest="mcl",
                      help="Use MCL output to select genes")
    parser.add_option('-n', '--nogs', dest='nogs', action='store_true', default=False,
                      help='These are NOGs. Name files based on NOG of cluster')
    (options, args) = parser.parse_args()
    
    if not options.input:
        print "No input"
        return
    if not options.nogs:
        if not options.output:
            print "No output"
            return
    
    print "Loading %s into memory" % options.input
    database, keys = load_database(options.input)
    
    assert database
    
    if options.selection:
        selection = get_selection(options.selection)
        if len(selection) == 0:
            print "No genes found in selection file."
            return
        else:
            print "Selection containts %i genes" % len(selection)
    elif options.mcl:
        clusters = get_mcl_clusters(options.mcl)
        print "Found %i clusters in MCL file" % len(clusters)
        out = 0
        if not options.nogs:
            match = re.match(r'(.*)\.(.*)', options.output)
            if match:
                (name, ext) = match.groups()
            else:
                name = options.output
                ext = 'fa'
        for i, cluster in enumerate(clusters):
            try:
                if options.nogs:
                    output = cluster.pop(0)+'.fa'
                else:
                    output = "%s_%i.%s"%(name, i, ext)
                #print "Output sequences to %s" % output
                select_genes(database, keys, output, cluster)
            except Exception as e:
                print e
                continue
            out += 1
        print "Wrote %i fasta files" % out
        
        return
    else:
        if args and len(args) > 0:
            selection = args
        else:
            print "No genes selected. Program will exit."
            sys.exit(1)
            
    assert options.output
    select_genes(database, keys, options.output, selection)
    
if __name__ == '__main__':
    main()

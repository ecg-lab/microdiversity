#!/usr/bin/env python

import sys
import imp
try:
    imp.find_module('Bio')
except ImportError:
    print "No Bio Module. Wrong Python verison?"
    sys.exit(1)

import numpy
from Bio import SeqIO
from itertools import combinations
from multiprocessing import Pool, TimeoutError
from optparse import OptionParser

def read_file(file, format='fasta'):
    with open(file, 'rU') as f:
        seqs = SeqIO.to_dict(SeqIO.parse(f, format))
        return seqs

def get_nt_pos(seqs={}):
    pos = {}
    for seq in seqs:
        pos[seq] = set([i for (i,x) in enumerate(list(str(seqs[seq].seq))) if x != '-'])
    return pos

def get_shared_pos(seq1, seq2):
    return len(seq1 & seq2)

def compare_seqs(file, format, out, threshold=50):
    seqs = read_file(file, format)
    print "Found %i sequences in %s" % (len(seqs.keys()), file)
    pos = get_nt_pos(seqs)
    coverage = {}
    for (seq1, seq2) in combinations(seqs.keys(), r=2):
        shared = get_shared_pos(pos[seq1], pos[seq2])
        if seq1 not in coverage:
            coverage[seq1] = {}
        if seq2 not in coverage:
            coverage[seq2] = {}
        coverage[seq1][seq2] = shared
        coverage[seq2][seq1] = shared

    for seq in sorted(coverage.keys(), key=lambda x: len([a for a in coverage[x].values() if a >= threshold])):
        if len(coverage[seq].values()) == 0 or min(coverage[seq].values()) < threshold:
            #print seq,
            #print len([a for a in coverage[seq].values() if a >= threshold])
            del seqs[seq]
            for seq2 in coverage[seq].keys():
                del coverage[seq2][seq]
    
    if len(seqs.values()) > 0:
        with open(out, 'wb') as w:
            SeqIO.write(seqs.values(), w, 'fasta')
            print "Wrote %i sequences to %s" % (len(seqs.keys()), out)
    else:
        print "Removed all sequences"
    
# def filter_seq(file, format, out, threshold=50):
    # seqs = read_file(file, format)
    # print "Found %i sequences in %s" % (len(seqs.keys()), file)
    # pos = get_nt_pos(seqs)
    # for (seq1, seq2) in combinations(seqs.keys(), r=2):
        # if get_shared_pos(pos[seq1], pos[seq2]) < threshold:
            # if seq1 in seqs: del seqs[seq1]
            # if seq2 in seqs: del seqs[seq2]
    # if len(seqs.values()) > 0:
        # with open(out, 'wb') as w:
            # SeqIO.write(seqs.values(), w, 'fasta')
            # print "Wrote %i sequences to %s" % (len(seqs.keys()), out)
    # else:
        # print "Removed all sequences"

def multi_wrapper(arguments):
    #filter_seq(*arguments)
    compare_seqs(*arguments)
        
def multithread(files, format, out, threshold=50, threads=1):

    if not (threads or type(threads) == int):
        print "Invalid thread count, set to 1"
        threads = 1
    
    if threads < 1:
        print "Will use all available CPUs"
        threads = None
    
    if threads == 1:
        for file in files:
            try:
                #filter_seq(file, format, '%s.%s'%(file, out), threshold)
                compare_seqs(file, format, '%s.%s'%(file, out), threshold)
            except KeyboardInterrupt:
                break
            except Exception as e:
                print e
    else:
        try:
            p = None
            p = Pool(threads)
            args = [(file, format, '%s.%s'%(file, out), threshold) for file in files]
            workers = p.map_async(multi_wrapper, args)
            while not workers.ready():
                try:
                    workers.get(5)
                except TimeoutError:
                    pass
                except KeyboardInterrupt:
                    if p: p.terminate()
                    break
                except Exception as e:
                    print e
                    if p: p.terminate()
                    break
            if p: p.close()
        except KeyboardInterrupt:
            if p: p.terminate()
        except Exception as e:
            print e
            if p: p.terminate()
 
def main():
    parser = OptionParser(usage='Usage: %prog [options] file[...]')
    parser.add_option('-t', '--threshold', dest='threshold', type='int', default=50,
                      help='Filter out sequences with less than this many non-gap overlaps, default is 50')
    parser.add_option('-o', '--out', dest='out', default='out',
                      help='Specify extension to append to output file(s)')
    parser.add_option('-f', '--format', dest='format', default='fasta',
                      help='Format of sequences, default is fasta')
    parser.add_option('-p', '--proc', dest='proc', type='int', default=1,
                      help='Specify # of processors to use')
    (options, args) = parser.parse_args()
    
    if not args:
        print "No files specified"
        return
    
    multithread(args, options.format, options.out, options.threshold, options.proc)
    
if __name__ == "__main__":
    main()

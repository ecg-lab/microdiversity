#!/usr/bin/env python
"""Script to run MUSCLE and Mothur analysis on FASTA files

MUSCLE can be optionally run prior to mothur if files aren't already aligned.
Mothur pipeline includes unique-ing, filtering, calculating distances, and
clustering of genes based on those distances.
"""

import sys, re, os, time

from multiprocessing import Pool, Manager, TimeoutError
from subprocess import check_call, list2cmdline
from optparse import OptionParser

MOTHUR_EXE = 'mothur'

BASENAME = re.compile(r"(.*)\.(.*)")

def call_and_print(cmd):
    """Calls a command and waits for it to return
    
    Raises an exception if return code was not 0
    """
    try:
        output = check_call(cmd)
    except KeyboardInterrupt:
        return
    except Exception as e:
        output = e
        
    return output

def __multi_call(args):
    
    cmd, q = args
    
    output = call_and_print(cmd)
    q.put(output)
    return output

def thread_cmds(cmds, threads=1):
    """Will multi-thread commands if threads > 1. Otherwise, normal for loop.
    """
    if not (cmds and len(cmds)):
        print "Must specify commands"
        return None
    
    if not (threads or type(threads) == int):
        print "Invalid thread count, set to 1"
        threads = 1
    
    if threads < 1:
        print "Will use all available CPUs"
        threads = None
    
    if threads == 1:
        for cmd in cmds:
            try:
                if call_and_print(cmd) != 0:
                    print "Found non-zero return value for MUSCLE"
                    break
            except KeyboardInterrupt:
                break
            except Exception as e:
                print e
                break
        
    else:
        p = Pool(threads)
        m = Manager()
        q = m.Queue()
        args = [(cmd, q) for cmd in cmds]
        muscle = p.map_async(__multi_call, args)
        i = 1
        while True:
            try:
                if muscle.ready():
                    break
                else:
                    if i % 60 == 0:
                        print "There are %i files left to align" % (len(cmds) - q.qsize())
                        sys.stdout.flush()
                    i += 1
                    time.sleep(1)
            except KeyboardInterrupt:
                p.terminate()
                break
            except Exception as e:
                print e
                p.terminate()
                break
        
        p.close()

def muscle_cmd(filename, maxiters=2, diags=False):
    """Builds MUSCLE command line with options for Popen
    """
    try:
        (name, ext) = BASENAME.match(filename).groups()
        aligned_file = "%s.aligned.%s" % (name, ext)
        
        cmd = ['muscle', '-in', filename, '-out', aligned_file,
               '-maxiters', str(maxiters), '-quiet']
               
        if diags:
            cmd.append('-diags')
        
        return cmd, aligned_file
    except Exception as e:
        print e
        return None

def test_files(files):
    """Returns the files that exist (might not be open-able though)
    """
    if files and type(files) == list and len(files) > 0:
        return filter(os.path.isfile, files)
    else:
        return []

def mothur_code(filename, cutoff=-1, prec=1000):
    """Builds mothur code for batch file
    """
    (name, ext) = BASENAME.match(filename).groups()
    
    #unique = "%s.unique.%s" % (name, ext)
    #filter = "%s.unique.filter.fasta" % name
    filter = "%s.filter.fasta" % name
    #dist = "%s.unique.filter.phylip.dist" % name
    dist = "%s.filter.phylip.dist" % name
    #names = "%s.names" % name
    
    #unique_seqs = "unique.seqs(fasta=%s);" % filename
    
    filter_seqs = "filter.seqs(fasta=%s);" % filename
    if cutoff >= 0 and cutoff < 1:
        dist_seqs = "dist.seqs(fasta=%s, output=lt, countends=F, cutoff=%g);" % (filter, cutoff)
        cluster = "cluster(phylip=%s, cutoff=%g, precision=%i, method=furthest);" % (dist, cutoff, prec)
    else:
        dist_seqs = "dist.seqs(fasta=%s, output=lt, countends=F);" % (filter)
        cluster = "cluster(phylip=%s, precision=%i, method=furthest);" % (dist, prec)
    
    #return '\n'.join([unique_seqs, filter_seqs, dist_seqs, cluster])
    return '\n'.join([filter_seqs, dist_seqs, cluster])

def mothur_code_all(filename, cutoff=-1, prec=1000):
    """Builds mothur code for batch file without uniques or filtering
    """
    (name, ext) = BASENAME.match(filename).groups()
    
    dist = "%s.phylip.dist" % name
    
    if cutoff >= 0 and cutoff < 1:
        dist_seqs = "dist.seqs(fasta=%s, output=lt, countends=F, cutoff=%g);" % (filename, cutoff)
        cluster = "cluster(phylip=%s, precision=%i, method=furthest);" % (dist, prec)
    else:
        dist_seqs = "dist.seqs(fasta=%s, output=lt, countends=F);" % (filename)
        cluster = "cluster(phylip=%s, precision=%i, method=furthest);" % (dist, prec)
    
    return '\n'.join([dist_seqs, cluster])
    
    
def make_mothur_batch_file(files, cutoff=-1, prec=1000, filtering=False, batchfile='mothur_batch.txt'):
    """Build file to run mothur in batch mode
    """
    
    if not (files and len(files) > 0):
        print "No files specified"
        return
    
    if prec < 0:
        prec = 1000
        
    if not batchfile:
        batchfile = 'mothur_batch.txt'
    
    with open(batchfile, 'wb') as w:
        for file in files:
            if filtering:
                w.write(mothur_code(file, cutoff, prec)+'\n')
            else:
                w.write(mothur_code_all(file, cutoff, prec)+'\n')
            
    return batchfile
    
def main():
    parser = OptionParser(usage="Usage: %prog [options] file[...]")
    
    parser.add_option("-c", "--cutoff", dest="cutoff", type="float", default=-1,
                      help="specify cutoff for clustering")
    parser.add_option("-p", "--prec", dest="prec", type="int", default=1000,
                      help="Specify precision for clustering to use (default=1000)")
    parser.add_option("-b", "--batch", dest="batch", action="store_true", default=False,
                      help="If flagged, will write batch file, but not run mothur")
    parser.add_option("-f", "--filter", dest="filter", action="store_true", default=False,
                      help="If flagged, will perform filtering step in mothur")
    parser.add_option("-m", "--muscle", dest="muscle", action='store_true', default=False,
                      help="If flagged, with perform alignment with muscle first")
    parser.add_option("-t", "--threads", dest="threads", type="int", default=1,
                      help="Specify number of threads to spawn MUSCLE instances (default=1)")
    parser.add_option("-i", "--maxiters", dest="maxiters", type="int", default=2,
                      help="Specifies maxiters for muscle alignment (default=2)")
    parser.add_option("-d", "--diags", dest="diags", action='store_true', default=False,
                      help="If flagged, will find diagonals for muscle alignment")
    
    (options, args) = parser.parse_args()
    
    if not (args and len(args) > 0):
        print "No files specified"
        return
    
    threads = options.threads
    cutoff = options.cutoff
    prec = options.prec
    
    if options.muscle:
        maxiters = options.maxiters
        diags = options.diags
        temp = [muscle_cmd(arg, maxiters, diags) for arg in args]
        cmds = [x[0] for x in temp]
        args = [x[1] for x in temp]
        thread_cmds(cmds, threads)
        
        args = test_files(args)
    
    batchfile = make_mothur_batch_file(args, cutoff, prec, filtering=options.filter)
    
    print "Wrote mothur batchfile to %s" % batchfile
    
    if not options.batch:
        call_and_print([MOTHUR_EXE, batchfile])
            
if __name__ == "__main__":
    main()
    


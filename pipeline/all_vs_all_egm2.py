#!/usr/bin/env python
"""Utility to compare sequences files using Afree & EGM2
"""

import os
import sys
import time

from optparse import OptionParser
from subprocess import Popen, list2cmdline, PIPE
from itertools import combinations, product

import multiprocessing

# YOU WILL NEED TO CHANGE THESE
AFREE_EXE = 'afree'
EGM2_EXE = 'egm2'

def all_vs_all_pairs(files):
    """Return all pairs of comparisons
    """
    return combinations(files, r=2)

def one_vs_all_pairs(one, files):
    """Return one vs all pairs
    """
    return [(file, one) for file in files]

def self_pairs(files):
    """Return pairs of self comparisons
    """
    return [(file, file) for file in files]
    
def dir_vs_dir_pairs(dir1, dir2):
    """Returns pairs comparing one directory to another
    """
    files1 = [os.path.join(dir1,f) for f in os.listdir(dir1) if os.path.isfile(os.path.join(dir1,f))]
    files2 = [os.path.join(dir2,f) for f in os.listdir(dir2) if os.path.isfile(os.path.join(dir2,f))]
    return product(files1, files2)
    
def list_cmds(pairs, c=40, k=5, hi=4, sd=10, afree=True):
    """Returns Afree or EGM2 commands with given arguments
    """
    if c > 100: c = 100
    if c < 0: c = 0
    if k < 1: k = 1
    if hi < 1: hi = 1
    if sd < 0: sd = 0
    
    c = str(c)
    k = str(k)
    hi = str(hi)
    sd = str(sd)
    
    all_cmds = []
    if afree:
        for o, pair in enumerate(pairs):
            all_cmds.append(get_afree_cmd(pair[0], pair[1], str(o), sd, k))
    else:
        for o, pair in enumerate(pairs):
            all_cmds.append(get_egm_cmd(pair[0], pair[1], str(o), c, k, hi, sd))
        
    return all_cmds

def get_egm_cmd(file1, file2, o, c='40', k='5', hi='4', sd='15'):
    """Returns EGM command with given input
    
    Suitable for input to Popen
    """
    if not (file1 and file2 and o):
        return None
    
    try:
        with open(file1, 'rb'):
            pass
        with open(file2, 'rb'):
            pass
    except Exception as e:
        print e
        return None
    
    cmd = [EGM2_EXE,
           '-g1', file1,
           '-g2', file2,
           '-o', o,
           '-c', c,
           '-k', k,
           '-hi', hi,
           '-sd', sd,
          ]
    
    return cmd

def get_afree_cmd(file1, file2, o, c='10', k='5'):
    """Returns Afree command with given input
    
    Suitable for input to Popen
    """
    if not(file1 and file2 and o):
        return None
        
    try:
        with open(file1, 'rb'):
            pass
        with open(file2, 'rb'):
            pass
    except Exception as e:
        print e
        return None
    
    cmd = [AFREE_EXE,
           '-g1', file1,
           '-g2', file2,
           '-o', o,
           '-c', c,
           '-k', k,
          ]

    return cmd

def call_and_print(cmd):
    """Call a command with Popen and print its output
    """
    try:
        proc = Popen(cmd, stdout=PIPE, stderr=PIPE)
        (output, error) = proc.communicate()
        if output:
            return "%s\n%s" % (list2cmdline(cmd), output)
        elif error:
            return error
    except KeyboardInterrupt:
        if proc: proc.terminate()
    except Exception as e:
        print e
        if proc: proc.terminate()

def single_thread_cmds(cmds):
    """Call each command one by one
    """
    for cmd in cmds:
        try:
            print call_and_print(cmd)
            sys.stdout.flush()
        except KeyboardInterrupt:
            break
        except Exception as e:
            print e
            break

  
def __multi_call(args):
    """Wrapper for call_and_print with a queue
    """
    cmd, q = args
    
    output = call_and_print(cmd)
    q.put(output)
    return output
    
def thread_cmds(cmds, threads=1):
    """Multi-threaded call commands
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
        print "Single thread specified. Will not multi-thread."
        single_thread_cmds(cmds)
    else:
        p = multiprocessing.Pool(threads)
        m = multiprocessing.Manager()
        q = m.Queue()
        args = [(cmd, q) for cmd in cmds]
        egm2 = p.map_async(__multi_call, args)
        i = 0
        results = None
        while True:
            try:
                if egm2.ready():
                    results = egm2.get()
                    p.close()
                    break
                else:
                    if i % 60 == 0:
                        print "There are %i remaining comparisons" % (len(cmds) - q.qsize())
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
    
        if results:
            for result in results:
                print result
                sys.stdout.flush()
    
    

def main():
    
    parser = OptionParser(usage="Usage: %prog [options] file [...]")
    parser.add_option('-1', '--vs1', dest='one', default=None,
                      help='If specified, will run EGM of all files versus specified file')
    parser.add_option('--dir1', dest='dir1', default=None,
                      help='With dir2, specify directory of files to compare')
    parser.add_option('--dir2', dest='dir2', default=None,
                      help='With dir1, specify other directory of files to compare')
    parser.add_option('--self', dest='self', action='store_true', default=False,
                      help='If flagged, will run egm/afree on each file vs itself')
    parser.add_option("-t", "--threads", dest="threads", type="int", default=1,
                      help="Specify number of threads to run")
    parser.add_option("-c", "--c", dest="c", type="int", default=40,
                      help="minimum protein coverage threshold: default c=40")
    parser.add_option("-k", "--k", dest="k", type="int", default=5,
                      help="k-mer length:default k=5")
    parser.add_option("-i", "--hi", dest="hi", type="int", default=4,
                      help="gene mapping iteration: default hi=4")
    parser.add_option("-s", "--sd", dest="sd", type="int", default=10,
                      help="Sorensen-Dice similarity index: default sd=10")
    parser.add_option("-a", "--afree", dest='afree', action='store_true', default=False,
                      help='Not used anymore. Still here for posterity.')
    parser.add_option('-e', '--egm', dest='egm', action='store_true', default=False,
                      help='Use EGM2, instead of Afree')
    parser.add_option("-d", "--debug", dest="debug", action="store_true", default=False,
                      help='Print commands. Do not run them.')
    (options, args) = parser.parse_args()

    if options.dir1 and options.dir2:
        pairs = dir_vs_dir_pairs(options.dir1, options.dir2)
    else:
        if not args:
            print 'No files specified'
            return
        if options.one:
            pairs = one_vs_all_pairs(options.one, args)
        elif options.self:
            pairs = self_pairs(args)
        else:
            pairs = all_vs_all_pairs(args)
    
    if not pairs:
        print 'No pairs found'
        return
    
    c = options.c
    k = options.k
    hi = options.hi
    sd = options.sd
    
    cmds = list_cmds(pairs, c, k, hi, sd, not options.egm)
    
    if options.debug:
        # print cmds for debugging
        for cmd in cmds: print list2cmdline(cmd)
        return
    
    thread_cmds(cmds, options.threads)
    
if __name__ == "__main__":
    main()

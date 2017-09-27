#!/usr/bin/env python
"""This script will create phylogenetic trees from a coalescent simulation
"""

import sys
import imp
try:
    imp.find_module('Bio')
except ImportError:
    print "No Bio Module. Wrong Python verison?"
    sys.exit(1)

import random
import time
import multiprocessing
import Bio.Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, _DistanceMatrix
from math import ceil
from itertools import combinations, chain
from optparse import OptionParser

def read_ntaxa_file(file):
    """Read file where size of tree defined by integers
    
    Can read white space delimited or one per line
    """
    if not file:
        return
    
    with open(file, 'rb') as f:
        ntaxa = []
        for line in f:
            try:
                ntaxa.extend([int(x) for x in line.split()])
            except Exception as e:
                print e
    
        return ntaxa

def seed(n=100):
    """Returns initial seed of n species
    """
    if n < 1:
        print "n < 1 for seed()"
        return []
    return [[i] for i in range(n)]
    
def next_generation(cur_gen=[], p_nothing=0.0):
    """Simulate one generation based on cur_gen
    """
    if not (cur_gen and type(cur_gen) == list and len(cur_gen) > 0):
        print "Cur_gen does not have taxa"
        return []

    if p_nothing > 0 and random.random() < p_nothing:
        return cur_gen
    
    max_int = len(cur_gen)-1
    die_off = random.randint(0, max_int)
    temp = die_off
    # make sure the die off and duplicate taxa are different
    while temp == die_off:
        temp = random.randint(0, max_int)
    duplicate = temp

    next_gen = []
    for i, item in enumerate(cur_gen):
        if i == die_off:
            # taxon dies, don't add to next_gen
            continue
        if i == duplicate:
            temp = item[:] # force copy of item
            temp.append(len(next_gen))
            next_gen.append(temp)

        item.append(len(next_gen))
        next_gen.append(item)

    return next_gen

def test_coalesce(gen=[]):
    """Test to see if gen has coalesced
    """
    if not (gen and type(gen) == list and len(gen) > 0):
        print "No generation to test for coalescence"
        return False
    test = gen[0][0]
    for taxon in gen:
        if taxon[0] != test:
            return False

    return True

def get_latest_coalescence(gen=[]):
    """Return the latest coalescence of a generation

    This will remove generational data from before the latest coalescence
    """
    if not (gen and type(gen) == list and len(gen) > 0):
        print "Gen was empty or not a list"
        return []

    consensus = gen[0]
    last = 0
    for x in gen:
        for i, y in enumerate(x):
            if i >= len(consensus):
                break
            if y != consensus[i]:
                consensus = consensus[:i]
                last = i
                break

    return [x[i-1:] for x in gen]



def simulate(n_taxa=100, n_gen=10000, p_nothing=0.0):
    """Simulate coalescence on n_taxa, checking for coalescence every n_gen
    generations.

    Can also sample from taxa (simulate incomplete sampling).
    """
    if not (n_taxa or type(n_taxa) == int or n_taxa > 0):
        print "n_taxa had invalid value"
        return []
    if not (n_gen or type(n_gen) == int or n_gen > 0):
        print "n_gen had invalid value"
        return []

    cur_gen = seed(n_taxa)
    while not test_coalesce(cur_gen):
        for i in xrange(n_gen):
            try:
                cur_gen = next_generation(cur_gen, p_nothing)
            except KeyboardInterrupt:
                return []
            except Exception as e:
                print e
                break

    return get_latest_coalescence(cur_gen)

def sample_gen(gen=[], n=None):
    """Sample n taxa from gen
    """
    if not (gen and type(gen) == list and len(gen) > 0):
        print "Gen was not a list or was empty"
        return []

    if n and type(n) == int and n < len(gen):
        return random.sample(gen, n)
    else:
        return gen
    
def make_samples(gen=[], n_samples=None, n_taxa=None):
    """Sub-sample n times from gen
    """
    
    if not (gen and type(gen) == list and len(gen) > 0):
        print "gen had invalid value"
        return []
    if not (n_taxa and type(n_taxa) == int and n_taxa > 0):
        print "n_taxa had invalid value"
        return gen
    if not (n_samples and type(n_samples) == int and n_samples > 0):
        print "n_samples had invalid value"
        return gen
    
    samples = []
    for i in xrange(n_samples):
        samples.append(sample_gen(gen, n_taxa))
    
    return samples
    
def get_dist(x=[], y=[]):
    """Get the distance (in generations) between x and y
    """
    if not (x and type(x) == list and len(x) > 0):
        print "x was invalid entry"
        return None
    if not (y and type(y) == list and len(y) > 0):
        print "y was invalid entry"
        return None

    if len(x) != len(y):
        print "Lengths of x and y were different. Will use smaller of two."

    if len(x) > len(y):
        length = len(y)
    else:
        length = len(x)

    for i in range(length):
        if x[i] != y[i]:
            return (length - i)*2

    return 0

def dist_matrix(gen=[], normalize=False, lt=True):
    """Build a lower triangle distance matrix of taxa in gen
    """
    if not (gen and type(gen) == list and len(gen) > 0):
        print "gen was empty or not a list"
        return [[]]

    gen_length = len(gen)
    indices = list(combinations(range(gen_length), 2))
    dists = []
    for i, (x, y) in enumerate(combinations(gen, 2)):
        dists.append(get_dist(x, y))

    if normalize:
        max_dist = float(max(dists))

    matrix = []
    for j in range(gen_length):
        temp = []
        for i in range(gen_length):
            if i == j:
                temp.append(0)
                continue
            if not lt and (j, i) in indices:
                x = j
                j = i
                i = x
            if (i, j) in indices:
                if normalize:
                    temp.append(dists[indices.index((i,j))]/max_dist)
                else:
                    temp.append(dists[indices.index((i,j))])
        matrix.append(temp)

    return matrix

def multiply_matrix(matrix, value):
    """Multiply all values in a matrix by the specified value
    """
    
    return [[x*value for x in row] for row in matrix]

def write_matrix(matrix, file='infile'):
    """Writes the matrix to a formatted file

    This file can be read in the program neighbor
    """
    if not (matrix and type(matrix) == list and len(matrix) > 0):
        print "Empty or invalid matrix value"
        return
    if not (file and type(file) == str):
        print "Invalid or empty file string"
        return

    with open(file, 'wb') as w:
        w.write('\t%i\n'%len(matrix))
        for i, line in enumerate(matrix):
            w.write('X%08i\t'%i)
            w.write('\t'.join(["%f" % x for x in line])+'\n')

def construct_tree(matrix, nj=True):
    """Build a tree from a distance matrix

    Can either use neighbor-joining (nj) or UPGMA.
    """

    if not (matrix and type(matrix) == list and len(matrix) > 0):
        print "matrix has invalid value"
        return

    dm = _DistanceMatrix(names=[str(i) for i in range(len(matrix))], matrix=matrix)
    
    constructor = DistanceTreeConstructor()
    if nj:
        tree = constructor.nj(dm)
    else:
        tree = constructor.upgma(dm)
    
    # this will remove the names from the inner nodes
    # this is critical for seq-gen to read in the tree
    for clade in tree.get_nonterminals():
        clade.name = ''
    
    return tree

def write_trees(trees, file, format='newick', individual=False):
    """Write trees to a file in a given format
    """
    if not (trees and (type(trees) == list or type(trees) == Bio.Phylo.BaseTree.Tree)):
        print "No trees found"
        return
    if not (file and type(file) == str):
        print "No file specified"
        return
    if not (format and type(format) == str):
        print "Invalid format. Will use newick"
        format = 'newick'

    if individual:
        for i, tree in enumerate(trees):
            Bio.Phylo.write(tree, '%i.%s' % (i, file), format)
    else:
        Bio.Phylo.write(trees, file, format)       

def run_simulation(n_sim=1, n_taxa=100, n_gen=10000, p_nothing=0.0, samples=None, frac=None, nj=True, dist=1.0, out='outtree', outgroup=False, individual=False, interactive=False, matrix=False):
    """Run n_sim simulations with n_taxa, checking every n_gen. Output to out.

    Can also sample. Can be run interactively (will ask about errors).
    """

    if not (n_sim and type(n_sim) == int and n_sim > 0):
        print "n_sim had invalid value"
        return
    if not (n_taxa and type(n_taxa) == int and n_taxa > 0):
        print "n_taxa had invalid value"
        return
    if not (n_gen and type(n_gen) == int and n_gen > 0):
        print "n_gen had invalid value"
        return
    if p_nothing and not (type(p_nothing) == float and p_nothing >= 0 and p_nothing < 1):
        print "p_nothing had invalid value. Will set to zero."
        p_nothing = 0.0
    if samples and not (type(samples) == int and samples > 0):
        print "sample had invalid value. Will not sample."
        samples = None
    if frac and not(type(frac) == float and frac > 0):
        print "frac had invalid value. Will not sample."
        samples = None
    if not (dist and (type(dist) == float or type(dist) == int)):
        print "invalid normalized distance"
        dist = 1.0
    if not (out and type(out) == str):
        print "invalid output filename"
        return
    
    i = 0
    trees = []
    while i < n_sim:
        try:
            sim = simulate(n_taxa, n_gen, p_nothing)
            if samples and frac:
                sampled_taxa = int(round(frac*n_taxa))
                if sampled_taxa < 1:
                    print "Frac set too small for # of taxa."
                    return
                sampled = make_samples(sim, samples, sampled_taxa)
                for j, sample in enumerate(sampled):
                    dm = dist_matrix(sample)
                    if outgroup:
                        tempdist = max(list(chain.from_iterable(dm)))
                        dm.append([tempdist]*(len(dm)+1))
                    if dist:
                        norm_value = 1. / max(list(chain.from_iterable(dm))) * dist
                        dm = multiply_matrix(dm, norm_value)
                    if matrix:
                        write_matrix(dm, '%s_%i_%i.dist' % (out, i, j))
                    # name = 'infile_%i_%i' % (i, j)
                    # write_matrix(dm, name)
                    # matrices.append(name)
                    tree = construct_tree(dm, nj)
                    trees.append(tree)
                
                dm = dist_matrix(sim)
                if dist:
                    norm_value = 1. / max(chain.from_iterable(dm)) * dist
                    dm = multiply_matrix(dm, norm_value)
                if matrix:
                    write_matrix(dm, 'true.dist')
                # if nj:
                    # dist_exe = '/home/tstraub/bin/neighbor'
                # else:
                    # dist_exe = '/home/tstraub/bin/kitsch/'
                # make_trees.process_files(files=['true.dist'], out='true.tre', distance=dist_exe)
                true = construct_tree(dm, nj)
                write_trees(true, 'true.%s'%out)
                print "Wrote 'true' tree to true.%s"%out
            else:
                dm = dist_matrix(sim)
                if outgroup:
                    tempdist = max(list(chain.from_iterable(dm)))
                    dm.append([tempdist]*(len(dm)+1))
                if dist:
                    norm_value = 1. / max(chain.from_iterable(dm)) * dist
                    dm = multiply_matrix(dm, norm_value)
                if matrix:
                    write_matrix(dm, '%s_%i.dist' % (out, i))
                # name = 'infile_%i' % i
                # write_matrix(dm, name)
                # matrices.append(name)
                tree = construct_tree(dm, nj)
                trees.append(tree)
            print "Simulation %i of %i complete." % (i + 1, n_sim)
        except KeyboardInterrupt:
            return
        except Exception as e:
            print e
            if interactive:
                s = raw_input('(R)etry? (A)bort? (C)ontinue?')
                if s[0] == 'R' or s[0] == 'r':
                    continue
                elif s[0] == 'A' or s[0] == 'a':
                    break
                elif s[0] == 'C' or s[0] == 'c':
                    pass
                else:
                    print "%s is not an option. Aborting" % s
                    break
            else:
                break

        i += 1

    if trees and len(trees) > 0:
        write_trees(trees, out, individual=individual)
        # make_trees.process_files(matrices, out)
        
        print "Wrote %i trees to %s" % (len(trees), out)
    else:
        print "No trees simulated"

def __sim_wrapper(args):
    """Multi-thread safe simulation.
    
    Returns tuple of (tree, matrix, sampled trees, sampled matrices)
    """
    
    (n_taxa, n_gen, p_nothing, samples, frac, nj, dist, outgroup, q) = args
    try:
        sim = simulate(n_taxa, n_gen, p_nothing)
        if samples and frac:
            trees = []
            matrices = []
            sampled_taxa = int(round(frac*n_taxa))
            if sampled_taxa < 1:
                return
            sampled = make_samples(sim, samples, sampled_taxa)
            for j, sample in enumerate(sampled):
                dm = dist_matrix(sample)
                if dist:
                    norm_value = 1. / max(list(chain.from_iterable(dm))) * dist
                    dm = multiply_matrix(dm, norm_value)
                matrices.append(dm)
                tree = construct_tree(dm, nj)
                trees.append(tree)
            dm = dist_matrix(sim)
            if outgroup:
                tempdist = max(list(chain.from_iterable(dm)))
                dm.append([tempdist]*(len(dm)+1))
            if dist:
                norm_value = 1. / max(chain.from_iterable(dm)) * dist
                dm = multiply_matrix(dm, norm_value)
            true = construct_tree(dm, nj)
            q.put((true, dm, trees, matrices))
            return (true, dm, trees, matrices)
        else:
            dm = dist_matrix(sim)
            if outgroup:
                tempdist = max(list(chain.from_iterable(dm)))
                dm.append([tempdist]*(len(dm)+1))
            if dist:
                norm_value = 1. / max(chain.from_iterable(dm)) * dist
                dm = multiply_matrix(dm, norm_value)
            tree = construct_tree(dm, nj)
            q.put((tree, dm, [], []))
            return (tree, dm, [], [])
    except KeyboardInterrupt:
        return
    except Exception as e:
        print e
        return
    

def multi_thread_sim(n_sim=1, n_taxa=100, n_gen=10000, p_nothing=0.0, samples=None, frac=None, nj=True, dist=1.0, out='outtree', outgroup=False, individual=False, interactive=False, matrix=False, threads=1):
    """Multithread simulation. Run n_sim simulations with n_taxa, checking every n_gen. Output to out.

    Can also sample. Can be run interactively (will ask about errors).
    """

    if not (n_sim and type(n_sim) == int and n_sim > 0):
        print "n_sim had invalid value"
        return
    if not (n_taxa and type(n_taxa) == int and n_taxa > 0):
        print "n_taxa had invalid value"
        return
    if not (n_gen and type(n_gen) == int and n_gen > 0):
        print "n_gen had invalid value"
        return
    if p_nothing and not (type(p_nothing) == float and p_nothing >= 0 and p_nothing < 1):
        print "p_nothing had invalid value. Will set to zero."
        p_nothing = 0.0
    if samples and not (type(samples) == int and samples > 0):
        print "sample had invalid value. Will not sample."
        samples = None
    if frac and not(type(frac) == float and frac > 0):
        print "frac had invalid value. Will not sample."
        samples = None
    if not (dist and (type(dist) == float or type(dist) == int)):
        print "invalid normalized distance"
        dist = 1.0
    if not (out and type(out) == str):
        print "invalid output filename"
        return
    
    if threads > 1:
        p = multiprocessing.Pool(threads)
        m = multiprocessing.Manager()
        q = m.Queue()
        args = [(n_taxa, n_gen, p_nothing, samples, frac, nj, dist, outgroup, q) for i in xrange(n_sim)]
        exp = p.map_async(__sim_wrapper, args)
        iter = 1
        results = None
        while True:
            try:
                if exp.ready():
                    results = exp.get()
                    break
                else:
                    if iter % 60 == 0:
                        print "There are %i trees left to simulate" % (n_sim - q.qsize())
                        sys.stdout.flush()
                    iter += 1
                    time.sleep(1)
            except KeyboardInterrupt:
                p.terminate()
                break
            except Exception as e:
                print e
                p.terminate()
                break
        
        p.close()
        
        if results and len(results) > 0:
            if samples and frac:
                for i, result in enumerate(results):
                    if result[0] and result[1] and result[2] and result[3]:
                        write_tree(result[0], '%s_true_%i.tre' % (out, i))
                        write_trees(result[2], '%s_sampled_%i.tre' % (out, i), individual)
                        if matrix:
                            write_matrix(result[1], '%s_true_%i.dist' % (out, i))
                            for j, dm in enumerate(result[3]):
                                write_matrix(dm, '%s_sampled_%i_%j.dist' % (i, j))
                    
            else:
                trees = []
                matrices = []
                for result in results:
                    if result:
                        if result[0]:
                            trees.append(result[0])
                        if result[1]:
                            matrices.append(result[1])
                if trees and len(trees) > 0:
                    write_trees(trees, out, individual=individual)
                if matrix and matrices and len(matrices) > 0:
                    for i, dm in enumerate(matrices):
                        write_matrix(dm, '%s_%i.dist' % (out, i))
    else:
        run_simulation(n_sim, n_taxa, n_gen, p_nothing, samples, frac, nj, dist, out, individual, interactive, matrix)

def _multi_thread_sim_ntaxa(n_taxa=[], n_gen=10000, p_nothing=0.0, samples=None, frac=None, nj=True, dist=1.0, out='outtree', outgroup=False, individual=False, interactive=False, matrix=False, threads=1):
    """Multithread simulation. Run simulations based on n_taxa, checking every n_gen. Output to out.

    Can also sample. Can be run interactively (will ask about errors).
    """

    if not (n_taxa and type(n_taxa) == list and len(n_taxa) > 0):
        print "n_taxa had invalid value"
        return
    if not (n_gen and type(n_gen) == int and n_gen > 0):
        print "n_gen had invalid value"
        return
    if p_nothing and not (type(p_nothing) == float and p_nothing >= 0 and p_nothing < 1):
        print "p_nothing had invalid value. Will set to zero."
        p_nothing = 0.0
    if samples and not (type(samples) == int and samples > 0):
        print "sample had invalid value. Will not sample."
        samples = None
    if frac and not(type(frac) == float and frac > 0):
        print "frac had invalid value. Will not sample."
        samples = None
    if not (dist and (type(dist) == float or type(dist) == int)):
        print "invalid normalized distance"
        dist = 1.0
    if not (out and type(out) == str):
        print "invalid output filename"
        return
    
    p = multiprocessing.Pool(threads)
    m = multiprocessing.Manager()
    q = m.Queue()
    args = [(n, n_gen, p_nothing, samples, frac, nj, dist, outgroup, q) for n in n_taxa]
    exp = p.map_async(__sim_wrapper, args)
    iter = 1
    results = None
    while True:
        try:
            if exp.ready():
                results = exp.get()
                break
            else:
                if iter % 60 == 0:
                    print "There are %i trees left to simulate" % (len(n_taxa) - q.qsize())
                    sys.stdout.flush()
                iter += 1
                time.sleep(1)
        except KeyboardInterrupt:
            p.terminate()
            break
        except Exception as e:
            print e
            p.terminate()
            break
    
    p.close()
    
    if results and len(results) > 0:
        if samples and frac:
            for i, result in enumerate(results):
                if result[0] and result[1] and result[2] and result[3]:
                    write_tree(result[0], '%s_true_%i.tre' % (out, i))
                    write_trees(result[2], '%s_sampled_%i.tre' % (out, i), individual)
                    if matrix:
                        write_matrix(result[1], '%s_true_%i.dist' % (out, i))
                        for j, dm in enumerate(result[3]):
                            write_matrix(dm, '%s_sampled_%i_%j.dist' % (i, j))
                
        else:
            trees = []
            matrices = []
            for result in results:
                if result:
                    if result[0]:
                        trees.append(result[0])
                    if result[1]:
                        matrices.append(result[1])
            if trees and len(trees) > 0:
                write_trees(trees, out, individual=individual)
            if matrix and matrices and len(matrices) > 0:
                for i, dm in enumerate(matrices):
                    write_matrix(dm, '%s_%i.dist' % (out, i))

def _experiment(taxa=100, gen=10000, p_nothing=0.0, sampling=[.1,.25,.5,.75,.9], n_per_sampling=100, nj=False, out='out.tre'):
    """Run a sampling experiment
    """

    sim = simulate(taxa, gen, p_nothing)
    dm = dist_matrix(sim, False)
    
    norm_value = 1. / max(list(chain.from_iterable(dm)))
    
    dm = multiply_matrix(dm, norm_value)
    
    true_tree = construct_tree(dm, nj)
    write_trees(true_tree, 'true.%s'%out)
    
    # write_matrix(dm, 'true.dist')
    # if nj:
        # dist_exe = '/home/tstraub/bin/neighbor'
    # else:
        # dist_exe = '/home/tstraub/bin/kitsch/'
    # make_trees.process_files(files=['true.dist'], out='true.tre', distance=dist_exe)
    
    for i, frac in enumerate(sampling):
        trees = []
        sampled_taxa = int(round(frac*taxa))
        sampled = make_samples(sim, n_per_sampling, sampled_taxa)
        for sample in sampled:
            dm = dist_matrix(sample, False)
            dm = multiply_matrix(dm, norm_value)
            # name = 'infile_%i' % i
            # write_matrix(dm, name)
            # matrices.append(name)
            tree = construct_tree(dm, nj)
            trees.append(tree)
        if matrices and len(matrices) > 0:
            write_trees(trees, 's%g.%s'%(frac,out))
            # return make_trees.process_files(files=matrices, out=out, distance=dist_exe)

def __exp_wrapper(arguments):
    return _experiment(*arguments)
    
def _experiment2(n=100, taxa=100, gen=10000, sampling=[.25,.5,.75], n_per_sampling=100, nj=False, out='tree', threads=1):
    """Run another sampling experiment
    """
    
    if threads > 1:
        p = multiprocessing.Pool(threads)
        args = [(taxa, gen, 0, sampling, n_per_sampling, nj, "%s.%i"%(out, i)) for i in range(n)]
        exp = p.map_async(__exp_wrapper, args)
        while not exp.ready():
            try:
                results = exp.get(5)
                #if results:
                    #for result in results:
                        #print result
                        #sys.stdout.flush()
            except multiprocessing.TimeoutError:
                pass
            except KeyboardInterrupt:
                p.terminate()
                break
            except Exception as e:
                print e
                p.terminate()
                break
    else:
        for i in range(n):
            print "Running trial %i of %i" % (i, n)
            experiment(taxa=taxa, gen=gen, sampling=sampling, n_per_sampling=n_per_sampling, nj=nj, out='%s.%i'%(out,i))

def __exp3(arguments):
    """Multi-thread safe simulation for experiment3
    """
    
    
    try:
        (taxa, gen, sampling, nj, dist) = arguments
        
        gen = sample_gen(simulate(int(taxa / sampling), gen), taxa)
        if gen:
            dm = dist_matrix(gen)
            if dm:
                if dist:
                    norm_value = 1. / max(chain.from_iterable(dm)) * dist
                    dm = multiply_matrix(dm, norm_value)
                return construct_tree(dm, nj)
    except KeyboardInterrupt:
        return
    except Exception as e:
        print e
        return
    

def _experiment3(n=100, taxa=100, gen=10000, sampling=[1.,.9,.8,.7,.6,.5,.4,.3], nj=True, dist=1.0, out='out', threads=1):
    """Generate unrelated trees of given # of taxa
    
    But with different levels of sampling
    """
    
    if threads > n:
        print "Will only use %i threads. (More does nothing)" % n
        threads = n
    
    for sample in sampling:
        print "Starting sampling of %g" % sample
        if threads > 1:
            p = multiprocessing.Pool(threads)
            args = [(taxa, gen, sample, nj, dist) for i in xrange(n)]
            exp = p.map_async(__exp3, args)
            while not exp.ready():
                try:
                    trees = exp.get(5)
                    if len(trees):
                        trees = filter(None, trees)
                        write_trees(trees, '%s_%g.tre' % (out, sample))
                        print "Wrote %i trees" % len(trees)
                except multiprocessing.TimeoutError:
                    pass
                except KeyboardInterrupt:
                    p.terminate()
                    break
                except Exception as e:
                    print e
                    p.terminate()
                    break
        else:
            trees = []
            for i in xrange(n):
                try:
                    temp = __exp3((taxa, gen, sample, nj, dist))
                    if temp:
                        trees.append(temp)
                except Exception as e:
                    print e
                    break
            if len(trees):
                write_trees(trees, '%s_%g.tre' % (out, sample))
                print "Wrote %i trees" % len(trees)
    
    
    

if __name__ == "__main__":
    parser = OptionParser(usage="Usage: %prog [options] -n TREES -t TAXA -o OUTPUT")
    parser.add_option('-n', '--trees', dest='n_trees', type='int', default=1,
                      help='# of independent trees to simulate')
    parser.add_option('-t', '--taxa', dest='taxa', type='int', default=100,
                      help='# of taxa to simulate')
    parser.add_option('--taxafile', dest='taxafile', default='',
                      help='Read # of taxa from file')
    parser.add_option('-g', '--gen', dest='gen', type='int', default=10000,
                      help='Test for coalescence every # of generations')
    parser.add_option('-p', '--prob', dest='p_nothing', type='float', default=0.0,
                      help='Probability of nothing happening in a generation')
    parser.add_option('-s', '--samples', dest='samples', type='int', default=None,
                      help='If specified, will only subsample n trees from "true" tree')
    parser.add_option('-f', '--frac', dest='frac', type='float', default=None,
                      help='If specified, will sample frac x # of taxa to form sampled trees')
    parser.add_option('-u', '--upgma', dest='upgma', action='store_true', default=False,
                      help='If flagged, will use UPGMA instead of NJ')
    parser.add_option('-d', '--dist', dest='dist', type='float', default = 1.0,
                      help='Normalize distances so that max is this value')
    parser.add_option('-o', '--out', dest='out', default='out.tre',
                      help='Output tree(s)to file')
    parser.add_option('-i', '--individual', dest='individual', action='store_true', default=False,
                      help='Output each tree to a separate newick file')
    parser.add_option('--outgroup', dest='outgroup', action='store_true', default=False,
                      help='Add an outgroup to the simulated tree, so as to root it')
    parser.add_option('-m', '--matrix', dest='matrix', action='store_true', default=False,
                      help='Output distance matrix to files')
    parser.add_option('--threads', dest='threads', type='int', default=1,
                      help='Specify # of threads to use in simulation')
    (options, args) = parser.parse_args()

    if options.taxafile:
        # code if you want to seed in n_taxa values from a file
        n_taxa = read_ntaxa_file(options.taxafile)
        if len(n_taxa) < options.n_trees:
            # if you specified more trees than seed values
            n_taxa = n_taxa * (int(options.n_trees / len(n_taxa))+1)
        if len(n_taxa) != options.n_trees:
            # sample to correct number of n_taxa values for # trees to simulate
            n_taxa = random.sample(n_taxa, options.n_trees)
        
        _multi_thread_sim_ntaxa(n_taxa=n_taxa, n_gen=options.gen, p_nothing=options.p_nothing, samples=options.samples, frac=options.frac, nj=not options.upgma, dist=options.dist, out=options.out, individual=options.individual, outgroup=options.outgroup, interactive=True, matrix=options.matrix, threads=options.threads)
    else:
        multi_thread_sim(n_sim=options.n_trees, n_taxa=options.taxa, n_gen=options.gen, p_nothing=options.p_nothing, samples=options.samples, frac=options.frac, nj=not options.upgma, dist=options.dist, out=options.out, individual=options.individual, outgroup=options.outgroup, interactive=True, matrix=options.matrix, threads=options.threads)

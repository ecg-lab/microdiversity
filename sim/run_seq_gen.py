#!/usr/bin/env python
"""Utility to help run Seq-Gen on multiple files
"""

import os
import sys
import random
import numpy
import re
import StringIO

from optparse import OptionParser
from subprocess import Popen, list2cmdline, PIPE

import imp
try:
    imp.find_module('Bio')
except ImportError:
    print 'No Bio Module. Wrong Python verison?'
    sys.exit(1)

from Bio import SeqIO, Phylo, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

# YOU MAY NEED TO CHANGE THIS
SEQ_GEN_EXE = 'seq-gen'

def get_seqgen_command(treefile, length=1000, model='HKY', freq=[.25, .25, .25, .25], a=None, d=None, exe=SEQ_GEN_EXE):
    """Return seq-gen command with given parameters
    """
    try:
        with open(treefile, 'rb'):
            pass
    except Exception as e:
        print e
        return
    cmd = [exe,
           '-m%s' % model,
           '-l', '%i'%(length),
           '-q',
           treefile,
          ]
    if a:
        cmd.insert(2, '-a')
        cmd.insert(3, '%.6f'%(a))
    if d:
        cmd.insert(2, '-s')
        cmd.insert(3, '%.6f'%(d))
    if test_freq(freq):
        cmd.insert(2, '-f')
        cmd.insert(3, '%.6f' % freq[0])
        cmd.insert(4, '%.6f' % freq[1])
        cmd.insert(5, '%.6f' % freq[2])
        cmd.insert(6, '%.6f' % freq[3])
    return cmd

def list_cmds(files, lengths, model, freq, a, scales, exe=SEQ_GEN_EXE):
    """Return seq-gen commands for files with given parameters
    """
    cmds = []
    for (file, length, d) in zip(files, lengths, scales):
        cmds.append(get_seqgen_command(file, length, model, freq, a, d, exe))
    
    return cmds

def test_freq(freq=[]):
    """Make sure nt frequencies make sense
    """
    try:
        if freq and type(freq) is str:
            temp = freq.split()
            if len(temp) is 4:
                temp = [float(x) for x in temp]
                if sum(temp) == 1:
                    return temp
    except Exception:
        pass
    
    return [.25, .25, .25, .25]

def sample_from_file(file, n, replace=False):
    """Sample n numeric values from a text file with or without replacement
    """
    if not file:
        print 'No file specified'
        return
    if not n:
        print 'Number of samples unspecified'
        return
    
    dist = []
    with open(file, 'rb') as f:
        for line in f:
            try:
                dist.extend([float(x) for x in line.split()])
            except:
                pass
    
    if not len(dist):
        print "No numbers found in file"
        return
    
    
    if len(dist) > n and not replace:
        return random.sample(dist, n)
    else:
        return list(numpy.random.choice(dist, size=n, replace=True))
        
def seed_gamma(n, alpha, beta, minimum = 0.0, maximum = None):
    """Returns n float values based on a gamma distribution
    """
    values = []
    for i in xrange(n):
        temp = 0.0
        while (minimum != None and temp < minimum) or (maximum != None and temp > maximum):
            temp = random.gammavariate(alpha, beta)
        values.append(temp)
            
    
    return values

def get_cor_lengths(files, m, b, sd=.065, min=100):
    """Returns sequence lengths correlated to the size the tree
    """
    # note this won't work if there's more than 1 tree
    lengths = []
    for file in files:
        try:
            with open(file, 'rb') as f:
                temp = None
                while Temp < min:
                    temp = int((m * Phylo.read(f, 'newick').count_terminals() + b) * random.gauss(mu=1, sigma=sd))
                lengths.append(temp)
        except Exception as e:
            print e
            lengths.append(1000)
    return lengths
                
def parse_output(file, output):
    """Parses seq-gen output (phylip) into fasta file(s)
    """
    assert file
    assert output
    
    strio = StringIO.StringIO(output)
    phylip = list(AlignIO.parse(strio, 'phylip-sequential'))
    
    if len(phylip) > 1:
        for i, msa in enumerate(phylip):
            AlignIO.write(msa, '%i.%s'%(i, file), 'fasta')
    else:
        AlignIO.write(phylip[0], file, 'fasta')

def call_and_write(cmd, out):
    """Use Popen to run seq-gen command and write output to file
    """
    assert cmd
    assert out
    try:
        proc = Popen(cmd, stdout=PIPE, stderr=PIPE)
        (output, error) = proc.communicate()        
        if output:
            parse_output(out, output) # supports multiple trees
            return 'Wrote to %s.' % out
        elif error:
            return error
    except KeyboardInterrupt:
        if proc: proc.terminate()
    except Exception as e:
        print e
        if proc: proc.terminate()

def run_cmds(cmds, outputs):
    """Run each command and write to each output
    """
    for (cmd, output) in zip(cmds, outputs):
        try:
            print list2cmdline(cmd)
            print call_and_write(cmd, output)
        except KeyboardInterrupt:
            break
        except Exception as e:
            print e
            break


def get_output_names(files=[]):
    """Return output names based on input files
    """
    assert files
    
    output = []
    basename = re.compile(r'(.*)\.')
    
    for file in files:
        match = basename.match(file)
        if match and match.groups():
            output.append('%s.fa' % match.groups()[0])
        else:
            output.append('%s.fa'%file)
    
    return output
            
def main():
    parser = OptionParser(usage='Usage: %prog [options] file[...]')
    parser.add_option('-m', '--model', dest='model', default='HKY',
                      help='Specify model for Seq-Gen to use (default=HKY). At this time only HKY with parameters set to a JC model is supported.')
    parser.add_option('-f', '--freq', dest='freq', default='.25 .25 .25 .25',
                      help='Nucleotide freqences: -f "#A #C #G #T"')
    parser.add_option('-a', '--alpha', dest='a', type='float', default=None,
                      help='Specify alpha for heterogeneity gamma distribution')
    parser.add_option('-d', '--scale', dest='scale', type='float',
                      help='Branch length scaling factor (default = use branch lengths)')
    parser.add_option('--sbeta', dest='sbeta', type='float',
                      help='If specified, will seed gamma distribution with alpha=D and beta=SBETA as branch length scaling factor')
    parser.add_option('--smin', dest='smin', type='float', default=.01,
                      help='Specify minimum scale value (Default=0.01)')
    parser.add_option('--scalefile', dest='scalefile', default=None,
                      help='Draw branch length scaling factors from this file of numbers')
    parser.add_option('-l', '--length', dest='length', type='float', default=1000,
                      help='Sequence length (default=1000)')
    parser.add_option('--lbeta', dest='lbeta', type='float',
                      help='If specified, will seed sequence length with gamma distribution with alpha=LENGTH and beta=LBETA')
    parser.add_option('--lcor', dest='lcor', type='float',
                      help='If specified, length will correlate to # taxa in phylip file (this is slope)')
    parser.add_option('--lint', dest='lint', type='float',
                      help='Intercept of correlation between length and # taxa')
    parser.add_option('--lsd', dest='lsd', type='float', default=0.065,
                      help='Specify the st. dev. of random component of correlation b/w length and # taxa')
    parser.add_option('--lmin', dest='lmin', type='int', default=100,
                      help='Specify minimum sequence length (default=100)')
    parser.add_option('--lengthfile', dest='lengthfile', default=None,
                      help='Draw lengths from this file')
    parser.add_option('-e', '--exe', dest='exe', default=SEQ_GEN_EXE,
                      help='If you want to use a different copy of Seq-Gen, specify executable here')
    (options, args) = parser.parse_args()
    
    if len(args) < 1:
        print "No files specified"
        return
    
    if options.lmin < 50:
        print "Minimum sequence length set to 50 (it was %i)" % options.lmin
        options.lmin = 50
    
    if options.lengthfile:
        lengths = sample_from_file(options.lengthfile, len(args))
    elif options.lbeta:
        lengths = [int(x) for x in seed_gamma(len(args), options.length, options.lbeta, options.lmin)]
    elif options.lcor != None and options.lint != None:
        lengths = get_cor_lengths(args, m=options.lcor, b=options.lint, sd=options.lsd, min=options.lmin)
    else:
        lengths = [options.length] * len(args)
    
    if options.scalefile:
        scales = sample_from_file(options.scalefile, len(args))
    elif options.scale and options.sbeta:
        scales = seed_gamma(len(args), options.scale, options.sbeta, options.smin, 0.5)
    elif options.scale:
        scales = [options.scale] * len(args)
    else:
        scales = [None] * len(args)
    
    cmds = list_cmds(args, lengths, options.model, test_freq(options.freq), options.a, scales, options.exe)
    outputs = get_output_names(args)
    run_cmds(cmds, outputs)
    
if __name__ == '__main__':
    main()

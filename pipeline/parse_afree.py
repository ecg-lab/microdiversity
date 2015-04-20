#!/usr/bin/env python
"""Parse output from Afree
"""

import re, numpy, sys, os
# import pandas
from itertools import chain
from multiprocessing import Pool, TimeoutError
from optparse import OptionParser

GI_MATCH = re.compile(r'[fg]id?\|(\d+)')
NOG_MATCH = re.compile(r'NOG(\d+)\.fasta')

def store_as_int(line):
    """This function stores a line as an integer.
    
    The order of the genes does not matter (sorts them).
    """
    try:
        temp = ''.join(sorted(line.split()[:2]))
        temp = int(temp)
        if temp:
            return temp
    except:
        pass
        
    return ''.join(sorted(line.split()[:2]))

def store_as_int_tuple(x):
    """Store a tuple as an integer
    
    Rearranged in increasing order
    """
    try:
        temp = ''.join([str(y) for y in sorted(x)])
        temp = int(temp)
        if temp:
            return temp
    except:
        pass
    
    return ''.join([str(y) for y in sorted(x)])

def split_to_n_files(files, n):
    """Split a list of files into lists of size n
    """
    for i in xrange(0, len(files), n):
        yield files[i:i+n]


def filter_line(line='', NOG=False, log=False, threshold=0):
    """Filter a line on a given threshold.
    
    Does different things if NOG.
    """
    if not line:
        return
        
    if type(threshold) != float:
        try:
            threshold = float(threshold)
        except Exception as e:
            print e
            return
    
    try:
        if type(line) == str:
            temp = line.split()
        else:
            temp = list(line)
    except Exception as e:
        print e
        return
    
    if len(temp) < 3:
        return
    
    if temp[0] == temp[1]:
        return
    
    if NOG:
        if ('fid' in temp[0] and 'fid' in temp[1]) or ('gi' in temp[0] and 'gi' in temp[1]):
            # this does not have a NOG in it
            return
        if 'fid' in temp[0] or 'gi' in temp[0]:
            fidmatch = GI_MATCH.search(temp[0])
            if fidmatch:
                gene1 = fidmatch.group(1)
            nogmatch = NOG_MATCH.search(temp[1])
            if nogmatch:
                gene2 = nogmatch.group(1)
            else:
                gene2 = temp[1]
        elif 'fid' in temp[1] or 'gi' in temp[1]:
            fidmatch = GI_MATCH.search(temp[1])
            if fidmatch:
                gene1 = fidmatch.group(1)
            nogmatch = NOG_MATCH.search(temp[0])
            if nogmatch:
                gene2 = nogmatch.group(1)
            else:
                gene2 = temp[0]
        else:
            # this was two NOGs, we don't care
            return
    else:
        # this is all against all
        temp0 = GI_MATCH.search(temp[0])
        temp1 = GI_MATCH.search(temp[1])
        if temp0:
            gene1 = temp0.group(1)
        else:
            gene1 = temp[0]
        if temp1:
            gene2 = temp1.group(1)
        else:
            gene2 = temp[1]
    
    
    score = float(temp[2])
    if threshold > 0:
        if score < threshold:
            return
    if log:
        score = numpy.log10(score)
    
    return (gene1, gene2, score)

def afree_unique(files=[], output='', NOG=False, log=False, threshold=0):
    """
    """
    if not output:
        return
    if not (files and len(files)):
        return
    with open(output, 'wb') as w:
        seen = set()
        if len(files) <= 100:
            update_every_n = 1
        else:
            update_every_n = len(files) / 100
        for i, file in enumerate(files):
            if i % update_every_n == 0:
                print "Parsing file %i (%i%%)" % ((i + 1), round(100.* i / len(files)))
                sys.stdout.flush()
            
            # data = pandas.read_csv(file, sep='\t', engine='c', lineterminator='\n', header=None)
            # if data.all().all():
                # file_lines = []
                # for line in data.as_matrix():
                    # try:
                        # temp = filter_line(line=line, NOG=NOG, log=log, threshold=threshold)
                        # if temp:
                            # temp2 = store_as_int_tuple(temp)
                            # if temp2 not in seen:
                                # seen.add(temp2)
                                # file_lines.append('\t'.join([str(x) for x in temp]))
                    # except KeyboardInterrupt:
                        # return
                    # except Exception as e:
                        # print e
                        # break
                # if file_lines:
                    # w.write('\n'.join(file_lines)+'\n')
            
            
            with open(file, 'rb') as f:
                file_lines = []
                for line in f:
                    try:
                        temp = filter_line(line=line, NOG=NOG, log=log, threshold=threshold)
                        if temp:
                            temp2 = store_as_int_tuple(temp)
                            if temp2 not in seen:
                                seen.add(temp2)
                                file_lines.append('\t'.join([str(x) for x in temp]))
                    except KeyboardInterrupt:
                        return
                    except Exception as e:
                        print e
                        break
                if file_lines:
                    w.write('\n'.join(file_lines)+'\n')
    
    
                    
def read_afree(input='', NOG=False, log=False, threshold=0, string=False):
    """Read in an Afree file and filter it
    
    If string is True, return list of strings
    """
    with open(input, 'rb') as f:
        matches = []
        for line in f:
            temp = filter_line(line=line, NOG=NOG, log=log, threshold=threshold)
            if temp:
                if string:
                    temp = '\t'.join([str(x) for x in temp])+'\n'
                matches.append(temp)
    
        return matches

def _get_unique(lines):
    """Return unique lines
    
    Order of the sequences does not matter
    """
    seen = set()
    unique = []
    for line in lines:
        temp = store_as_int(line)
        if temp not in seen:
            seen.add(temp)
            unique.append(line)
    
    return unique
    
def __multi_wrapper(arguments):
    return _get_unique(read_afree(*arguments))
    
def multithreaded_afree(files=[], output='', NOG=False, log=False, threshold=0, threads=1):
    """Multi-thread Afree with filtering
    
    Can do NOGs. Can do threshold filtering. 
    This function is probably unnecessary as file I/O is the limit. Use afree_unique instead.
    """
    if not output:
        return
    if not (files and len(files)):
        return
    
    # sort files by file size
    # this tries to optimize time by grouping blocks into similar sized files
    try:
        temp = sorted(files, key=os.path.getsize)
        files = temp
    except Exception as e:
        print e
    
    # use blocks to conserve on memory!
    # write values out at end of each block
    # Block size may need to be optimized...
    blocks = split_to_n_files(files, threads)
    
    final = set()
    p = Pool(threads)
    for i, block in enumerate(blocks):
        try:
            args = [(file, NOG, log, threshold, True) for file in block]
            results = None
            afree = p.map_async(__multi_wrapper, args)
            while not results:
                try:
                    results = afree.get(10)
                except TimeoutError:
                    pass
                except KeyboardInterrupt:
                    p.terminate()
                    return
                except Exception as e:
                    print e
                    p.terminate()
                    return
            
            for result in results:
                final.update(result)
        except Exception as e:
            print e
            break
    p.close()

    if final:
        with open(output, 'wb') as w:
            # write line by line, in case something goes wrong
            # have fail safe, really want to try to save results!
            for line in final:
                for attempt in range(10):
                    try:
                        w.write(line)
                    except:
                        pass
                    else:
                        break
                else:
                    "Failed to write to file"
                    break
                    
   
    
def parse_afree(files=[], output='', NOG=False, log=False, threshold=0):
    """Read in files, single-threaded. Can filter as well.
    
    Does not remove duplicate values.
    """
    with open(output, 'wb') as w:
        for file in files:
            try:
                temp = read_afree(input=file, NOG=NOG, log=log, threshold=threshold, string=True)
                if temp:
                    w.write(''.join(temp))
            except KeyboardInterrupt:
                break
            except Exception as e:
                print e

def unique_results(lines=[], file='', output=''):
    """Return unique results.
    
    Either enter lines in lines, or enter a file to be read.
    """
    if not output:
        if file:
            output = '%s.unique' % file
        else:
            output = 'unique.mcl'
    
    if file:
        with open(file, 'rb') as r:
            with open(output, 'wb') as w:
                seen = set()
                for line in r:
                    temp = store_as_int(line)
                    if temp not in seen:
                        seen.add(temp)
                        w.write(line)
    elif lines:
        lines = set(lines)
        with open(output, 'wb') as w:
            if lines[0][-1] == '\n':
                w.write(''.join(lines))
            else:
                w.write('\n'.join(lines)+'\n')
    else:
        print "In unique_results: Nothing to do"
        return
    

def filter_mcl_file(file, output, threshold=40):
    """Utility function to filter out lines in an MCL file based on a threshold
    """
    with open(file, 'rb') as r:
        with open(output, 'wb') as w:
            i = 0
            for line in r:
                if float(line.split('\t')[2]) >= threshold:
                    w.write(line)
                    i += 1
            print "Wrote %i lines to %s" % (i, output)
    
def main():
    parser = OptionParser(usage='Usage: %prog [options] -o OUTPUT file[...]')
    parser.add_option('-o', '--out', dest='out', default='afree.mcl',
                      help='Output file for MCL input')
    parser.add_option('-t', '--threshold', dest='threshold', type='float', default=0,
                      help='Keep pairs only >= threshold (default=0)')
    parser.add_option('-n', '--nog', dest='nog', action='store_true', default=False,
                      help='Genome vs. NOG search')
    parser.add_option('-l', '--log', dest='log', action='store_true', default=False,
                      help='Log transform match scores')
    parser.add_option('-u', '--unique', dest='unique', action='store_true', default=False,
                      help='No longer used.')
    parser.add_option('--threads', type='int', default=1,
                      help='Specify threads to use')
    (options, args) = parser.parse_args()
    
    if not args:
        print "No files to parse"
        return
    
    if options.threads > 1:
        multithreaded_afree(files=args, output=options.out, NOG=options.nog, log=options.log, threshold=options.threshold, threads=options.threads)
    else:
        afree_unique(files=args, output=options.out, NOG=options.nog, log=options.log, threshold=options.threshold)
    
if __name__ == "__main__":
    main()

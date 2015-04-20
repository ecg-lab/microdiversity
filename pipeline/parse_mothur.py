#!/usr/bin/env python

from optparse import OptionParser

def parse_file(file):
    assert file
    
    with open(file, 'rb') as f:
        parsed = {}
        for line in f:
            temp = line.split()[:2]
            if temp[0] == 'unique':
                parsed['0'] = temp[1]
            else:
                parsed[temp[0]] = temp[1]
    
    return parsed

def normalize_parsed(parsed={}):
    
    if len(parsed.keys()) < 2:
        return {}
    
    floatdist = sorted(parsed.keys(), key=float)
    distmin = float(floatdist[0])
    distmax = float(floatdist[-1])
    distdenom = distmax - distmin
    
    
    floatclust = sorted(parsed.values(), key=float)
    clustmin = float(floatclust[0])
    clustmax = float(floatclust[-1])
    clutsdenom = clustmax - clustmin
    
    norm = {}
    
    for dist in parsed.keys():
        norm['%.3f'%((float(dist) - distmin)/distdenom)] = '%.3f'%((float(parsed[dist]) - clustmin) / clutsdenom)
        
    return norm

def get_all_distances(all_parsed={}):
    
    distances = []
    
    for file in all_parsed.keys():
        distances.extend(all_parsed[file].keys())
    
    return list(set(distances))
           
def combine_parsed(all_parsed={}):
    
    combined = {}

    for distance in get_all_distances(all_parsed):
        clusters = []
        for file in sorted(all_parsed.keys()):
            if distance in all_parsed[file].keys():
                clusters.append(all_parsed[file][distance])
            else:
                clusters.append('')
        combined[distance] = clusters
    
    return combined
        
def write_file(file, combined={}, header=[], delim='\t'):
    
    assert file
    assert delim
    
    # fix in case delimiter in names
    header = [ x.replace(delim,'') for x in header ]
    
    with open(file, 'wb') as w:
        
        # write header
        w.write('Distance'+delim+delim.join(header)+'\n')

        # write rows
        for dist in sorted(combined.keys(), key=float):
            w.write(dist+delim+delim.join(combined[dist])+'\n')
    
            

def parse_files(files=[], normalized=False, out='', delim='\t', minimum=3):
    
    if not (files and len(files)):
        print "No files found to parse"
        return
    
    if not out:
        print "No output file specified"
        return
        
    if not delim:
        delim = '\t'
    
    all_parsed = {}
    for file in files:
        
        if normalized:
            temp = normalize_parsed(parse_file(file))
        else:
            temp = parse_file(file)
            
        if len(temp.keys()) >= minimum:
            all_parsed[file] = temp
    
    write_file(out, combine_parsed(all_parsed), sorted(all_parsed.keys()), delim)

def main():
    
    parser = OptionParser(usage="Usage: %prog [options] -o OUTPUT file[...]")
    parser.add_option('-o', '--out', dest='out', default='mothur.txt',
                      help='Specify output filename')
    parser.add_option('-n', '--norm', dest='norm', action='store_true', default=False,
                      help='Normalize output to 0-1 scale')
    parser.add_option('-d', '--delim', dest='delim', default='\t',
                      help='Delimiter for output file (Default is \\t)')
    parser.add_option('-m', '--minimum', dest='minimum', default=3,
                      help='Minimum number of entries in file for output (default is 3)')
    (options, args) = parser.parse_args()
    
    parse_files(args, options.norm, options.out, options.delim, options.minimum)
    
if __name__ == "__main__":
    main()

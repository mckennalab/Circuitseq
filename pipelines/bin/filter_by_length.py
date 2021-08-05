#!/usr/bin/python
import gzip
import sys
import itertools
import statistics
import numpy as np

lengths = []
with gzip.open(sys.argv[1],'rt') as f:
    for name,sequence,strand,quals in itertools.zip_longest(*[f]*4):
        lengths.append(len(sequence))
 
print(len(lengths))
mean = statistics.mean(lengths)
std = statistics.stdev(lengths)
print(mean)
print(std)

# https://gist.github.com/tammoippen/4474e838e969bf177155231ebba52386
def crappyhist(a, bins=50, width=140):
    h, b = np.histogram(a, bins)

    for i in range (0, bins):
        print('{:12.5f}  | {:{width}s} {}'.format(
            b[i], 
            '#'*int(width*h[i]/np.amax(h)), 
            h[i], 
            width=width))
    print('{:12.5f}  |'.format(b[bins]))

crappyhist(lengths,20,50)

lengths = []
output = gzip.open(sys.argv[2], 'wt')
with gzip.open(sys.argv[1],'rt') as f:
    for name,sequence,strand,quals in itertools.zip_longest(*[f]*4):
        if len(sequence) < (mean + (4 * std)):
            output.write(name + sequence + strand + quals) # they already have endlines...
output.close()

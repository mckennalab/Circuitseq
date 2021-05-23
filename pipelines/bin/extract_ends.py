#!/usr/bin/python
import sys
import gzip

with gzip.open(sys.argv[1],'rt') as f:
    toggle = False
    for line in f:
        if toggle:
            print(line[0:120]+line[-120:-1])
        else:
            print(line[ 0:-1])
        toggle = not toggle

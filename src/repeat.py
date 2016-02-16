#!/usr/bin/python
# Running the same simulation script for n times.

import os
import sys

# print sys.argv

print sys.argv

nLenArg = len(sys.argv)
if nLenArg == 1:
    print 'SPECIFY the number of repeat runs and the code to repeat'
    print 'Example: repeat.py 100 geopip_test'
else:
    nRun = int(sys.argv[1])
    # print nRun
    codeToRun = sys.argv[2]
    if nLenArg == 2:
        codeParam = ''
    else:
        codeParam = ' '.join(sys.argv[3:])

code = 'python ' + codeToRun + ' ' + codeParam

# print code

for run in xrange(nRun):
    os.system(code)

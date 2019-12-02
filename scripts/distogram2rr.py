#!/usr/bin/python
## this program generates contact map (rr) file from distance map histogram (distogram) rawdistpred.current generated by DMPfold with user defined threshold and top contact selection cutoffs
import operator
import os
import math
import sys
import optparse    # for option sorting
import random
from decimal import *
getcontext().prec = 4
parser = optparse.OptionParser()
parser.add_option('-d', dest='dist',
        default = '',    # default empty
        help = 'name of raw distance map (DMPfold distogram)')
parser.add_option('-a', dest='fasta',
        default = '',    # default empty
        help = 'name of fasta file')

parser.add_option('-t', dest='threshold',
        default = '12',    # default 12
        help = 'distance threshold for distance to contact map')

parser.add_option('-r', dest='rr',
        default = 'out',    # default out
        help = 'name of output contact map')
parser.add_option('-x', dest='L',
        default = '12',    # default 12
        help = 'contact cutoff xL, where L is the sequence length')

(options,args) = parser.parse_args()
dist = options.dist
fasta = options.fasta
threshold = options.threshold
rr = options.rr
L = options.L
try:
        ffasta = open(fasta, 'r')
        flines = ffasta.readlines()
except IOError:
        sys.exit()
        print ('skipping..')

frr = open(rr, 'w')
frr.write(flines[1])
frr.write('\n')



try:
        f = open(dist, 'r')
        lines = f.readlines()
        rrList = []
        for line in lines:
                line = line.split()

                try:
                        res1 = int(line[0])
                        res2 = int(line[1])
                except ValueError:
                        print ('value error..')
                        continue

                count = 0
                for i in range(3, int(threshold)+2):

                        count += (float(line[i])) #Decimal(float(line[i]))
                if (count < 0.0001):
                        count = 0.0001
                rrList.append([str(res1), str(res2), str(threshold), str(count)])
        f.close()
        dict1 = {}
        for x in rrList:
                dict1[(x[0], x[1])] = x[3]
        sorted_rr = sorted(dict1.items(), key=operator.itemgetter(1))
        sorted_rr.reverse()
        count = 0
        for sr in sorted_rr:
                if (count >= (len(flines[1]) * float(L))):
                        break
                (i, j), p = sr[0], sr[1]

                frr.write(str(i))
                frr.write(' ')
                frr.write(str(j))
                frr.write(' ')
                frr.write('0')
                frr.write(' ')
                frr.write(threshold)
                frr.write(' ')
                frr.write(str(p))
                frr.write('\n')
                count += 1
except IOError:

        print ('skippig..')


ffasta.close()
frr.close()

#!/usr/bin/python

import numpy as np

def processbins(filechromdict, chromlist, start, segments):
    """ process chromosome dict and chromlist to produce segments 
        that are binned and averaged over bin
    """
    chrombinaverage = []
    for key, value in filechromdict.iteritems():
        startchrom = []
        segmentchrom = []
        for i, row in enumerate(chromlist):
            if key == row:
                startchrom.append(start[i])
                segmentchrom.append(segments[i])
        startchrom = np.array(startchrom).astype(float)
        segmentchrom = np.array(segmentchrom).astype(float)
        bins = np.linspace(0, value, 5)
        digitized = np.digitize(startchrom, bins)
        means = []
        keylist = []
        for i in range(1, len(bins)):
            temp = []
            keylist.append(key)
            for index, value in enumerate(digitized):
                if i == value:
                    temp.append(segmentchrom[index])
            if temp:
                means.append([bins[i]] + [np.array(temp).mean()])
            else:
                means.append([bins[i]] + [np.nan])
        print keylist
        print len(means)
        for i, item in enumerate(keylist):
            chrombinaverage.append([item] + means[i])        
    return chrombinaverage

testfile = '/Users/iManda/Desktop/final_project/GAMED_p_TCGA_B_312_313_314_NSP_GenomeWideSNP_6_E10_1361620.hg18.seg.txt'
officialchrom = '/Users/iManda/Desktop/final_project/chrom_bps.txt'

with open(officialchrom, 'r') as f:
    filechrom = [item.split('\t') for item in f.read().splitlines()]

filechromdict = dict(zip([row[0] for row in filechrom], [float(row[1]) for row in filechrom]))

with open(testfile, 'r') as f:
    cnafile = [item.split('\t') for item in f.read().splitlines()]
	
del cnafile[0]

		
chromlist = [row[1] for row in cnafile]
start = [row[2] for row in cnafile]
segments = [row[-1] for row in cnafile]

processed = processbins(filechromdict, chromlist, start, segments)

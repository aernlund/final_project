#!/usr/bin/python
#sys args should include path to files
#workingdir = "/Users/iManda/Desktop/final_project/cancer_files"


import sys
import glob
import numpy as np

def processbins(filechromdict, chromlist, start, segments):
    """ process chromosome dict and chromlist to produce segments 
        that are binned and averaged over bin, returns chromosome,
        bin, mean of bin as list of lists
    """
    chrombinaverage = []
    for key, value in filechromdict.iteritems():
        startchrom = []
        segmentchrom = []
        templist = []
        for i, row in enumerate(chromlist):
            if key == row:
                startchrom.append(start[i])
                segmentchrom.append(segments[i])
        if startchrom:
            startchromfull = np.array(startchrom).astype(float)
            segmentchromfull = np.array(segmentchrom).astype(float)
            bins = np.linspace(0, value, 5)
            digitized = np.digitize(startchromfull, bins)
            means = []
            keylist = []
            for i in range(1, len(bins)):
                temp = []
                keylist.append(key)
                for index, value in enumerate(digitized):
                    if i == value:
                        temp.append(segmentchromfull[index])
                if temp:
                    means.append([str(bins[i])] + [str(np.array(temp).mean())])
                else:
                    means.append([str(bins[i])] + [str(np.nan)])
            for i, item in enumerate(keylist):
                chrombinaverage.append([item] + means[i])        
    return chrombinaverage

def process_inputfile(inputfile, filechromdict):
    """ take input cnafile, process it, and return cancer, processed
    """
    start = [row[1] for row in inputfile]
    chromlist = [row[0] for row in inputfile]
    segments = [row[2] for row in inputfile]
    cancername = inputfile[3][3]
    processed = processbins(filechromdict, chromlist, start, segments)
    cancer = [cancername] * len(processed)
    process_matrix = []
    for i in range(len(processed)):
        process_matrix.append([cancer[i]] + processed[i])
    return process_matrix
    
    
    
#workingdir = str(sys.argv[1])
workingdir = "/Users/iManda/Desktop/final_project/cancer_files"
masterlist = glob.glob('{0}/*.txt'.format(workingdir))
officialchrom = '/Users/iManda/Desktop/final_project/chrom_bps.txt'
    
with open(officialchrom, 'r') as f:
    filechrom = [item.split('\t') for item in f.read().splitlines()]

filechromdict = dict(zip([row[0] for row in filechrom], [float(row[1]) for row in filechrom]))
	
matrix = []
for thing in masterlist:
    with open(thing, 'r') as f:
        inputfile = [item.split('\t') for item in f.read().splitlines()]
    processedfile = process_inputfile(inputfile, filechromdict)
    for row in processedfile:
        matrix.append(row)
        
with open('{0}/all_cancer.txt'.format(workingdir), 'w') as f:
    for row in matrix:
        f.write('\t'.join(row) + '\n')

    


#%prun barcode_file(manifestfile)
#%prun processfile(copynumfile)




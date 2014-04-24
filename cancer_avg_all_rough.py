#!/usr/bin/python
"""This script processes all cancers to give average bin copy number"""

import sys
import glob
import numpy as np

def barcode_file(manifestfile):
    """ (manifest_file) -> dictionary of filenames: barcodes, names of files (excluding .txt)
    """
    cnafilename = [item[-1].split('.')[0] for item in manifestfile]
    tcgabarcode = ['-'.join(item[-2].split('-')[0:3]) for item in manifestfile]
    barcodename = dict(zip(cnafilename, tcgabarcode))
    return barcodename
    
        

def processbins(filechromdict, chromlist, start, segments):
    """ process chromosome dict and chromlist to produce segments 
        that are binned and averaged over bin
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
    
def process_file(cnafile, barcode, filechromdict):
    """ take input cnafile, process it, and return patient, cancer, processed
    """
    samplename = cnafile[0][0]
    tcgacode = barcode[samplename]
    chromlist = [row[1] for row in cnafile]
    start = [row[2] for row in cnafile]
    segments = [row[-1] for row in cnafile]
    processed = processbins(filechromdict, chromlist, start, segments)
    cancer = [cancername] * len(processed)
    patientcode = [tcgacode] * len(processed)
    process_matrix = []
    for i in range(len(processed)):
        process_matrix.append([patientcode[i]] + [cancer[i]] + processed[i])
    return process_matrix            


masterlist = glob.glob('./*.txt')
#cancername = str(sys.argv[1])
#workingdir = str(sys.argv[2])

cancername = 'ACC'
#workingdir = '/Users/iManda/Desktop/final_project/ACC_test'
officialchrom = '/Users/iManda/Desktop/final_project/chrom_bps.txt'

with open(officialchrom, 'r') as f:
    filechrom = [item.split('\t') for item in f.read().splitlines()]

filechromdict = dict(zip([row[0] for row in filechrom], [float(row[1]) for row in filechrom]))

#with open('{0}/{1}/file_manifest.txt'.format(workingdir, cancername), 'r') as f:
    #manifestfile = [item.split('\t') for item in f.read().splitlines()]

with open('/Users/iManda/Desktop/final_project/ACC_test/file_manifest.txt') as f:
    manifestfile = [item.split('\t') for item in f.read().splitlines()]

del manifestfile[0]
del manifestfile[0]
del manifestfile[0]

barcode = barcode_file(manifestfile)


#fix so that hg18 vs hg19 are averaged

matrix = []
for file in masterlist:
    temp = file.split('.')[2]
    if temp[0:5] != 'nocnv':
        with open(file, 'r') as f:
            cnafile = [item.split('\t') for item in f.read().splitlines()]
        del cnafile[0]
        #print file
        if file[-12:] == 'hg18.seg.txt':
            processed_file = process_file(cnafile, barcode, filechromdict)
            reference = ['hg18'] * len(processed_file)
            for i in range(len(processed_file)):
                matrix.append(processed_file[i] + [reference[i]])
        if file[-12:] == 'hg19.seg.txt':
            processed_file = process_file(cnafile, barcode, filechromdict)
            reference = ['hg19'] * len(processed_file)
            for i in range(len(processed_file)):
                matrix.append(processed_file[i] + [reference[i]])





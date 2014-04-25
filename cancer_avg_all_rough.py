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
        that are binned and averaged over bin, returns chromosome,
        bin,  and mean of bin as list of lists
    """
    chrombinaverage = []
    for key, value in filechromdict.iteritems():
        startchrom = []
        segmentchrom = []
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
    """ take input cnafile, returns patient, cancer, processed
        as list of lists
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


cancername = str(sys.argv[1])
workingdir = str(sys.argv[2])
masterlist = glob.glob('{0}/{1}/CNV_SNP_Array/BI__Genome_Wide_SNP_6/Level_3/*.txt'.format(workingdir, cancername))
officialchrom = '/Users/iManda/Desktop/final_project/chrom_bps.txt'


with open(officialchrom, 'r') as f:
    filechrom = [item.split('\t') for item in f.read().splitlines()]

filechromdict = dict(zip([row[0] for row in filechrom], [float(row[1]) for row in filechrom]))

with open('{0}/{1}/file_manifest.txt'.format(workingdir, cancername), 'r') as f:
    manifestfile = [item.split('\t') for item in f.read().splitlines()]

del manifestfile[0]
del manifestfile[0]
del manifestfile[0]

barcode = barcode_file(manifestfile)


matrix = []
for file in masterlist:
    temp = file.split('.')[1]
    if temp[0:5] != 'nocnv' and temp[-4:] == 'hg19':
        with open(file, 'r') as f:
            cnafile = [item.split('\t') for item in f.read().splitlines()]
        del cnafile[0]
        processed= process_file(cnafile, barcode, filechromdict)
        for i in range(len(processed)):
            matrix.append(processed[i])

with open('{0}/{1}_cnabinned_matrix.txt'.format(workingdir, cancername), 'w') as f:
    for row in matrix:
        f.write('\t'.join(row) + '\n')



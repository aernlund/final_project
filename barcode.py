#!/usr/bin/python
#run with pylint for one file
#run pylint python ACC /Users/iManda/Desktop/final_project
#import sys
#import glob
"""This script process the copy number array files to put into the database"""

workingdir = "/Users/iManda/Desktop/final_project"
cancername = "ACC"
testfile = "CNV_SNP_Array/BI__Genome_Wide_SNP_6/Level_3/GAMED_p_TCGA_B_312_313_314_NSP_GenomeWideSNP_6_E10_1361620.hg18.seg.txt"


def barcode_file(manifestfile):
    """ (manifest_file) -> list of lists of barcodes names of files (excluding .txt)
    """
    cnafilename = [item[-1].split('.')[0] for item in manifestfile]
    tcgabarcode = ['-'.join(item[-2].split('-')[0:3]) for item in manifestfile]
    barcodename = zip(cnafilename, tcgabarcode)
    return barcodename
    
def processfile(copynumfile):
    """(extract columns from hg18 files) -> list of lists with all file info
    """
    format = []
    samplefilename = [item[0] for item in copynumfile]
    shortbarcode = []
    for item in samplefilename:
        for row in barcode:
            if item == row[0]:
                shortbarcode.append(row[1])		 
    chr = ['chr'+item[1] for item in copynumfile]
    startend = [item[2:4] for item in copynumfile]
    mean = [item[-1] for item in copynumfile]
    cancer = [cancername] * len(mean)
    for i in range(len(mean)):
        format.append([shortbarcode[i]] + [chr[i]] + startend[i] + [mean[i]] + [cancer[i]])
    return format

with open('{0}/{1}/file_manifest.txt'.format(workingdir, cancername), 'r') as f:
    manifestfile = [item.split('\t') for item in f.read().splitlines()]
    
with open('{0}/{1}/{2}'.format(workingdir, cancername, testfile), 'r') as f:
    copynumfile = [item.split('\t') for item in f.read().splitlines()]
    

del manifestfile[0]
del manifestfile[0]
del manifestfile[0]

del copynumfile[0]
	
barcode = barcode_file(manifestfile)
processed_file = processfile(copynumfile)
#reference = ['hg18'] * len(mean)


print barcode[0]
print processed_file[0]



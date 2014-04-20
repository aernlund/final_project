#!/usr/bin/python
#run with pylint for one file
#run pylint python ACC /Users/iManda/Desktop/final_project
#import sys
#import glob
"""This script process the copy number array files to put into the database"""

workingdir = "/Users/iManda/Desktop/final_project"
cancername = "ACC"


def barcode_file(manifestfile):
    """ (manifest_file) -> list of lists of barcodes names of files (excluding .txt)
    """
    cnaname = [item[-1].split('.')[0] for item in manifestfile]
    tcganame = ['-'.join(item[-2].split('-')[0:3]) for item in manifestfile]
    barcodename = zip(cnaname, tcganame)
    return barcodename

with open('{0}/{1}/file_manifest.txt'.format(workingdir, cancername), 'r') as f:
    manifestfile = [item.split('\t') for item in f.read().splitlines()]


del manifestfile[0]
del manifestfile[0]
del manifestfile[0]
	
barcode = barcode_file(manifestfile)

print barcode[0]

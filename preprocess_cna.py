#!/usr/bin/python
#run with pylint for one file
#run pylint python ACC /Users/iManda/desktop/Database_methods
import sys
import glob
#import re

masterlist = glob.glob('./*.txt')
cancername = str(sys.argv[1])
workingdir = str(sys.argv[2])
#print cancername

with open('{0}/{1}/file_manifest.txt'.format(workingdir, cancername), 'r') as f:
	file = [item.split('\t') for item in f.read().splitlines()]


del file[0]
del file[0]
del file[0]
	
cnaname = [item[-1].split('.')[0] for item in file]
TCGAname = ['-'.join(item[-2].split('-')[0:3]) for item in file]
barcode = zip(cnaname, TCGAname)
barcodeuniq = []
for item in barcode:
	if item not in barcodeuniq:
		barcodeuniq.append(item)


cna18list = []
cna19list = []

for file in masterlist:
	if file[-12:] == 'hg18.seg.txt':
		with open(file, 'r') as f:
			file18 = [item.split('\t') for item in f.read().splitlines()]
		del file18[0]
		sample18name = [item[0] for item in file18]
		TCGA18 = []
		for item in sample18name:
			for row in barcode:
				if item == row[0]:
					TCGA18.append(row[1])		 
		chr = ['chr'+item[1] for item in file18]
		startend = [item[2:4] for item in file18]
		mean = [item[-1] for item in file18]
		cancer = [cancername] * len(mean)
		reference = ['hg18'] * len(mean)
		for i in range(len(mean)):
			cna18list.append([TCGA18[i]] + [chr[i]] + startend[i] + [mean[i]] + [cancer[i]] + [reference[i]])
		#print 'finished with file hg18'
 	if file[-12:] == 'hg19.seg.txt':
		with open(file, 'r') as f:
			file19 = [item.split('\t') for item in f.read().splitlines()]
		del file19[0]
		sample19name = [item[0] for item in file19]
		TCGA19 = []
		for item in sample19name:
			for row in barcode:
				if item == row[0]:
					TCGA19.append(row[1])		 
		chr = ['chr'+item[1] for item in file19]
		startend = [item[2:4] for item in file19]
		mean = [item[-1] for item in file19]
		cancer = [cancername] * len(mean)
		reference = ['hg19'] * len(mean)
		for i in range(len(mean)):
			cna19list.append([TCGA19[i]] + [chr[i]] + startend[i] + [mean[i]] + [cancer[i]] + [reference[i]])
		#print 'finished with file hg19'

			

			
header = ['PatientID', 'chrom', 'loc.start', 'loc.end', 'seg.mean', 'cancer', 'reference']
with open('{0}/{1}_cna_matrix.txt'.format(workingdir, cancername), 'w') as f:
	f.write('\t'.join(header) + '\n')
	for row in cna18list:
		f.write('\t'.join(row) + '\n')
	for row in cna19list:
		f.write('\t'.join(row) + '\n')


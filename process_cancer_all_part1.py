#!/usr/bin/python
#workingdir = "/Users/iManda/Desktop/final_project"
#cancername = "ACC"


import sys
import glob

def file_info(copynumfile):
    """(extract columns from hg18 files) -> list of lists with all file info
    """
    format = []
    chr = [item[1] for item in copynumfile]
    start = [item[2] for item in copynumfile]
    mean = [item[-1] for item in copynumfile]
    cancer = [cancername] * len(mean)
    for i in range(len(mean)):
        format.append([chr[i]] + [start[i]] + [mean[i]] + [cancer[i]])
    return format
    
cancername = str(sys.argv[1])
workingdir = str(sys.argv[2])
masterlist = glob.glob('{0}/{1}/CNV_SNP_Array/BI__Genome_Wide_SNP_6/Level_3/*.txt'.format(workingdir, cancername))


	
matrix = []
for file in masterlist:
    temp = file.split('.')[1]
    if temp[0:5] != 'nocnv' and temp[-4:] == 'hg19':
        with open(file, 'r') as f:
            cnafile = [item.split('\t') for item in f.read().splitlines()]
        del cnafile[0]
        processed_file = file_info(cnafile)
        for i in range(len(processed_file)):
            matrix.append(processed_file[i])

with open('{0}/{1}_cna_matrix.txt'.format(workingdir, cancername), 'w') as f:
    for row in matrix:
        f.write('\t'.join(row) + '\n')

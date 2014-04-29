from cancer_avg_all_functions import processbins
import unittest

#should pass original file
incomingfile = '/Users/iManda/Desktop/final_project/patient_input_files/ACC_patient_hg19_cna.txt'
with open(incomingfile ,'r') as f:
	inputfile = [item.split('\t') for item in f.read().splitlines()]

chromlist = [item[0] for item in inputfile]
start = [item[1] for item in inputfile]
segments = [item[2] for item in inputfile]

officialchrom = '/Users/iManda/Desktop/final_project/chrom_bps.txt'

with open(officialchrom, 'r') as f:
    filechrom = [item.split('\t') for item in f.read().splitlines()]

filechromdict = dict(zip([row[0] for row in filechrom], [float(row[1]) for row in filechrom]))

#file without y chromosome - this should pass
incomingfile2 = '/Users/iManda/Desktop/final_project/patient_input_files/BLCA_patient_hg19_cna_noy.txt'
with open(incomingfile ,'r') as f:
    inputfile2 = [item.split('\t') for item in f.read().splitlines()]

chromlist2 = [item[0] for item in inputfile2]
start2 = [item[1] for item in inputfile2]
segments2 = [item[2] for item in inputfile2]
filechromdict2 = dict(zip([row[0] for row in filechrom], [float(row[1]) for row in filechrom]))

class TestProcessBins(unittest.TestCase):
    """ Test class for function python_hw2.is_palindrome """
    def test_processbins(self):
        self.assertTrue(processbins(filechromdict, chromlist, start, segments))
    
    def test_processbins(self):
        self.assertTrue(processbins(filechromdict2, chromlist2, start2, segments2))

if __name__ == '__main__':
    unittest.main()

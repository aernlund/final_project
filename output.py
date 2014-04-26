import MySQLdb as MS
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

#input file
#modify to accept command line file input
incomingfile = 'ACC_patient_hg19_cna.txt'
with open(incomingfile ,'r') as f:
	inputfile = [item.split('\t') for item in f.read().splitlines()]
	
def process_inputfile(inputfile, filechromdict):
    """(extract columns from hg18 files) -> list of lists with 
        chromosome, bin, means
    """
    chromlist = [item[0] for item in inputfile]
    start = [item[1] for item in inputfile]
    segments = [item[3] for item in inputfile]
    process_matrix = processbins(filechromdict, chromlist, start, segments)
    return process_matrix


def processbins(filechromdict, chromlist, start, segments):
    """ process chromosome dict and chromlist to produce segments 
        that are binned and averaged over bin, returns chromosome,
        bin, mean of bin as list of lists
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
    
def patient_rank(df, uniq, inputmean):
    """cancer dataframe, list of unique cancers, mean column from inputdata
       return an array containing the top cancer and correlation value
    """
    corr = 0
    for item in uniq:
        inputdf['chromosome'] = inputdf['chromosome'].values.astype('str')
        temp = df.ix[df['PatientID']== item,:]
        newtemp = df.ix[df['chromosome'] != 'Y', :]	
        sorted_cancer = newtemp.sort(['chromosome','bin'], ascending=[True,True])
        cancer_mean = sorted_cancer['mean'].values
        new_matrix = np.vstack([inputmean, cancer_mean]).transpose()
        means = pd.DataFrame(new_matrix, columns = ['input_mean', 'cancer_mean'])
        means2 = means.dropna(how='any')
        newcorr = np.square(np.corrcoef(means2['input_mean'], means2['cancer_mean'])[0][-1])
        if newcorr > corr:
            corr = newcorr
            patient = item
            top_patient = [patient, corr]
            top_patient_means = means2
    return top_patient, top_patient_means
    
def cancer_rank(df, uniq, inputmean, ytest):
    """cancer dataframe, list of unique cancers, mean column from inputdata
       return an array containing the top cancer and correlation value
    """
    corr = 0
    for item in uniq:
        inputdf['chromosome'] = inputdf['chromosome'].values.astype('str')
        temp = df.ix[df['cancer']== item,:]
        if ytest:
            sorted_cancer = temp.sort(['chromosome','bin'], ascending=[True,True])
            cancer_mean = sorted_cancer['mean'].values
            new_matrix = np.vstack([inputmean, cancer_mean]).transpose()
            means = pd.DataFrame(new_matrix, columns = ['input_mean', 'cancer_mean'])
            means2 = means.dropna(how='any')
            newcorr = np.square(np.corrcoef(means2['input_mean'], means2['cancer_mean'])[0][-1])
            if newcorr > corr:
                corr = newcorr
                cancer = item
                top_cancer = [cancer, corr]
                top_means = means2
        else:
            newtemp = df.ix[df['chromosome'] != 'Y', :]	
            cancer_mean = sorted_cancer['mean'].values
            new_matrix = np.vstack([inputmean, cancer_mean]).transpose()
            means = pd.DataFrame(new_matrix, columns = ['input_mean', 'cancer_mean'])
            means2 = means.dropna(how='any')
            newcorr = np.square(np.corrcoef(means2['input_mean'], means2['cancer_mean'])[0][-1])
            if newcorr > corr:
                corr = newcorr
                cancer = item
                top_cancer = [cancer, corr]
                top_means = means2
    return top_cancer, top_means


         


officialchrom = 'chrom_bps.txt'
    
with open(officialchrom, 'r') as f:
    filechrom = [item.split('\t') for item in f.read().splitlines()]

#input file processing
filechromdict = dict(zip([row[0] for row in filechrom], [float(row[1]) for row in filechrom]))
processedfile = process_inputfile(inputfile, filechromdict)
inputdf = pd.DataFrame(processedfile, columns=['chromosome', 'bin', 'mean'])
y = inputdf['chromosome'].values.astype('str')
ytest = 'Y' in temp
sortedinput = inputdf.sort(['chromosome','bin'], ascending=[True,True])
inputmean = sortedinput['mean'].values.astype('float')


#cancer matrix
conn = MS.connect(host="localhost", user="awe220", passwd="coco2828", db="cancer_copy_number")
df = pd.read_sql("SELECT * FROM all_cancer_copy_number;", conn)
uniq = df['cancer'].unique()
top_cancer, top_means = cancer_rank(df, uniq, inputmean)

#Graph of top cancer
m, b = np.polyfit(top_means['input_mean'], top_means['cancer_mean'], 1)
plt.clf()
plt.plot(top_means['input_mean'], top_means['cancer_mean'], '.')
plt.plot(top_means['input_mean'], m*top_means['input_mean'] + b, '-')
plt.xlabel('input_mean')
plt.ylabel('{0}_mean'.format(top_cancer[0])
plt.title('Top Cancer Correlation')
#plt.legend(loc = 'lower right', scatterpoints = 0)
#plt.text(xytext='axes fraction','R^2= %.2f'%corr, verticalalignment='bottom', horizontalalignment='right')
plt.savefig('myfig')

#specific cancer matrix
spdf = pd.read_sql("SELECT * FROM {0}_copy_number;".format(top_cancer[0]),conn)
spdf.save("ACC_specific_matrix.txt")
df = pd.load("ACC_specific_matrix.pkl")
uniq_patient = spdf['PatientID'].unique()
new_inputmean_mat = sortedinput.ix[sortedinput['chromosome'] != 'Y',:] 
new_inputmean = new_inputmean_mat['mean'].values.astype('float')
top_patient, top_patmeans = patient_rank(spdf, uniq_patient, new_inputmean) 



































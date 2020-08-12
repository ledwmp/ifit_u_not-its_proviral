import sys
import numpy as np
from scipy.stats import ranksums,ttest_ind
from sklearn.neighbors import KernelDensity
#from sklearn.grid_search import GridSearchCV
#from sklearn.cross_validation import LeaveOneOut
from collections import Counter
import matplotlib.pyplot as plt

FPKM_dict = {}

with open("/home/mitch/Desktop/IFIT2_CLIP_temp/alt_pvalue/protein_coding/meta_intersect/chi2_geneontology/FPKM_group_novogene_filterwFPKM.txt") as r:
	for line in r:
		FPKM_dict[line.split("\t")[0]] = float(line.split("\t")[1].strip())

CLIP_dict = {}
for item in sys.argv[1:-1]:
	with open(item) as r:
		for line in r:
			if line.strip() in FPKM_dict.keys():
				CLIP_dict[line.strip()] = FPKM_dict[line.strip()]

log_FPKM = [np.log10(i+1) for i in FPKM_dict.values()]
log_CLIP = [np.log10(i+1) for i in CLIP_dict.values()]

def calc_gc(nuc_string):
	len_list = float(len(nuc_string))
	count_nuc = Counter(nuc_string)
	return (count_nuc["G"]+count_nuc["C"])/len_list

def kde_sklearn(x,x_grid):
	#bandwidth = find_bandwidth(x)
	#print bandwidth
	kde_sk1 = KernelDensity(kernel="gaussian",bandwidth=0.1).fit(x[:,np.newaxis])
	log_pdf = kde_sk1.score_samples(x_grid)
	return np.exp(log_pdf)

def return_xgrid(tmp_list):
	new_list = sorted(tmp_list)
	high = new_list[0]
	low = new_list[-1]
	points = len(tmp_list)
	#return np.linspace(low,high,1000)[:, np.newaxis]
	return np.linspace(0.0,4.5,1000)[:,np.newaxis]

#def find_bandwidth(x):
	bandwidths = 10**np.linspace(-1,1,100)
	grid = GridSearchCV(KernelDensity(kernel="gaussian"),{"bandwidth":bandwidths},cv=LeaveOneOut(len(x)))
	grid.fit(x[:,None]);
	return grid.best_params
ranksums,pvalue = ranksums(log_FPKM,log_CLIP)
print str(pvalue)
ttestind,pvalue = ttest_ind(log_FPKM,log_CLIP,equal_var=False)
print str(pvalue)

pdf_1 = kde_sklearn(np.array(log_CLIP),return_xgrid(log_CLIP))
pdf_2 = kde_sklearn(np.array(log_FPKM),return_xgrid(log_FPKM))

fit1,ax1 = plt.subplots()
ax1.plot(return_xgrid(log_CLIP),pdf_1,color='gray',alpha=0.5,lw=3,label="CLIP")
ax1.plot(return_xgrid(log_FPKM),pdf_2,color='blue',alpha=0.5,lw=3,label="Total")
plt.title(sys.argv[-1]+"\nWelch's unpaired ttest p="+str(pvalue)+"\nCLIP_n="+str(len(log_CLIP))+" total_n="+str(len(log_FPKM)),fontsize=10)
plt.legend(loc="upper left")
plt.xlabel("log10(FPKM+1)")
plt.ylabel("Probability")
plt.tick_params(direction="out")
plt.xlim(0.0,)
plt.ylim(0.0,)
#plt.savefig(sys.argv[-1]+".svg")
plt.show()

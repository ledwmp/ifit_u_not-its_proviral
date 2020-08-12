#!/usr/bin python
_author_ = "Mitch Ledwith"
"""
arguments fasta_peak_to_gc_kde.py CLIP_peaks.fasta sub_region.fasta out_file
"""
from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
import sys 
from scipy.stats import ranksums,ttest_ind
from sklearn.neighbors import KernelDensity
#from sklearn.grid_search import GridSearchCV
#from sklearn.cross_validation import LeaveOneOut

def calc_gc(nuc_string):
	len_list = float(len(nuc_string))
	count_nuc = Counter(nuc_string)
	return (count_nuc["G"]+count_nuc["C"])/len_list

def kde_sklearn(x,x_grid):
	#bandwidth = find_bandwidth(x)
	#print bandwidth
	kde_sk1 = KernelDensity(kernel="gaussian",bandwidth=0.03).fit(x[:,np.newaxis])
	log_pdf = kde_sk1.score_samples(x_grid)
	return np.exp(log_pdf)

def return_xgrid(tmp_list):
	new_list = sorted(tmp_list)
	high = new_list[0]
	low = new_list[-1]
	points = len(tmp_list)
	return np.linspace(0.0,1.0,1000)[:,np.newaxis]
	#return np.linspace(low,high,1000)[:, np.newaxis]

#def find_bandwidth(x):
	#bandwidths = 10**np.linspace(-1,1,100)
	#grid = GridSearchCV(KernelDensity(kernel="gaussian"),{"bandwidth":bandwidths},cv=LeaveOneOut(len(x)))
	#grid.fit(x[:,None]);
	#return grid.best_params
	
	
	
plot_list = []

for item in sys.argv[1:-1]:
	temp_list = [] 
	for record in SeqIO.parse(item,"fasta"):
		nucleotide_str = str(record.seq).upper()
		temp_list.append(calc_gc(nucleotide_str))
	plot_list.append(temp_list)

ranksums,pvalue = ranksums(plot_list[0],plot_list[1])
print str(pvalue)
#ttestind,pvalue = ttest_ind(plot_list[0],plot_list[1],equal_var=False)
print str(pvalue)


pdf_1 = kde_sklearn(np.array(plot_list[0]),return_xgrid(plot_list[0]))
pdf_2 = kde_sklearn(np.array(plot_list[1]),return_xgrid(plot_list[1]))



fit1,ax1 = plt.subplots()
ax1.plot(return_xgrid(plot_list[0]),pdf_1,color='blue',alpha=0.65,lw=3,label="CLIP n="+str(len(plot_list[0])))
ax1.plot(return_xgrid(plot_list[1]),pdf_2,color='gray',alpha=0.65,lw=3,label="total n="+str(len(plot_list[1])))
plt.text(0.6,3.5,"Mann-Whitney U-test\np="+str(pvalue))
plt.legend(loc="upper right")
plt.xlabel("GC-content (%)",fontsize=12)
plt.tick_params(direction="out")
plt.ylabel("Probability",fontsize=12)
plt.xlim(0,1)
plt.ylim(0,5.0)
plt.savefig(sys.argv[-1]+".svg")
plt.show()



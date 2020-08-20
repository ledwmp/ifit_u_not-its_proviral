#!/usr/bin python
_author_ = "Mitch Ledwith"
"""
arguments cumulative_distribution_CLIP_combinepeaks_delta_final DE_file.txt FPKM_file.txt
Plots cumulative distribution curves
"""
import sys
import numpy as np 
import matplotlib.pyplot as plt
import scipy.stats as stats
import random
from collections import defaultdict

random.seed(30)

#files: pause, rpkm

gene_id = np.loadtxt(sys.argv[1],dtype=str,delimiter="\t",skiprows=1, usecols=(0,))
pvalue = np.loadtxt(sys.argv[1],dtype=str,delimiter="\t",skiprows=1, usecols=(4,))
logfc = np.loadtxt(sys.argv[1],dtype=float,delimiter="\t",skiprows=1, usecols=(1,))
gene_idd = np.loadtxt(sys.argv[2],dtype=str,delimiter="\t", usecols=(0,))
logrpkm = np.loadtxt(sys.argv[2],dtype=float,delimiter="\t", usecols=(1,))

dict_rpkm = {k.strip("+"):v for k,v in zip(gene_idd,np.log2(logrpkm))}
filt_rpkm = {j:v for k,v in dict_rpkm.iteritems() for j in k.split(",")}
dict_p = {k:v for k,v in zip(gene_id,pvalue)}
dict_fc = {k:v for k,v in zip(gene_id,logfc)}

temp_set = set(filt_rpkm.keys())

for item in [i.split(":")[0].strip('"') for i in gene_id]:
	if temp_set.isdisjoint([item]) == True:
		filt_rpkm[item] = 0.0 


dict_c = {k:(dict_p[k],dict_fc[k]) for k in dict_p.keys() if filt_rpkm[k.split(":")[0].strip('"')] > 3.32}


total_p_dict = defaultdict(list)
total_fc_dict = defaultdict(list)
total_ribonorm_dict = defaultdict(list)
individual_ribonorm_dict = defaultdict(list)

for k,v in dict_c.iteritems():
	new_k = k.split(":")[0].strip('"')
	total_p_dict[new_k].append(v[0])
	total_fc_dict[new_k].append(v[1])
	total_ribonorm_dict[new_k].append(v[1])
	individual_ribonorm_dict[k] = v[1]


max_fc_dict = {k:np.log2(2**np.max(v)/2**np.min(v)) for k,v in total_ribonorm_dict.iteritems()}


for i,j in max_fc_dict.iteritems():
	print i,j

plt.figure(figsize=(5,5))

CLIP_list = np.loadtxt("CLIP_genes.txt",dtype=str,delimiter="\t", usecols=(0,))

CLIP_zip = {}
for item in CLIP_list:
	for k,v in max_fc_dict.iteritems():
		if item in k:
			CLIP_zip[k] = v
	
dict_b = {k:v for k,v in max_fc_dict.iteritems() if k not in CLIP_zip.keys()}

sorted_total = sorted(dict_b.values())
sorted_CLIP = sorted(CLIP_zip.values())

total_y = []

for z in range(0,len(sorted_total)):
	total_y.append(float(1.0/len(sorted_total))*z)
print len(total_y)

CLIP_y = []

for z in range(0,len(sorted_CLIP)):
	CLIP_y.append(float(1.0/len(sorted_CLIP))*z)
print len(CLIP_y)


statistic,pvalue = stats.ranksums(sorted_total,sorted_CLIP)
print pvalue

params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)

plt.scatter(sorted_total,total_y,s=1,color='k',alpha=0.5,label="not CLIP n="+str(len(sorted_total)))
plt.scatter(sorted_CLIP,CLIP_y,s=1,color='r',alpha=0.5,label="IFIT2 CLIP n="+str(len(sorted_CLIP)))

lgnd = plt.legend(loc="upper left",scatterpoints=1,fontsize=10)
lgnd.legendHandles[0]._sizes = [12] 
lgnd.legendHandles[1]._sizes = [12]

plt.text(0.05,0.75,"M-W U p-value=\n"+str(pvalue))

plt.xlabel('$log_{2}FC  \Delta^{max}_{min}${'+"$peak^{IFIT2^{-/-}}/peak^{WT}$"+"} / ($RPF^{IFIT2^{-/-}}/RPF^{WT})$",fontsize=10)
plt.ylabel('Cumulative Fraction',fontsize=10)
plt.xlim(0,3)
plt.ylim(0,1)


plt.savefig("ixndeltapause.svg")
plt.show()

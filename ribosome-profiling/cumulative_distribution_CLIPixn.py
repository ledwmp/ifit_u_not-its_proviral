#!/usr/bin python
_author_ = "Mitch Ledwith"
"""
arguments cumulative_distribution_CLIPixn.py DE_file.txt
Plots cumulative distribution curves
"""

import sys
import numpy as np 
import matplotlib.pyplot as plt
import scipy.stats as stats
import random

random.seed(30)


gene_id = np.loadtxt(sys.argv[1],dtype=str,delimiter="\t",skiprows=1, usecols=(0,))
logfc = np.loadtxt(sys.argv[1],dtype=float,delimiter="\t",skiprows=1, usecols=(1,))
logcpm = np.loadtxt(sys.argv[1],dtype=float,delimiter="\t",skiprows=1, usecols=(2,))

dict_a = {k:v for k,v,j in zip(gene_id,logfc,logcpm) if j > 5}

plt.figure(figsize=(5,5))

CLIP_list = np.loadtxt("IFIT2ixn_CLIPgenes.txt",dtype=str,delimiter="\t", usecols=(0,))

CLIP_zip = {}
for item in CLIP_list:
	for k,v in dict_a.iteritems():
		if item in k:
			CLIP_zip[k] = v

dict_b = {k:v for k,v in dict_a.iteritems() if (k not in CLIP_zip.keys() and "WSN" not in k)}
WSN_dict = {k:v for k,v in dict_a.iteritems() if "WSN" in k}

sorted_total = sorted(dict_b.values())
sorted_CLIP = sorted(CLIP_zip.values())
sorted_WSN = sorted(WSN_dict.values())
sorted_random = sorted(random.sample(dict_b.values(),len(CLIP_zip.values())))

total_y = []

for z in range(0,len(sorted_total)):
	total_y.append(float(1.0/len(sorted_total))*z)

CLIP_y = []

for z in range(0,len(sorted_CLIP)):
	CLIP_y.append(float(1.0/len(sorted_CLIP))*z)

random_y = []

for z in range(0,len(sorted_random)):
	random_y.append(float(1.0/len(sorted_random))*z)

WSN_y = []

for z in range(0,len(sorted_WSN)):
	WSN_y.append(float(1.0/len(sorted_WSN))*z)

statistic,pvalue_CLIP = stats.ranksums(sorted_total,sorted_CLIP)
print pvalue_CLIP


params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)

plt.scatter(sorted_total,total_y,s=1,color='k',alpha=0.5,label="not CLIP n="+str(len(sorted_total)))
plt.scatter(sorted_CLIP,CLIP_y,s=1,color='r',alpha=0.5,label="IFIT2 CLIP n="+str(len(sorted_CLIP)))
#plt.scatter(sorted_random,random_y,s=1,color='g',alpha=0.5)
lgnd = plt.legend(loc="upper left",scatterpoints=1,fontsize=10)
lgnd.legendHandles[0]._sizes = [12] 
lgnd.legendHandles[1]._sizes = [12] 
plt.text(-1.45,0.75,"M-W U p-value=\n"+str(pvalue_CLIP))

plt.xlabel('$log_{2}FC$ '+"$TE^{IFIT2^{-/-}}$"+"/ "+"$TE^{WT}$  infection",fontsize=10)
plt.ylabel('Cumulative Fraction',fontsize=10)
plt.axvline(0.0,c="k",linewidth=1,linestyle="--",alpha=0.5)

plt.xlim(-1.5,1.5)
#plt.xlim(-1.0,1.0)
plt.ylim(0,1)

plt.savefig("ixnTE.svg")
plt.show()

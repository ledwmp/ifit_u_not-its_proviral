#!/usr/bin python
_author_ = "Mitch Ledwith"

"""
arguments count_CLIP_metabiotype_new.py intersect1.bed intersect2.bed out_file
"""
import sys
from collections import Counter
import matplotlib.pyplot as plt


final_list = []
len_list = []
for item in sys.argv[1:-1]:
	tmp_list = []
	with open(item) as r:
		for line in r:
			new_line = line.split("\t")
			if "protein_coding" in new_line[9]:
				tmp_list.append("protein_coding")
			elif "lincRNA" in new_line[9]:
				tmp_list.append("lncRNA")
			elif "pseudogene" in new_line[9]:
				tmp_list.append("pseudogene")
			else:
				tmp_list.append("other")
		blah_blah = Counter(tmp_list)
	len_list.append(len(tmp_list))
	final_list.append([blah_blah["protein_coding"]/float(len(tmp_list)),blah_blah["lncRNA"]/float(len(tmp_list)),blah_blah["pseudogene"]/float(len(tmp_list)),blah_blah["other"]/float(len(tmp_list))])

print final_list

barwidth = 0.4
r = [0.25,0.75]
plt.bar(r,[i[0] for i in final_list],color='gray',label="protein_coding",alpha=0.5,width=barwidth,edgecolor="white")
plt.bar(r,[i[1] for i in final_list],bottom=[i[0] for i in final_list], color='blue',label="lincRNA",alpha=0.5,width=barwidth,edgecolor="white")
plt.bar(r,[i[2] for i in final_list],bottom=[i[0]+i[1] for i in final_list],color='red',label="pseudogene",alpha=0.7,width=barwidth,edgecolor="white")
plt.bar(r,[i[3] for i in final_list],bottom=[i[0]+i[1]+i[2] for i in final_list],color='yellow',label="other",alpha=0.7,width=barwidth,edgecolor="white")
plt.legend(loc="upper left")
plt.ylabel("%")
plt.xticks([0.25,0.75],["ixn_IFIT2\nn="+str(len_list[0])+"_clusters","IFN_IFIT2\nn="+str(len_list[1])+"_clusters"])
plt.xlim(0,1)
plt.ylim(0,1)
plt.tick_params(direction="out")
plt.savefig(sys.argv[-1]+".svg")


		

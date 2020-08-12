#!/usr/bin python
_author_ = "Mitch Ledwith"
"""
arguments count_CLIP_meta.py intersect_file.bed out_file
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
			tmp_list.append(new_line[11])
		blah_blah = Counter(tmp_list)
	len_list.append(len(tmp_list))
	final_list.append([blah_blah["five_prime_utr"]/float(len(tmp_list)),blah_blah["CDS"]/float(len(tmp_list)),blah_blah["three_prime_utr"]/float(len(tmp_list))])

barwidth = 0.4
r = [0.25,0.75]
plt.bar(r,[i[0] for i in final_list],color='gray',label="5'UTR",alpha=0.5,width=barwidth,edgecolor="white")
plt.bar(r,[i[1] for i in final_list],bottom=[i[0] for i in final_list], color='blue',label="CDS",alpha=0.5,width=barwidth,edgecolor="white")
plt.bar(r,[i[2] for i in final_list],bottom=[i[0]+i[1] for i in final_list],color='red',label="3'UTR",alpha=0.7,width=barwidth,edgecolor="white")
plt.legend(loc="upper left")
plt.ylabel("%")
plt.xticks([0.25,0.75],["ixn_IFIT2\nn="+str(len_list[0])+"_clusters","IFN_IFIT2\nn="+str(len_list[1])+"_clusters"])
plt.xlim(0,1)
plt.ylim(0,1)
plt.tick_params(direction="out")
plt.savefig(sys.argv[-1]+".svg")
plt.show()


		

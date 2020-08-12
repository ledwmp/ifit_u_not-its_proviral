import sys
from collections import Counter
import numpy as np 
import matplotlib.pyplot as plt
import scipy.stats as stats
#CLIP1,CLIP2,SMInput1,SMInput2,CLIP1,CLIP2

coverage_dict = {}
with open("/media/mitch/TopSeqret_NTFS/CLIP_Lib1_16/fastq/Alignments/pcr_dedup/IFIT2_CLIP/coverage/coverage_total.txt") as r:
	for line in r:
		coverage_dict[line.split("\t")[0]] = float(line.split("\t")[1].strip())

def calc_coverage(input_file):
	tmp_counter = Counter()
	with open(input_file) as r:
		for line in r:
			tmp_counter["\t".join(line.split("\t")[0:6])] += float(line.split("\t")[7].strip())
	r.close()
	sort_tmp = [(i,tmp_counter[i]) for i in sorted(tmp_counter.keys())]
	return sort_tmp
total_list = []

for item in sys.argv[1:-1]:
	cov_file = calc_coverage(item)
	cov_key = item.split("/")[-1].split("_")[0]
	file_norm = [(i[0],i[1]/coverage_dict[cov_key]) for i in cov_file]
	total_list.append(file_norm)

pearson_list = []

for item in total_list:
	for thing in total_list:
		pearsonsr,pvalue = stats.pearsonr([i[1] for i in item],[i[1] for i in thing])
		pearson_list.append(pearsonsr**2)

print pearson_list

pearson_array = np.reshape(pearson_list, [len(total_list),len(total_list)])
print pearson_array

	

plt.figure(figsize=(5,5))

fig,ax = plt.subplots()
im = ax.pcolormesh(np.linspace(0,1.2,len(total_list)+1),np.linspace(0,1.2,len(total_list)+1),pearson_array,cmap="Reds")
fig.colorbar(im)

name_list = ["ixnCLIP1","ixnCLIP2","ixnSMIn1","ixnSMIn2","IFNCLIP1","IFNCLIP2"]
#plt.tick_params(length=0)
ax.set_xticks(np.linspace(0.1,1.1,len(total_list)))
ax.set_yticks(np.linspace(0.1,1.1,len(total_list)))
ax.set_xticklabels(name_list,fontsize=8)
#ax.set_yticklabels(name_list,rotation="vertical",fontsize=8)
ax.set_yticklabels(name_list,fontsize=8)
ax.tick_params(length=0)
ax.axis("tight")
ax.yaxis.set_label_position("right")
plt.ylabel("(Pearson's R)^2",rotation=270,labelpad=90)
plt.xlabel("Read Density in ixn Clusters")
plt.savefig(sys.argv[-1]+".svg")

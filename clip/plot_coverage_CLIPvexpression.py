import sys
from collections import Counter
import numpy as np 
import matplotlib.pyplot as plt
import scipy.stats as stats

#11	62841604	62841776	.	1000	-	10	1671
coverage_dict = {}
with open("/media/mitch/TopSeqret_NTFS/CLIP_Lib1_16/fastq/Alignments/pcr_dedup/IFIT2_CLIP/coverage/coverage_total.txt") as r:
	for line in r:
		coverage_dict[line.split("\t")[0]] = float(line.split("\t")[1].strip())
r.close()
		

def calc_coverage(input_file):
	tmp_counter = Counter()
	with open(input_file) as r:
		for line in r:
			tmp_counter["\t".join(line.split("\t")[0:6])] += float(line.split("\t")[7].strip())
	r.close()
	sort_tmp = [(i,tmp_counter[i]) for i in sorted(tmp_counter.keys())]
	return sort_tmp

def norm(input_file):
	tmp_file = calc_coverage(input_file)
	tmp_key = input_file.split("/")[-1].split("_")[0]
	tmp_file_norm = [(i[0],i[1]/coverage_dict[tmp_key]) for i in tmp_file]
	return tmp_file_norm

def fold(list_clip,list_input):
	zip_lists = zip(list_clip,list_input)
	fold_list = [(a[0],a[1],b[1]) for a,b in zip_lists]
	return fold_list

a_file = norm(sys.argv[1]) #CLIP_1
b_file = norm(sys.argv[2]) #CLIP_2
c_file = norm(sys.argv[3]) #Input_1
d_file = norm(sys.argv[4]) #Input_2

CLIPs = [(i[0],(i[1]+j[1])/2.0) for i,j in zip(a_file,b_file)]
inputs = [(i[0],(i[1]+j[1])/2.0) for i,j in zip(c_file,d_file)]

FC = fold(CLIPs,inputs)
zip_CLIPs_FC = zip(FC,CLIPs)

filter_FC = [(a[0],a[1]/a[2]) for a,b in zip_CLIPs_FC if a[2] > 0.0 and b[1] > 0.0]
filter_CLIPs = [(a[0],b[1]) for a,b in zip_CLIPs_FC if a[2] > 0.0 and b[1] > 0.0]

for i,j in zip(filter_FC,filter_CLIPs):
	if np.log2(i[1]) >= 2.0 and np.log2(j[1]) >= -18.0: #ixn cut at 2 and -18
		print i[0]

plt.scatter([np.log2(i[1]) for i in filter_CLIPs],[np.log2(i[1]) for i in filter_FC],s=5,color='k',alpha=0.1) 

plt.ylabel('log2(CLIP_ave/SMIn_ave)')
plt.xlabel('log2(CLIP_peak_density)')
plt.tick_params(direction="out")
plt.axhline(y=2,color="r")
plt.axvline(x=-18,color="r")
plt.show()




	


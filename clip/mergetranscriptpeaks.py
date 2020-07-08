#!/usr/bin python
_author_ = "Mitch Ledwith"

"""
arguments mergetranscriptpeaks.py file1,file2,file3...
CLIPPER output is transcript-based, so need to convert transcript-based CLIPPER output to genomic-based gene output
output of inputnorm.bed is not sorted. This script will re-sort input norm bed files first
"""



import sys
print sys.argv
import subprocess

def find_overlap(peak1,peak2,sortcollapsebed):
	peak1_split = peak1.split("\t")
	peak2_split = peak2.split("\t")
	if peak1_split[0] == peak2_split[0]:
		if (int(peak1_split[2]) > int(peak2_split[1])) and (peak1_split[12] == peak2_split[12]): #overlap and on same strand
			#then these peaks overlap, keep only larger FC as newpeak1. Discard smaller FC. Read in next line as newpeak2.
			if float(peak1_split[10]) >= float(peak2_split[10]):
				newpeak1 = peak1 #keep peak1 as newpeak1
			else:
				newpeak1 = peak2 #reset peak2 as newpeak1	
		else:
			sortcollapsebed.write(peak1) #write peak1 to file
			newpeak1 = peak2
			#these peaks do not overlap, write peak1 to file. Reset peak2 as newpeak1. Read in next line as newpeak2.  	
	else:
		sortcollapsebed.write(peak1)
		newpeak1 = peak2
		#write peak1 to file. Reset peak2 to newpeak1 and read next line as newpeak2. 
	return newpeak1

for item in sys.argv[1:]:
	sorted_item = item.strip("_inputnorm.bed")+"_inputnormsort.bed"
	collapsed_item = item.strip("_inputnorm.bed")+"_inputnormsortcollapse.bed"
	cmd_line = "sort -k 1,1 -k2,2n "+item+" > "+sorted_item
	child=subprocess.call(cmd_line,shell=True)
	with open(sorted_item) as r:
		write_file = open(collapsed_item,"a")
		old_line = next(r)
		for line in r:
			out_line = find_overlap(old_line,line,write_file)
			old_line = out_line
		write_file.close()
	r.close()
			
			
				
			
			

#!/usr/bin python
_author_ = "Mitch Ledwith"

"""
arguments merge_adjacent_peaks.py bed1,bed2,...
CLIPPER overfits our peaks like crazy, this connects adjacent peaks

"""
import sys

def merge_peaks(peak1,peak2):
	new_peak = []
	peak1_split = peak1.split("\t") #first peak
	peak2_split = peak2.split("\t") #second peak
	if float(peak1_split[4]) >= float(peak2_split[4]): #pvalues
		name = peak2_split[3]
		pval = peak2_split[4]
	else:
		name = peak2_split[3]
		pval = peak1_split[4]
	new_peak.append(peak1_split[0])
	new_peak.append(peak1_split[1])
	new_peak.append(peak2_split[2])
	new_peak.append(name)
	new_peak.append(pval)
	new_peak.append(peak1_split[5])
	return "\t".join(new_peak)
	
def iterate_peaks(peak1,peak2,write_file): #peak1 = oldpeak, peak2 = newpeak
	peak1_split = peak1.split("\t")
	peak2_split = peak2.split("\t")
	if peak1_split[0] == peak2_split[0]: #they're on the same chromosomes
		if (int(peak1_split[2]) == int(peak2_split[1])) and (peak1_split[5] == peak2_split[5]): #overlap and on same strand, cuz of 0-based indexing
			newpeak1 = merge_peaks(peak1,peak2)	
			print newpeak1
		else:
			write_file.write("\t".join(peak1.split("\t")[0:6])+"\n") #write peak1 to file
			newpeak1 = peak2
			#these peaks do not overlap, write peak1 to file. Reset peak2 as newpeak1. Read in next line as newpeak2.  	
	else:
		write_file.write("\t".join(peak1.split("\t")[0:6])+"\n")
		newpeak1 = peak2
		#write peak1 to file. Reset peak2 to newpeak1 and read next line as newpeak2. 
	return newpeak1


for item in sys.argv[1:]:
	merge_file = open(item.strip(".bed")+"_filter.bed","a")
	with open(item) as r:
		old_line = next(r)
		for line in r:
			out_line = iterate_peaks(old_line,line,merge_file)
			old_line = out_line
	merge_file.close()
	r.close()
			
	

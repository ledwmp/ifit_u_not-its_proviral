#!/usr/bin python
_author_ = "Mitch Ledwith"

"""
arguments scan_coverage_for_peak.py depth_file
where depth file is output of concat_coverage.py
but,looks like: ensemblid\tstrand\tcsv of read depth along transcript for each replicate
This script takes a concatenated coverage file and finds positions of putative pausing on transcript based on shape of an RPF
"""
import sys
import numpy as np
import matplotlib.pyplot as plt 
import window_func
from collections import defaultdict
from scipy.signal import find_peaks
from operator import itemgetter
from itertools import *
			
def string_to_bed12(bed_list,strand):
	tmp_list = []
	for k,g in groupby(enumerate(bed_list), lambda (i,x):i-x):
		tmp_list.append(map(itemgetter(1),g))
	tmp_list = [i for i in tmp_list if i[-1]-i[0] > 0]
	start = tmp_list[0][0]
	end = tmp_list[-1][-1]+1
	blocks = len(tmp_list)
	lengths = ",".join([str(i[-1]-i[0]+1) for i in tmp_list])
	starts = ",".join([str(i[0]-start) for i in tmp_list])
	return str(start),str(end),str(blocks),lengths,starts
	
	
bed_dict_int = defaultdict(list)
bed_dict_chrom = {}
bed_dict_strand = {}
with open("Homo_sapiens.GRCh38.95_CDS_filter_exp_collapse_WSN.bed") as r:
	for line in r:
		new_line = line.split("\t")
		bed_dict_int[new_line[3]].append(range(int(new_line[1]),int(new_line[2])))
		bed_dict_chrom[new_line[3]] = new_line[0]
		bed_dict_strand[new_line[3]] = new_line[5]
r.close()

rep_list = [(0,6),(1,7),(2,8),(3,9),(4,10),(5,11)]
comp_list = [(0,1),(0,2),(1,2)]

write_file = open("pause.bed","a")


with open(sys.argv[1]) as r:
	for line in r:
		new_line = line.strip().split("\t")
		ensembl_id = new_line[0]
		strand = new_line[1]
		print ensembl_id
		if strand == "+":
			coverage = [[float(j) for j in i.split(",")] for i in new_line[2].strip(":").split(":")]
			bed_array = np.array([i for j in bed_dict_int[ensembl_id] for i in j])
		else:
			coverage = [[float(j) for j in i.split(",")][::-1] for i in new_line[2].strip(":").split(":")]
			bed_array = np.array([i for j in bed_dict_int[ensembl_id] for i in j][::-1])
		coverage_array = np.array(coverage)
		median_array = np.empty([len(coverage_array[0]),])
		for x in range(0,24):
			temp_array = coverage_array[x]/np.nanmedian(coverage_array[x])
			median_array = np.vstack((median_array,temp_array))
		median_array = median_array[1:]
		log_array = np.log2(median_array)[12:]

		new_log = log_array
		not_nan = new_log[:,np.isfinite(log_array).any(axis=0)]
		filtered_coverage = coverage_array[:,np.isfinite(log_array).any(axis=0)] #filtered based on log transformation of ribo array	
		filtered_median = median_array[:,np.isfinite(log_array).any(axis=0)]
		filtered_bed = bed_array[np.isfinite(log_array).any(axis=0)]
		divide = float(len(not_nan[0]))/float(len(log_array[0]))
		if divide > 0.75: #more than 75% of the sequence is covered by reads
			average = np.mean(np.mean(filtered_coverage[12:],axis=1)) 
			if average > 8: #average of 8 reads piled at each position
				len_array = len(filtered_coverage[0])
				delta_slope_array = np.empty([len_array,])
				for z in range(12,24):
					slope = []
					delta_slope = []
					med = []
					for x in range(0,len(filtered_coverage[0])):
						#tmp.append(window_func.average_normalize(filtered_coverage,z,x,135))
						med.append(window_func.pick_median(filtered_coverage[z],x,25))
						#slope.append(window_func.sliding_slope(filtered_coverage,z,x,5))
					for x in range(0,len(filtered_coverage[0])):
						slope.append(window_func.sliding_slope(med,"blah",x,5))
					for x in range(0,len(filtered_coverage[0])):
						delta_slope.append(window_func.delta_slides(slope,x,50,28))
					delta_slope_array = np.vstack((delta_slope_array,np.array(delta_slope).reshape(len_array,)))
			
		
				delta_slope_array = delta_slope_array[1:]

				orig_peaks = []
				comb_peaks = []
				
				for x in rep_list:
					peaks_a,properties_a = find_peaks(delta_slope_array[x[0]],distance=35,prominence=(1.0,None))
					peaks_b,properties_b = find_peaks(delta_slope_array[x[1]],distance=35,prominence=(1.0,None))
					orig_peaks.append((peaks_a,peaks_b))
					overlap_peaks = window_func.overlap_windows(peaks_a,peaks_b,17,17)
					comb_peaks.append(overlap_peaks)
				
				sup_peaks = []
				fun_peaks = []
				for x in comp_list:
					sup_peaks.append(comb_peaks[x[0]]+comb_peaks[x[1]])
					fun_peaks.append(window_func.combine_peak_make_windows(comb_peaks[x[0]],comb_peaks[x[1]],17,17))

				for item in fun_peaks[1]: #list of lists of peaks
					chrom = bed_dict_chrom[ensembl_id]
					if strand == "+":
						peak_coord = list(filtered_bed[item])
						#need to make bed12 out of list of coordinates
						bed_info = string_to_bed12(peak_coord,strand)
						write_file.write("\t".join([chrom]+[bed_info[0]]+[bed_info[1]]+[ensembl_id]+["."]+[strand]+[bed_info[0]]+[bed_info[1]]+["."]+[bed_info[2]]+[bed_info[3]]+[bed_info[4]])+"\n")

						
					
					else:
						peak_coord = list(filtered_bed[item])[::-1]
						#peak_coord = list(filtered_bed[item])
						bed_info = string_to_bed12(peak_coord,strand)
						write_file.write("\t".join([chrom]+[bed_info[0]]+[bed_info[1]]+[ensembl_id]+["."]+[strand]+[bed_info[0]]+[bed_info[1]]+["."]+[bed_info[2]]+[bed_info[3]]+[bed_info[4]])+"\n")

					
write_file.close()							
r.close()


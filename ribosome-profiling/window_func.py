#!/usr/bin python
_author_ = "Mitch Ledwith"
"""
a bunch of auxiliary functions to scan_coverage_for_peak.py
"""

from scipy import stats
import numpy as np

def find_median(tmp_list):
	len_list = len(tmp_list)
	if len_list % 2 == 1:
		return sorted(tmp_list)[len_list/2]
	else:
		return sum(sorted(tmp_list)[(len_list/2)-1:(len_list/2)+1])/2.0

def pick_median(nuc_list,current_index,median_size):
	len_list = len(nuc_list)
	if current_index < (median_size/2): #first n nucleotides
		return find_median(nuc_list[0:current_index+(median_size/2)+1])
	elif current_index > (len_list - median_size/2): #last n nucleotides
		return find_median(nuc_list[(current_index-(median_size/2)):len_list])
	else: #middle
		return find_median(nuc_list[(current_index-(median_size/2)):current_index+(median_size/2)+1])

def average_normalize(master_list,sample_number,current_index,window_size):
	cov_list = master_list[sample_number]	
	len_list = len(cov_list)
	current_coverage = float(cov_list[current_index])
	if current_index < (window_size/2): #first n nucleotides
		return current_coverage/(sum(cov_list[0:current_index+(window_size/2)+1])/float(current_index+(window_size/2)))
	elif current_index > (len_list - window_size/2):
		return current_coverage/(sum(cov_list[(current_index-(window_size/2)):len_list])/float(len_list-(current_index-(window_size/2))))
	else:
		return current_coverage/(sum(cov_list[(current_index-(window_size/2)):current_index+(window_size/2)+1])/float(window_size))

def sliding_slope(master_list,sample_number,current_index,slope_size):
	if sample_number == "blah":
		cov_list = master_list
	else:	
		cov_list = master_list[sample_number]
	len_list = len(cov_list)
	if current_index < (slope_size/2):
		y = cov_list[0:current_index+(slope_size/2)+1]
		x = range(0,len(y))
	elif current_index > (len_list-slope_size/2):
		y = cov_list[(current_index-(slope_size/2)):len_list]
		x = range(0,len(y))
	else:
		y = cov_list[(current_index-(slope_size/2)):current_index+(slope_size/2)+1]
		x = range(0,len(y))	
	slope,intercept,r_value,p_value,std_err = stats.linregress(x,y)
	return slope

def sliding_slope_nozero(master_list,sample_number,current_index,slope_size):
	if sample_number == "blah":
		cov_list = master_list
	else:	
		cov_list = master_list[sample_number]
	len_list = len(cov_list)
	if current_index < (slope_size/2):
		y = cov_list[0:current_index+(slope_size/2)+1]
		x = range(0,len(y))
	elif current_index > (len_list-slope_size/2):
		y = cov_list[(current_index-(slope_size/2)):len_list]
		x = range(0,len(y))
	else:
		y = cov_list[(current_index-(slope_size/2)):current_index+(slope_size/2)+1]
		x = range(0,len(y))
	if 0.0 in y:
		return np.nan
	else:
		slope,intercept,r_value,p_value,std_err = stats.linregress(x,y)
		return slope

def delta_slides(master_list,current_index,window_size,buffer_size):
	#lets set window_size ~ 50 usually, buffer_size = 28 usually	
	slope_list = master_list
	len_list = len(slope_list)
	if current_index < (window_size/2):
		return 0.0
	elif current_index > (len_list-window_size/2):
		return 0.0
	else:
		search_area = (window_size-buffer_size)/2
		tmp_up = slope_list[current_index-window_size/2:current_index-search_area] 
		tmp_down = slope_list[current_index+search_area:current_index+window_size/2+1]
		max_list = max(tmp_up)
		min_list = min(tmp_down)
		if max_list > 0.0 and min_list < 0.0:
			return max_list - min_list
		else:
			return 0.0

def overlap_windows(list_a,list_b,shift_l,shift_r): #where list_a/b are lists of coordinates
	tmp_list = []
	for i_a in list_a:
		for i_b in list_b:
			if abs(i_a - i_b) < (shift_l+shift_r):
				tmp_list.append((i_a+i_b)/2)
	return tmp_list

def combine_peak_make_windows(list_a,list_b,shift_l,shift_r):
	tot_list = sorted(list_a+list_b)
	tot_list = [range(i-shift_l,i+shift_r) for i in tot_list]
	tmp_list = []
	search = range(0,1)
	for item in tot_list:
		query = item
		if search[-1] > query[0]:
			search = sorted(list(set(search+query)))
	
		else:
			if len(search) > 1:
				tmp_list.append(search)
				search = query
			else:
				search = query
	if len(tot_list) > 0:
		tmp_list.append(search)
		return tmp_list
	else:
		return []
	
def sum_over_window(master_list,sample_number,current_index,window_size):
	cov_list = master_list[sample_number]	
	len_list = len(cov_list)
	current_coverage = float(cov_list[current_index])
	if current_index < (window_size/2): #first n nucleotides
		return sum(cov_list[0:current_index+(window_size/2)+1])
	elif current_index > (len_list - window_size/2):
		return sum(cov_list[(current_index-(window_size/2)):len_list])
	else:
		return sum(cov_list[(current_index-(window_size/2)):current_index+(window_size/2)+1])
				 




 

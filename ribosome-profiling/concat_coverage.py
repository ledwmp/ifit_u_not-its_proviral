#!/usr/bin python
_author_ = "Mitch Ledwith"

"""
arguments concat_coverage.py file1.depth,file2.depth,file3.depth,...
As it turns out, bedtools coverage of a bam over a bed produces a massive and difficult to read file. This script takes the output of many bedtools coverages and concatenates it together to make it easier to parse. 
"""
import sys
import numpy as np
import itertools
from collections import defaultdict

#read in order of names
write_file = open("/".join(sys.argv[1].split("/")[:-1])+"_CDScoverage.depth","a")

t = ()
for x in range(1,25):
	exec("open_%s=open('%s','r')" % (x,sys.argv[x]))
	exec("t = t + (open_%s,)" % (x))

last_id = ""
temp_list = defaultdict(list)
	
for i in itertools.izip(*t):
	ensembl_id = i[0].strip().split("\t")[3]
	strand = i[0].strip().split("\t")[5]
	if ensembl_id != last_id:
		if len(temp_list) > 0:
			#print len(temp_list[0])
			write_file.write(last_id+"\t"+last_strand+"\t")
			for x in range(0,24):
				write_file.write(",".join(temp_list[x])+":")
			write_file.write("\n")
			temp_list = defaultdict(list)
			for x in range(0,24):			
				temp_list[x].append(i[x].strip().split("\t")[7])
		else:
			for x in range(0,24):			
				temp_list[x].append(i[x].strip().split("\t")[7])
	else:
		for x in range(0,24):	
			temp_list[x].append(i[x].strip().split("\t")[7])
	last_id = ensembl_id
	last_strand = strand
write_file.write(last_id+"\t"+strand+"\t")
for x in range(0,24):
	write_file.write(",".join(temp_list[x])+":")
write_file.write("\n")
write_file.close()
for x in range(1,25):
	exec("open_%s.close()" % (x))
		



#!/usr/bin python
_author_ = "Mitch Ledwith"
"""
arguments split_CLIP_meta.py intersect_file.bed
"""
import sys



for item in sys.argv[1:]:
	open_3UTR = open(item.strip(".bed")+"_3UTR.bed","a")
	open_CDS = open(item.strip(".bed")+"_CDS.bed","a")
	open_5UTR = open(item.strip(".bed")+"_5UTR.bed","a")
	with open(item) as r:
		for line in r:
			new_line = line.split("\t")
			if "three" in line:
				open_3UTR.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (new_line[0],new_line[1],new_line[2],new_line[10],new_line[11],new_line[5]))
			if "five" in line:
				open_5UTR.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (new_line[0],new_line[1],new_line[2],new_line[10],new_line[11],new_line[5]))
			if "CDS" in line:
				open_CDS.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (new_line[0],new_line[1],new_line[2],new_line[10],new_line[11],new_line[5]))
	r.close()
	open_3UTR.close()
	open_5UTR.close()
	open_CDS.close()

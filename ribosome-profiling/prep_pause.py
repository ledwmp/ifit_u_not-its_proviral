#!/usr/bin python
_author_ = "Mitch Ledwith"

"""
arguments prep_pause.py depth_file1 depth_file2
where depth file is output of bedtools coverage against pause bed
"""

import sys
import numpy as np

out_file = sys.argv[1].split("/")[0]+sys.argv[-1]+".txt"

name_list = ["wt_mo1","wt_IFN1","wt_ixn1","D8_mo1","D8_IFN1","D8ixn1","wt_mo2","wt_IFN2","wt_ixn2","D8_mo2","D8_IFN2","D8_ixn2"]

column_number = 12 #count pauses

#1	33036742	33036776	ENSG00000004455	.	-	33036742	33036776	.	1	34	0	27

start,stop,gene = np.loadtxt(sys.argv[1],dtype=str,delimiter="\t",usecols=(1,2,3),unpack = True)
gene_id = [i+":"+j+"-"+k for i,j,k in zip(gene,start,stop)]
gene_id = np.concatenate((np.array(["Sample"]),np.array(gene_id))).reshape(-1,1)

for item,name in zip(sys.argv[1:-1],name_list):
	sample_1 = np.loadtxt(item,dtype=float,delimiter="\t", usecols=(column_number,))
	sample_1 = np.concatenate((np.array([name]),sample_1)).reshape(-1,1)
	gene_id = np.hstack((gene_id,sample_1))
print gene_id
print gene_id.shape

np.savetxt(out_file,gene_id,delimiter="\t",fmt='%s')

	

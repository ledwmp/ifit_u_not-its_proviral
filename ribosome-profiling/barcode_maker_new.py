#!/usr/bin python
_author_ = "Mitch Ledwith"
"""
adding this because it's a fun little problem
lets make 24 unique 5-nt barcodes, out of 1024 options
all requiring at least 2 mutations to get to another barcode (i.e. one mutation won't do it)
all between 25-75% GC ???? (no homopolymers) (avoid 3-identical bases in a row)
"""


from itertools import product
from itertools import combinations
from collections import Counter
import random

def onebpdifference(seq):
	base_list = ["A","T","C","G"]
	temp_list = []
	for x in range(0,len(seq)):
		for item in base_list:
			my_new_creation = seq[:x] + item + seq[x+1:]
			if my_new_creation != seq:
				temp_list.append(my_new_creation)
	return temp_list

def allperm(A_list,T_list,C_list,G_list):
	allproductset = [g for g in product([i for i in combinations(A_list,3)],[i for i in combinations(T_list,3)],[i for i in combinations(C_list,3)],[i for i in combinations(G_list,3)])]
	print "okay"
	return set(allproductset)
			
bases = "ATCG"
permutations = ["".join(i) for i in product(bases,repeat=5)]

new_perm_list = []

for item in permutations:
	if "AAA" not in item and "GGG" not in item and "CCC" not in item and "TTT" not in item: #gets rid of homopolymers
		my_count = Counter(item)
		if 1 < (float(my_count["A"])+float(my_count["T"])) < 4: #gets rid of anything that is < 40% GC or >60% GC
			new_perm_list.append(item)

random.Random(6).shuffle(new_perm_list) #needs to be shuffled, otherwise favors the first members of the list
one_base_library = []
random.seed() #reseed random

for item in new_perm_list:
	new_list = onebpdifference(item)
	if set(new_list).isdisjoint(set(one_base_library)) == True:
		one_base_library.append(item)

barcode_dict = {}
for x in range(0,5):
	A_list = []
	T_list = []
	C_list = []
	G_list = []
	for item in one_base_library:
		if item[x] == "A":
			A_list.append(item)
		if item[x] == "T":
			T_list.append(item)
		if item[x] == "C":
			C_list.append(item)
		if item[x] == "G":
			G_list.append(item)	
	print "yes"
	allpermutations = allperm(A_list,T_list,C_list,G_list)
	barcode_dict["Position_"+str(x)] = allpermutations
print barcode_dict
#for item in barcode_dict.values()[1:]:
	#super_set = 
	






	

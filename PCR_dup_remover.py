import sys
from collections import Counter
#arguments PCR_dup_remover.py logfile,file1,file2,file3,
#7nt random barcode, 16384 unique possibilities
#first build list 
#input a samtools view -F 4 sam file

#samflags F:0, R:16, F,secondary_alignment:256, R,secondary_alignment:272
out_log = open(sys.argv[1],"a")
def pcr_collapser(input_sam):
	output_sam = open(input_sam.strip(".sam")+"_pcrdedup.sam","a")
	out_log.write(input_sam+"\n")
	all_fragments = Counter()
	with open(input_sam) as r:
		for line in r:
			if line[0] == "@":
				output_sam.write(line)
				continue
			else:
				#readname,strand=FLAG,chromosome,startcoordinate
				new_line = line.split("\t")
				if new_line[1] in ("0","256"):
					seq_start = new_line[3]
				elif new_line[1] in ("16","272"):
					seq_start = str(int(new_line[3]) + len(new_line[9]))	
				else:
					continue
				umi = new_line[0].split("umi:")[1]
				all_fragments[str(umi+"_"+new_line[1]+"_"+new_line[2]+"_"+seq_start)] += 1
	single_use_set = set(dict(all_fragments).keys())
	deduplicated_umis = 0
	removed_umis = 0
	with open(input_sam) as r:
		for line in r:
			if line[0] == "@":
				continue
			else:
				new_line = line.split("\t")
				if new_line[1] in ("0","256"):
					seq_start = new_line[3]
				elif new_line[1] in ("16","272"):
					seq_start = str(int(new_line[3]) + len(new_line[9]))	
				else:
					continue
				umi = new_line[0].split("umi:")[1]
				new_set = set([str(umi+"_"+new_line[1]+"_"+new_line[2]+"_"+seq_start)])
				if new_set.isdisjoint(single_use_set) == False: #if the umi hasn't been used yet, returns False
					output_sam.write(line)
					single_use_set.discard(list(new_set)[0]) #retires the umi
					deduplicated_umis = deduplicated_umis+1
				else:	#else True, or the sets are disjoint, and the umi has already been used
					removed_umis = removed_umis+1
	out_log.write("Deduplicated_umis:"+str(deduplicated_umis)+"\n"+"Removed_umis:"+str(removed_umis)+"\n")
	print "Deduplicated_umis:"+str(deduplicated_umis)
	print "Removed_umis:"+str(removed_umis)

for item in sys.argv[2:]: #WHO IS JACK SPRAT?
	pcr_collapser(item)
	
print "I'm done! :)" #if this works first try, take a shot
				
		
	

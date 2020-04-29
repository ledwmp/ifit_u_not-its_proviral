import sys
#arguments debarcody.py fastq_file barcode_file
print sys.argv

sample_index = {}

with open(sys.argv[2]) as r:
	for line in r:
		new_line = line.split("\t")
		sample_index[new_line[2].strip()] = new_line[1].strip()
		
print sample_index
for key,value in sample_index.iteritems():
	exec("open_%s=open('/media/mitch/TopSeqret_NTFS/CLIP_Lib1_16/fastq/%s.fastq','a')" % (key,value))

write_none = open("/media/mitch/TopSeqret_NTFS/CLIP_Lib1_16/fastq/nobarcode.fastq",'a') #bad barcodes
write_log = open("/media/mitch/TopSeqret_NTFS/CLIP_Lib1_16/fastq/debarcode.log","a") #log

result_dict = {}

for item in sample_index.keys():
	result_dict[item] = 0.0 #populate a barcode counter
result_dict["Unidentified_barcode"] = 0.0 #only the sequencer knows

with open(sys.argv[1]) as r:
	for line in r:
		read_name = line.strip().split(" ")[0]
		read = next(r).strip()
		qual_name = next(r).strip()
		qual = next(r).strip()
		current_bc = read[5:9]
		current_umi = read[0:5]+read[9:11]
		if read[5:9] in sample_index.keys():
			exec("open_%s.write('%s_%s')" % (current_bc,read_name,"umi:"+current_umi))
			exec("open_%s.write('\\n')" % (current_bc))
			exec("open_%s.write('%s')" % (current_bc,read[11:]))
			exec("open_%s.write('\\n')" % (current_bc))
			exec("open_%s.write('%s')" % (current_bc,qual_name))
			exec("open_%s.write('\\n')" % (current_bc))
			exec("open_%s.write('%s')" % (current_bc,qual[11:]))
			exec("open_%s.write('\\n')" % (current_bc))
			result_dict[current_bc] = result_dict[current_bc] + 1.0 
				
		else:
			#print "CAUTION: unbarcoded read identified"
			"write_none.write(%s_%s)" % (read_name,"umi:"+current_umi)
			write_none.write('\n')
			"write_none.write('%s')" % (read[11:])
			write_none.write('\n')
			"write_none.write('%s')" % (qual_name)
			write_none.write('\n')
			"write_none.write('%s')" % (qual[11:])
			write_none.write('\n')
			result_dict["Unidentified_barcode"] = result_dict["Unidentified_barcode"]+ 1.0

for key,value in result_dict.iteritems():
	write_log.write(sample_index[key] + "\t" + key + "\t" + str(value) + "\n")

write_log.close()

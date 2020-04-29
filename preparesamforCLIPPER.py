#arguments preparesamforCLIPPER.py file1,file2,file3,protoheader_new.sam
#where protoheader_new.sam is sam header with chromosome names that match CLIPPER
#requires samtools
import sys
import subprocess
print sys.argv

for item in sys.argv[1:-1]:
	cmd_line = "samtools view -bS "+item+" | samtools sort -o"+item.strip(".sam")+"_sorted.bam" #sam to sorted bam
	child = subprocess.call(cmd_line,shell=True)
	cmd_line = "samtools view -H "+item+" > "+item.strip(".sam")+"_header.sam" #just the header
	child = subprocess.call(cmd_line,shell=True)
	with open(item.strip(".sam")+"_header.sam") as r:
		for line in r:
			if line[0:2] == "@PG":
				new_line = line
			else:
				continue
	new_header = open(item.strip(".sam")+"_header.sam","w")
	with open(sys.argv[-1]) as r:
		for line in r:
			new_header.write(line)
	new_header.write(new_line)
	cmd_line = "samtools reheader "+item.strip(".sam")+"_header.sam"+" "+item.strip(".sam")+"_sorted.bam"+" > "+item.strip(".sam")+"_sorted_reheader.bam"
	child = subprocess.call(cmd_line,shell=True)



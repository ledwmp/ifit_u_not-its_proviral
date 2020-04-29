import sys
import subprocess
from collections import defaultdict
from scipy.stats import chi2_contingency,fisher_exact
from collections import Counter
from math import log
#this script relies on bedtools intersect to determine reads in CLIP and input files that overlap with bedpeak file
#arguments input_norm.py file1,file2,file3,file4,file5,file6,...
#where file1=CLIP1bed,file2=input1bed,file3=bedpeak1,file4=CLIP2bed,file5=input2bed,file6=bedpeak2...

#CLIPPER bed8: chromosome, genomic_start, genomic_stop, cluster_name, min_pval, strand, thick_start, thick_stop
#now a bed6, after merge transcript peaks
#bamtobed bed6: chromosome, genomic_start, genomic_stop, readname, qual_score, strand

grouped_files = zip(sys.argv[1::3],sys.argv[2::3],sys.argv[3::3]) #zipped (CLIP_file,input_file,bedpeak_file)

for item in grouped_files:
	cmd_line = "bedtools intersect -a "+item[2]+" -b "+item[0]+" "+item[1]+" -wa -wb -s -names CLIP input > "+item[0].strip(".bed")+"_overlappeaks.bed"
	child = subprocess.call(cmd_line,shell=True)

	CLIP_reads = float(sum(1 for line in open(item[0])))
	input_reads = float(sum(1 for line in open(item[1])))

	peak_dict = defaultdict(list)
	with open(item[0].strip(".bed")+"_overlappeaks.bed") as r:
		for line in r:
			original_peak = ",".join(line.split("\t")[:6])
			overlap_read = line.split("\t")[6]
			peak_dict[original_peak].append(overlap_read)
#chromosome	start	end	name (colon separated region)	reads in CLIP	reads in INPUT	p-value	chi value or (F)isher	(F)isher or (C)hi square test	enriched or depleted	negative log10p value (400 if above certain threshold)	log2 fold change	entropy	

#2by2 contingency table: 	
				#	peak	notpeak
				#input	|a	|b
				#CLIP	|c	|d
	out_file = open(item[0].strip(".bed")+"_inputnorm.bed","a")
	for key,value in peak_dict.iteritems():
		count_reads = Counter(value)
		a = float(count_reads["input"]) #input reads in peak
		c = float(count_reads["CLIP"]) #CLIP reads in peak
		b = input_reads-a #input reads not in peak
		d = CLIP_reads-c #CLIP reads not in peak
		if a >= 1:
			fold_over_input = float(c)/float(a)	#Fold-change
			entropy = (float(c)/CLIP_reads)*log((float(c)/CLIP_reads)/(float(a)/input_reads),2) #from biorxiv paper, IPi * log2(IPi/Qi) where IPi = fraction IP reads in region and Qi = fraction input reads in region
		else:
			fold_over_input = float(c)/(float(a)+1) #FC,throws a 1 in denominator if no input reads detected in peak
			entropy = (float(c)/CLIP_reads)*log((float(c)/CLIP_reads)/((float(a)+1)/input_reads),2) #from biorxiv paper, IPi * log2(IPi/Qi) where IPi = fraction IP reads in region and Qi = fraction input reads in region 
		log2fc = log(fold_over_input,2)	#log2FC
		contin_array = [[a,b],[c,d]]
		if a>=1 and c>=1:
			statistic,pvalue,dof,exp = chi2_contingency(contin_array, correction=True) #chi2 test with Yates' correction
			test_type = "CHI2"
		else:
			statistic,pvalue = fisher_exact(contin_array) #fisher's exact test, two-sided
			test_type = "FISHER"
		if pvalue == 0.0:
			neglog10pvalue = "400.0"
		else:
			neglog10pvalue = -1.0*log(pvalue,10)
		new_key = key.split(",")
		#output = chrom	start	stop	namepeak	readCLIP	readinput	pvalue	chi2statistic	Chi2	neglog10pvalue	log2FC	entropy strand
		out_file.write(new_key[0]+"\t"+new_key[1]+"\t"+new_key[2]+"\t"+new_key[3]+"\t"+str(c)+"\t"+str(a)+"\t"+str(pvalue)+"\t"+str(statistic)+"\t"+test_type+"\t"+str(neglog10pvalue)+"\t"+str(log2fc)+"\t"+str(entropy)+"\t"+new_key[5]+"\n")
	out_file.close()

print "PEACE AND LOVE, WE'RE DONE"
			
			
			
		
			
			
		
			
	

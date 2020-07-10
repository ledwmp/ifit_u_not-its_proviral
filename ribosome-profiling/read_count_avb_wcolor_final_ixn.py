import sys
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib
import scipy.stats as stats

#rna,ribo,TE,name

gene_ida = np.loadtxt(sys.argv[1],dtype=str,delimiter="\t",skiprows=1, usecols=(0,))
gene_idb = np.loadtxt(sys.argv[2],dtype=str,delimiter="\t",skiprows=1, usecols=(0,))
sample_1 = np.loadtxt(sys.argv[1],dtype=float,delimiter="\t",skiprows=1, usecols=(1,))
sample_2 = np.loadtxt(sys.argv[2],dtype=float,delimiter="\t",skiprows=1, usecols=(1,))
gene_idc = np.loadtxt(sys.argv[3],dtype=str,delimiter="\t",skiprows=1, usecols=(0,))
pvalue = np.loadtxt(sys.argv[3],dtype=float,delimiter="\t",skiprows=1, usecols=(4,))
logfc = np.loadtxt(sys.argv[3],dtype=float,delimiter="\t",skiprows=1, usecols=(1,))
logcpm = np.loadtxt(sys.argv[3],dtype=float,delimiter="\t",skiprows=1, usecols=(2,))

dict_a = {k.strip('"'):v for k,v in zip(gene_ida,sample_1)}
dict_b = {k.strip('"'):v for k,v in zip(gene_idb,sample_2)}
dict_p = {k.strip('"'):v for k,v in zip(gene_idc,pvalue)}
dict_fc = {k.strip('"'):v for k,v in zip(gene_idc,logfc)}
dict_cpm = {k.strip('"'):v for k,v in zip(gene_idc,logcpm)}


comb_list = list(set(dict_a.keys()) & set(dict_b.keys()) & set(dict_p.keys()))

print len(comb_list)

dict_c = {k:(dict_a[k],dict_b[k],dict_p[k],dict_fc[k]) for k in comb_list if dict_cpm[k] > 5.0}

for k,v in dict_c.iteritems():
	if v[1]-v[0] > 1:
		print k.strip('"')


sizes = [24.0 if -1.0*np.log10(i[2]) > 2.0 else 12.0 if -1.0*np.log10(i[2]) > 1.3 else 2.0 for i in dict_c.values()]
sizes_wsn = [24.0 if -1.0*np.log10(v[2]) > 2.0 else 12.0 if -1.0*np.log10(v[2]) >1.3 else 2.0 for k,v in dict_c.iteritems() if "WSN" in k]
color = ["firebrick" if (i[1] - i[0] < -0.58 and -1.0*np.log10(i[2]) > 1.3) else "orange" if (i[1]-i[0] > 0.58 and -1.0*np.log10(i[2]) > 1.3) else "k" for i in dict_c.values()]
plt.figure(figsize=(5,5))
print "down"
down = [k for k,v in dict_c.iteritems() if (v[1] - v[0] < -0.58 and -1.0*np.log10(v[2]) > 1.3)]
for i in down:
	print i
print "up"
up = [k for k,v in dict_c.iteritems() if (v[1] - v[0] > 0.58 and -1.0*np.log10(v[2]) > 1.3)]
for i in up:
	print i

down = sum([1 if (i[1] - i[0] < -0.58 and -1.0*np.log10(i[2]) > 1.3) else 0 for i in dict_c.values()])
up = sum([1 if (i[1] - i[0] > 0.58 and -1.0*np.log10(i[2]) > 1.3) else 0 for i in dict_c.values()])

params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)

plt.scatter(0,0,s=24,c="k",label="p < 0.01")
plt.scatter(0,0,s=12,c="k",label="p < 0.05")
plt.scatter(0,0,s=2,c="k",label ="p > 0.05")
plt.scatter(0,0,s=24,c="firebrick",label="$TE^{DOWN}$"+" n="+str(down),alpha=0.6)
plt.scatter(0,0,s=24,c="orange",label="$TE^{UP}$"+" n="+str(up),alpha=0.6)

		

plt.scatter([i[0] for i in dict_c.values()],[i[1] for i in dict_c.values()],s=sizes,c=color,alpha=0.15)
#for k,v in dict_c.iteritems():
	#if (v[1]-v[0] < -1) or  (v[1] - v[0] > 1):
		#plt.annotate(k,(v[0],v[1]))
plt.scatter([v[0] for k,v in dict_c.iteritems() if "WSN" in k],[v[1] for k,v in dict_c.iteritems() if "WSN" in k],s=sizes_wsn,edgecolors="r",facecolors='none',label="WSN") 
lgnd = plt.legend(loc="upper left",scatterpoints=1,fontsize=8)
lgnd.legendHandles[0]._sizes = [24]
lgnd.legendHandles[1]._sizes = [12] 
lgnd.legendHandles[2]._sizes = [2]
lgnd.legendHandles[3]._sizes = [24]

plt.ylabel('$log_{2}FC$ '+"$RPF^{IFIT2^{-/-}}$"+"/ "+"$RPF^{WT}$  infection",fontsize=10)
plt.xlabel('$log_{2}FC$ '+"$RNA^{IFIT2^{-/-}}$"+"/ "+"$RNA^{WT}$  infection",fontsize=10)
x=np.linspace(-10,10,50)
y1=x+0.58
y2=x-0.58
plt.plot(x,y1,c="k",linewidth=1,linestyle="--",alpha=0.5)
plt.plot(x,y2,c="k",linewidth=1,linestyle="--",alpha=0.5)
plt.xlim(-4,4)
#plt.xlim(-2,2)
plt.ylim(-4,4)
#plt.savefig("rnavribo_ixn.svg")
#plt.savefig("rnavribo_ixn.png")
plt.show()

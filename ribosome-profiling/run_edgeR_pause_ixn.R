#Calculate log2FC(KO_pause/wt_pause) accounting for mock differences of RNA-seq or ribo-seq
#Mitchell Ledwith
library("edgeR")

read_in <- read.delim("ribo_counts.txt",row.names="Sample")
read_in = read_in[,c(1,3,4,6,7,9,10,12)]
group.ribo <- factor(c("wtmo","wtixn","D8mo","D8ixn","wtmo","wtixn","D8mo","D8ixn"))

y.ribo <- DGEList(counts=read_in,group=group.ribo)
y.ribo <- calcNormFactors(y.ribo)

ribo.lib.size <- y.ribo$samples$lib.size
ribo.norm.factors <- y.ribo$samples$norm.factors


read_pause <- read.delim("pauses.txt",row.names="Sample")
read_pause = read_pause[,c(1,3,4,6,7,9,10,12)]
group.pause <- factor(c("wtmo","wtixn","D8mo","D8ixn","wtmo","wtixn","D8mo","D8ixn"))
design.pause <- model.matrix(~0+group.pause)
colnames(design.pause) <- levels(group.pause)

read_offset <- read.delim("pause_ribo_offsetsscaled.txt") #takes offsets scaled to ribo-seq fold-change
read_offset <- read_offset[,c(2,3,4,5,6,7,8,9)]
read_offset <- as.matrix(read_offset)

colnames(read_offset) <- c("wtmo","wtixn","D8mo","D8ixn","wtmo","wtixn","D8mo","D8ixn")


y.pause <- DGEList(counts=read_pause,norm.factors=ribo.norm.factors,lib.size = ribo.lib.size,group=group.pause)
plotMDS(y.pause)
keep <- filterByExpr(y.pause)
y.pause <- y.pause[keep,,keep.lib.sizes=TRUE]
read_offset <- read_offset[keep,]

y.pause <- estimateDisp(y.pause,design.pause)
plotBCV(y.pause)
y.pause <- scaleOffset.DGEList(y.pause,read_offset)



fit <- glmQLFit(y.pause,design.pause,prior.count=5)
fit$offset
my.contrasts <- makeContrasts(D8vwtmo = D8mo-wtmo,D8vwtixn = D8ixn-wtixn,wtixnnorm=wtixn-wtmo,D8ixnnorm = D8ixn-D8mo,D8vwtixnnorm = (D8ixn-D8mo)-(wtixn-wtmo),levels=design.pause)


qlf.total <- glmQLFTest(fit,contrast=my.contrasts[,"D8vwtixnnorm"])
qlf.total
tab.total <- topTags(qlf.total,n=Inf)


write.table(tab.total,file="pause_KOvwtixnnorm.txt",append=FALSE,sep="\t")





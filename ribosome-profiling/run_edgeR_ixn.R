#Calculate log2FC(KO/wt) accounting for mock differences of RNA-seq or ribo-seq
#Mitchell Ledwith
library("edgeR")

read_in <- read.delim("counts.txt",row.names="Sample")
read_in <- read_in[,c(1,3,4,6,7,9,10,12)]
group <- factor(c("wtmo","wtixn","D8mo","D8ixn","wtmo","wtixn","D8mo","D8ixn"))
targets <- data.frame(group)

y <- DGEList(counts=read_in,group=group)
y <- calcNormFactors(y)


design <- model.matrix(~0+group)
colnames(design) <- levels(group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- estimateDisp(y,design)

fit <- glmQLFit(y,design)
my.contrasts <- makeContrasts(D8vwtmo = D8mo-wtmo,D8vwtixn = D8ixn-wtixn,wtixnnorm=wtixn-wtmo,D8ixnnorm = D8ixn-D8mo,D8vwtixnnorm = (D8ixn-D8mo)-(wtixn-wtmo),levels=design)

qlf <- glmQLFTest(fit,contrast=my.contrasts[,"D8vwtixnnorm"])
tab <- topTags(qlf,n=Inf)

write.table(tab,file="/home/mitch/Desktop/ribo/counts/DE/D8vswtixnnorm.txt",append=FALSE,sep="\t")


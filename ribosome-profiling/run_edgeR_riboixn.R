#Calculate log2FC(TE_KO/TE_wt) accounting for mock TE differences
#Mitchell Ledwith
library("edgeR")

read_in_rna <- read.delim("rna_counts.txt",row.names="Sample")
read_in_ribo <- read.delim("ribo_counts.txt",row.names="Sample")

read_in_rna = read_in_rna[,c(1,3,4,6,7,9,10,12)]
read_in_ribo = read_in_ribo[,c(1,3,4,6,7,9,10,12)]

group_rna <- factor(c("wtmo_rna","wtixn_rna","D8mo_rna","D8ixn_rna","wtmo_rna","wtixn_rna","D8mo_rna","D8ixn_rna"))
group_ribo <- factor(c("wtmo_ribo","wtixn_ribo","D8mo_ribo","D8ixn_ribo","wtmo_ribo","wtixn_ribo","D8mo_ribo","D8ixn_ribo"))
group <- factor(c("wtmo_rna","wtixn_rna","D8mo_rna","D8ixn_rna","wtmo_rna","wtixn_rna","D8mo_rna","D8ixn_rna","wtmo_ribo","wtixn_ribo","D8mo_ribo","D8ixn_ribo","wtmo_ribo","wtixn_ribo","D8mo_ribo","D8ixn_ribo"))

targets.rna <- data.frame(group_rna)
targets.ribo <- data.frame(group_ribo)
targets.comb <- data.frame(group)

y.rna <- DGEList(counts=read_in_rna,group=group_rna)
y.rna <- calcNormFactors(y.rna)

y.ribo <- DGEList(counts=read_in_ribo,group=group_ribo)

y.ribo <- calcNormFactors(y.ribo)

design.rna <- model.matrix(~0+group_rna)
design.ribo <- model.matrix(~0+group_ribo)
design.comb <- model.matrix(~0+group)
colnames(design.rna) <- levels(group_rna)
colnames(design.ribo) <- levels(group_ribo)
colnames(design.comb) <- levels(group)
colnames(read_in_rna) <- c("wtmo_rna1","wtixn_rna1","D8mo_rna1","D8ixn_rna1","wtmo_rna2","wtixn_rna2","D8mo_rna2","D8ixn_rna2")
colnames(read_in_ribo) <- c("wtmo_ribo1","wtixn_ribo1","D8mo_ribo1","D8ixn_ribo1","wtmo_ribo2","wtixn_ribo2","D8mo_ribo2","D8ixn_ribo2")



y.rna <- estimateDisp(y.rna,design.rna)
y.ribo <- estimateDisp(y.ribo,design.ribo)
plotBCV(y.rna)

master_counts <- cbind(read_in_rna,read_in_ribo)

master_norm <- c(y.rna$samples$norm.factors,y.ribo$samples$norm.factors)

rna.ribo <- DGEList(counts = master_counts, norm.factors = master_norm, group = group)
keep <- filterByExpr(rna.ribo)
rna.ribo <- rna.ribo[keep,,keep.lib.sizes=TRUE]
plotMDS(rna.ribo)
rna.ribo <- estimateDisp(rna.ribo,design.comb)
plotBCV(rna.ribo)
rna.ribo
fit <- glmQLFit(rna.ribo,design.comb)
fit
my.contrasts <- makeContrasts(D8vwtixnnorm = ((D8ixn_ribo-D8ixn_rna)-(D8mo_ribo-D8mo_rna))-((wtixn_ribo-wtixn_rna)-(wtmo_ribo-wtmo_rna)),levels=design.comb)

print("blah")
qlf <- glmQLFTest(fit,contrast=my.contrasts[,"D8vwtixnnorm"])
qlf
tab <- topTags(qlf,n=Inf)

write.table(tab,file="KOvwtixn_norm_TE.txt",append=FALSE,sep="\t")


library(Biobase)
library(GEOquery)
library(data.table)
library(tidyverse)
library(RColorBrewer)
# creating heatmaps
library(gplots)
# for the jackstraw tests
library(jackstraw)
# calcuating q-values and FDRs
library(qvalue)
# combining ggplots
library(sva)
library(broom)
library(patchwork)
library(limma)
library(edge)
library(gridExtra)

#Loading a GDS file with GEOquery
dat <- getGEO('GSE71956', destdir=".", GSEMatrix= TRUE)

eset <- dat[[1]]
edata <- exprs(eset)
pheno <- pData(eset)

class(eset)
mode(eset)

dim(edata)
dim(pheno)

edata[1:5,1:10]

colnames(pheno)

head(pheno)

# head(pheno, 5)
# pheno %>% distinct(`tissue:ch1`)

pheno <- pheno %>% select("title" = title, "age" = `age:ch1`,
                 "cell_type" = `cell type:ch1`, "diagnosis" = `diagonsis:ch1`)
head(pheno)

# remove the rows with missing values
length(edata)
rows.na <- apply(edata ,1,function(x) sum(is.na(x)))
length(edata[rows.na != 0,])

edata.gathered <- as.data.frame(edata) %>% gather("sample", "value")


ggplot(edata.gathered, aes(x=sample, y=value)) +
    geom_boxplot() + 
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    coord_flip()


ggplot(edata.gathered, aes(x=sample, y=value)) +
    geom_violin() + 
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    coord_flip()

hist(edata[,1])

design = model.matrix(~pheno$age)
arrayw <- arrayWeights(edata, design)
barplot(arrayw, xlab="Array", ylab="Weight", col="white", las=2)
abline(h=1, lwd=1, lty=2)
fitw <- lmFit(edata, design, weights=arrayw)
fitw <- eBayes(fitw)

limma_pvals = topTable(fitw,number=dim(edata)[1])$P.Value
limma_pvals_adj = topTable(fitw,number=dim(edata)[1])$adj.P.Val
hist(limma_pvals,col=4)
hist(limma_pvals_adj,col=4)

# remove the first eigenmatrix from the raw data
raw.svd = svd(edata)
ec = edata - (raw.svd$d[1] * raw.svd$u[,1] %*% t(raw.svd$v[,1]))

# take the log transformation
elv = log(ec^2)
elv.svd = svd(elv)

# remove the first eigenmatrix from the log-transformed data
eclv = elv - (elv.svd$d[1] * elv.svd$u[,1] %*% t(elv.svd$v[,1]))
en = sign(ec) * sqrt(exp(eclv))

# center each row
edata.norm = t(scale(t(en), center=TRUE, scale=TRUE))

edata.norm.gathered <- as.data.frame(edata.norm) %>% gather("sample", "value")


ggplot(edata.norm.gathered, aes(x=sample, y=value)) +
    geom_boxplot() + 
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    coord_flip()


ggplot(edata.norm.gathered, aes(x=sample, y=value)) +
    geom_violin() + 
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    coord_flip()

hist(edata.norm[,1])

mod = model.matrix(~as.factor(diagnosis), data=pheno)
num.sv(edata, mod, method="be")
num.sv(edata, mod, method="leek")

mod = model.matrix(~as.factor(diagnosis), data=pheno)
mod0 = model.matrix(~1, data=pheno)
sva_output = sva(edata.norm, mod, mod0, n.sv=10)

sva_batch <- tibble(SV1=sva_output$sv[,1], 
                    SV2=sva_output$sv[,2],
                    age=as.factor(pheno$age),
                    diagnosis=as.factor(pheno$diagnosis),
                    cell_type=as.factor(pheno$cell_type))



ggplot(sva_batch) + geom_point(aes(x=SV1,y=SV2, col=age))


ggplot(sva_batch) + geom_point(aes(x=SV1,y=SV2, col = diagnosis))


ggplot(sva_batch) + geom_point(aes(x=SV1,y=SV2, col=cell_type))


sva_batch <- tibble(sv1 = sva_output$sv[,1], sv2 = sva_output$sv[,2], cell_type = as.factor(pheno$cell_type))
sva_batch <- gather(sva_batch,"sv","value",-cell_type)
ggplot(sva_batch) + geom_violin(aes(x=cell_type,y=value)) + geom_jitter(aes(x=cell_type,y=value,col=cell_type)) + facet_wrap(~sv, ncol = 1)

m = dim(edata)[1]
n = dim(edata)[2]

geneexp.jackstraw = jackstraw_pca(edata.norm, r=10,
                                  s=round(m*.1), B=10)

# the m p-values of association tests between variables
# and their principal components
hist(geneexp.jackstraw$p.value,10)


# the observed F-test statistics
hist(geneexp.jackstraw$obs.stat,10)

qplot(geneexp.jackstraw$p.value, geom="histogram")

qplot(geneexp.jackstraw$obs.stat, geom="histogram", xlim=c(0,100))


qplot(as.vector(geneexp.jackstraw$null.stat), geom="histogram", xlim=c(0,30))


# lets combine 2 histograms, severly limiting the x-axis
obs.hist <- qplot(geneexp.jackstraw$obs.stat, geom="histogram", xlim=c(0,30))
null.hist <- qplot(as.vector(geneexp.jackstraw$null.stat), geom="histogram", xlim=c(0,30))
print(obs.hist / null.hist)

sum(qvalue(geneexp.jackstraw$p.value)$qvalue < .01)


pval.pc1 = jackstraw_pca(edata, r1=1, r=2, s=round(m*.1), B=10)$p.value

pval.pc2 = jackstraw_pca(edata, r1=2, r=2, s=round(m*.1), B=10)$p.value

# lets combine 2 histograms, severly limiting the x-axis
pc1.hist <- qplot(pval.pc1, geom="histogram")
pc2.hist <- qplot(pval.pc2, geom="histogram")
print(pc1.hist + pc2.hist)

q.pc1 <- qvalue(pval.pc1)
q.pc2 <- qvalue(pval.pc2)

q.pc1$pi0
q.pc2$pi0

sum(q.pc1$qvalue < .01)
sum(q.pc2$qvalue < .01)

tout = rowttests(x = edata, fac = as.factor(pheno$cell_type))

ttidy <- gather(tout)
ggplot(ttidy) + geom_histogram(bins = 30,aes(x=value)) + facet_wrap(~ key, scales="free")

q.obj <- qvalue(tout$p.value)
plot(q.obj)


jackstraw.q <- qvalue(geneexp.jackstraw$p.value)
plot(jackstraw.q)

tout$q.value <- q.obj$qvalues
ttidy <- gather(tout[,c("p.value","q.value")])
ggplot(ttidy) + geom_histogram(bins = 30,aes(x=value)) + facet_wrap(~ key, scales="free")

plot1 <- ggplot(tout) + geom_histogram(bins = 30,aes(x=p.value))
plot2 <- qplot(geneexp.jackstraw$p.value, geom="histogram")
grid.arrange(plot1, plot2, nrow=1)


plot1 <- ggplot(tout) + geom_histogram(bins = 30,aes(x=q.value))
plot2 <- qplot(jackstraw.q$qvalues, geom="histogram")
grid.arrange(plot1, plot2, nrow=1)


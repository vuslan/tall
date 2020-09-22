rm(list=ls())
cat("\014")

library(ALL)
data(ALL)

cat("\014")
ist = grep("T", as.character(ALL$BT))
length(ist)
T = ALL[, ist]

cat("\014")
samplesREF = c('01003',
               '01007',
               '02020',
               '04018',
               '09002',
               '10005',
               '11002',
               '12008',
               '15006',
               '16002',
               '16007',
               '17003',
               '18001',
               '19002',
               '19008',
               '19014',
               '19017',
               '20005',
               '24006',
               '26009',
               '28008',
               '28009',
               '31015',
               '37001',
               '43006',
               '43015',
               '44001',
               '49004',
               '56007',
               '64005',
               '65003',
               '83001')

T32 = T[,samplesREF]
T32
pdata = pData(T32)
summary(pdata)
dfdata = as.data.frame(pData(T32))

cat("\014")
setClass = lapply(T32$CR, function(t) (if ( ("REF" %in% t) || ("DEATH IN INDUCTION" %in% t) ) { return(1) } else { return(0) } ))
dfClass = as.data.frame(setClass)
names(dfClass)=sampleNames(T32)
table(unlist(dfClass))
dfdata['class'] = t(dfClass)

cat("\014")
nGene = 12625L
genes <- featureNames(T32)
head(genes)
matT32exprs = t(exprs(T32)[,])

library(iterativeBMA)
genesel <- BssWssFast(matT32exprs, t(dfClass), 2)

geneidx = genesel$ix[1:nGene]
genes <- featureNames(T32)[geneidx]
head(genes)

generatio = genesel$x[1:nGene]
head(generatio)
tail(generatio)

cat("\014")
library(hgu95av2.db)
?hgu95av2ENTREZID
D = T32[genes,]
symbolsD = select(hgu95av2.db,
                  featureNames(D),
                  c("SYMBOL","ENTREZID", "GENENAME"))

cat("\014")
genesize = 17
D17 = D[c(1:genesize),]
featureNames(D17)
symbolsD17 = select(hgu95av2.db,
                    featureNames(D17),
                    c("SYMBOL","ENTREZID", "GENENAME"))

# HIPK3
cat("\014")

featureNames(D17)[1]
select(hgu95av2.db,
       featureNames(D17)[1],
       c("SYMBOL","ENTREZID", "GENENAME"))

D1 = D[1,]
expdata = cbind(as.data.frame(t(round(exprs(D1),digits=2))),as.numeric(dfClass))
genesize = 1
colnames(expdata)[genesize+1] = "class"

cat("\014")
exprs(D17)[1,]
round(median(exprs(D17)[1,]),2)
round(min(exprs(D17)[1,dfdata['class']==0]),2)
round(mean(exprs(D17)[1,dfdata['class']==0]),2)
round(max(exprs(D17)[1,dfdata['class']==0]),2)
round(min(exprs(D17)[1,dfdata['class']==1]),2)
round(mean(exprs(D17)[1,dfdata['class']==1]),2)
round(max(exprs(D17)[1,dfdata['class']==1]),2)

# FABP4
cat("\014")

featureNames(D17)[4]
select(hgu95av2.db,
       featureNames(D17)[4],
       c("SYMBOL","ENTREZID", "GENENAME"))

D1 = D[4,]
expdata = cbind(as.data.frame(t(round(exprs(D1),digits=2))),as.numeric(dfClass))
genesize = 1
colnames(expdata)[genesize+1] = "class"

exprs(D17)[4,]
round(median(exprs(D17)[4,]),2)
round(min(exprs(D17)[4,dfdata['class']==0]),2)
round(mean(exprs(D17)[4,dfdata['class']==0]),2)
round(max(exprs(D17)[4,dfdata['class']==0]),2)
round(min(exprs(D17)[4,dfdata['class']==1]),2)
round(mean(exprs(D17)[4,dfdata['class']==1]),2)
round(max(exprs(D17)[4,dfdata['class']==1]),2)

# CCL20
cat("\014")

featureNames(D17)[7]
select(hgu95av2.db,
       featureNames(D17)[7],
       c("SYMBOL","ENTREZID", "GENENAME"))

D1 = D[7,]
expdata = cbind(as.data.frame(t(round(exprs(D1),digits=2))),as.numeric(dfClass))
genesize = 1
colnames(expdata)[genesize+1] = "class"

exprs(D17)[7,]
round(median(exprs(D17)[7,]),2)
round(min(exprs(D17)[7,dfdata['class']==0]),2)
round(mean(exprs(D17)[7,dfdata['class']==0]),2)
round(max(exprs(D17)[7,dfdata['class']==0]),2)
round(min(exprs(D17)[7,dfdata['class']==1]),2)
round(mean(exprs(D17)[7,dfdata['class']==1]),2)
round(max(exprs(D17)[7,dfdata['class']==1]),2)



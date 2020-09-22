rm(list=ls())
cat("\014")

library(ALL)
data(ALL)

cat("\014")
ist = grep("T", as.character(ALL$BT))
length(ist)
T = ALL[, ist]

cat("\014")
samplesCR = c('01003',
              '01007',
              '04018',
              '09002',
              '10005',
              '11002',
              '15006',
              '16002',
              '16007',
              '19002',
              '19014',
              '19017',
              '20005',
              '24006',
              '26009',
              '28008',
              '28009',
              '37001',
              '43015',
              '44001',
              '49004',
              '56007',
              '65003',
              '83001')

CR = T[,samplesCR]
CR
pdata = pData(CR)
summary(pdata)

cat("\014")
T22 = CR[,-c(2,4)]
pdata = pData(T22)
summary(pdata)
dfdata = as.data.frame(pData(T22))

cat("\014")
setClass = lapply(T22$f.u, function(t) (if ( ("REL" %in% t) ) { return(1) } else { return(0) } ))
dfClass = as.data.frame(setClass)
names(dfClass)=sampleNames(T22)
table(unlist(dfClass))
dfdata['class'] = t(dfClass)

cat("\014")
nGene = 12625L
genes <- featureNames(T22)
head(genes)
matT22exprs = t(exprs(T22)[,])

library(iterativeBMA)
genesel <- BssWssFast(matT22exprs, t(dfClass), 2)

geneidx = genesel$ix[1:nGene]
genes <- featureNames(T22)[geneidx]
head(genes)

generatio = genesel$x[1:nGene]
head(generatio)
tail(generatio)

cat("\014")
library(hgu95av2.db)
?hgu95av2ENTREZID
D = T22[genes,]
symbolsD = select(hgu95av2.db,
                  featureNames(D),
                  c("SYMBOL","ENTREZID", "GENENAME"))

cat("\014")
genesize = 25
D25 = D[c(1:genesize),]
featureNames(D25)
symbolsD25 = select(hgu95av2.db,
                    featureNames(D25),
                    c("SYMBOL","ENTREZID", "GENENAME"))

# IFITM1
cat("\014")

featureNames(D25)[1]
select(hgu95av2.db,
       featureNames(D25)[1],
       c("SYMBOL","ENTREZID", "GENENAME"))

D1 = D[1,]
expdata = cbind(as.data.frame(t(round(exprs(D1),digits=2))),as.numeric(dfClass))
genesize = 1
colnames(expdata)[genesize+1] = "class"

cat("\014")
exprs(D25)[1,]
round(median(exprs(D25)[1,]),2)
round(min(exprs(D25)[1,dfdata['class']==0]),2)
round(mean(exprs(D25)[1,dfdata['class']==0]),2)
round(max(exprs(D25)[1,dfdata['class']==0]),2)
round(min(exprs(D25)[1,dfdata['class']==1]),2)
round(mean(exprs(D25)[1,dfdata['class']==1]),2)
round(max(exprs(D25)[1,dfdata['class']==1]),2)

# AHNAK
cat("\014")

featureNames(D25)[10]
select(hgu95av2.db,
       featureNames(D25)[10],
       c("SYMBOL","ENTREZID", "GENENAME"))

D1 = D[10,]
expdata = cbind(as.data.frame(t(round(exprs(D1),digits=2))),as.numeric(dfClass))
genesize = 1
colnames(expdata)[genesize+1] = "class"

exprs(D25)[10,]
round(median(exprs(D25)[10,]),2)
round(min(exprs(D25)[10,dfdata['class']==0]),2)
round(mean(exprs(D25)[10,dfdata['class']==0]),2)
round(max(exprs(D25)[10,dfdata['class']==0]),2)
round(min(exprs(D25)[10,dfdata['class']==1]),2)
round(mean(exprs(D25)[10,dfdata['class']==1]),2)
round(max(exprs(D25)[10,dfdata['class']==1]),2)

# PFKFB3
cat("\014")

featureNames(D25)[19]
select(hgu95av2.db,
       featureNames(D25)[19],
       c("SYMBOL","ENTREZID", "GENENAME"))

D1 = D[19,]
expdata = cbind(as.data.frame(t(round(exprs(D1),digits=2))),as.numeric(dfClass))
genesize = 1
colnames(expdata)[genesize+1] = "class"

exprs(D25)[19,]
round(median(exprs(D25)[19,]),2)
round(min(exprs(D25)[19,dfdata['class']==0]),2)
round(mean(exprs(D25)[19,dfdata['class']==0]),2)
round(max(exprs(D25)[19,dfdata['class']==0]),2)
round(min(exprs(D25)[19,dfdata['class']==1]),2)
round(mean(exprs(D25)[19,dfdata['class']==1]),2)
round(max(exprs(D25)[19,dfdata['class']==1]),2)



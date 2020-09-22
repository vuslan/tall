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

featureNames(D17) = c("HIPK3","	HNRNPA3","SNRPA","FABP4","CDC25B","FABP5","CCL20","EFHD1","COL1A2","PDLIM1","SLC15A2","SLC18A2","MEF2C","CKS1B","PSD3","RNASEH2A","ARNTL")
featureNames(D17)

D17
exprs(D17)

?heatmap
colormap = function(t) (if ( ("REF" %in% t) || ("DEATH IN INDUCTION" %in% t) ) { return("#FF0000") } else { return("#0000FF") } )
patientColors = unlist(lapply(D17$CR, colormap))
heatmap(exprs(D17), ColSideColors = patientColors, Rowv = NA)

sampleNames(D17)[dfClass==1]
sampleNames(D17)[dfClass==0]


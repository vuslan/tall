rm(list=ls())
cat("\014")
library(affy)

celpath = "C:/e-mtab-5035"
eSetread <- readExpressionSet("C:/e-mtab-5035/esetmatrix.txt",sep="\t")
summary(eSetread)

expr = t(exprs(eSetread))
summary(expr)

genes <- sapply(featureNames(eSetread), function(x) gsub("\"","",x))
genes

cat("\014")
library(hgu133plus2.db)
?hgu133plus2.db
?hgu133plus2ENTREZID
x <- hgu133plus2ENTREZID
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes][1:300])
if(length(xx)>0) {
  xx[1:5]
  xx[[1]]
}
symbolsD = select(hgu133plus2.db,
                  keys=genes[1:54613],
                  columns=c("SYMBOL","ENTREZID", "GENENAME"))

symbolsC = symbolsD[complete.cases(symbolsD),]
                  


efdata = as.data.frame(expr)
# 1...54613 geneset 54614...54675 AFFX
efdata = efdata[,1:54613]
efdata[1:5,1:5]

sdrf <- read.table(file.choose(),sep="\t", header=TRUE)
colnames(sdrf)
efdata['Source.Name'] = sdrf['Source.Name']
efdata['Characteristics.disease.'] = sdrf['Characteristics.disease.']
efdata['Characteristics.disease.']
# 1...54613 geneset 54614...54615 Characteristics
efdata[1:5,54614:54615]

setClass = lapply(efdata$Characteristics.disease., function(t) (if ( ("T-cell acute lymphoblastic leukemia" %in% t) ) { return(1) } else { return(0) } ))
efClass = as.data.frame(setClass)
table(unlist(efClass))
efdata['class'] = t(efClass)
efdata[1:5,54614:54616]

cat("\014")
nGene = 54613L
head(genes)
library(iterativeBMA)
genesel <- BssWssFast(efdata[,1:54613], t(efClass), 2)

geneidx = genesel$ix[1:nGene]
genes <- genes[geneidx]
head(genes)

generatio = genesel$x[1:nGene]
head(generatio)
tail(generatio)
round(generatio[1:25],2)

t(efdata[1:5,1:5])
D = t(efdata)[genes,]
D[1:5,1:5]
symbolsD = select(hgu133plus2.db,
                  rownames(D)[1:3],
                  c("SYMBOL","ENTREZID", "GENENAME"))

ist = grep("T-cell", as.character(efdata$Characteristics.disease.))
length(ist)
T = efdata[ist,]

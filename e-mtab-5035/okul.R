rm(list=ls())
cat("\014")
library(affy)

celpath = "C:/e-mtab-5035"
celdata = ReadAffy(celfile.path=celpath)
summary(celdata)

eset = rma(celdata)
summary(eset)

expr = exprs(eset)
summary(expr)

pdata = pData(eset)
summary(pdata)

efdata = as.data.frame(expr)
write.table(expr,"C:/e-mtab-5035/esetmatrix.txt",sep="\t")
summary(efdata)

genes <- featureNames(eset)
summary(genes)

expr[5,]

image(celdata[,1])

sdrf <- read.table(file.choose(), sep="\t", header=TRUE)



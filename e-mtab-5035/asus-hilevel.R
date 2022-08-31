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

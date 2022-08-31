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

cat("\014")

nGene = 44700L
head(genes)

library(iterativeBMA)
genesel <- BssWssFast(efdata[,1:nGene], t(efClass), 2)

geneidx = genesel$ix[1:nGene]
genes <- genes[geneidx]
head(genes)

generatio = genesel$x[1:nGene]
head(generatio)
tail(generatio)
round(generatio[1:25],2)

symbolsGENESEL = select(hgu133plus2.db,
                   genes[1:25],
                   c("SYMBOL","ENTREZID", "GENENAME"))

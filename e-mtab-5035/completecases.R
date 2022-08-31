symbolsRaw = select(hgu133plus2.db,
                  keys=genes,
                  columns=c("SYMBOL","ENTREZID", "GENENAME"))

symbolsCC = symbolsRaw[complete.cases(symbolsRaw),]

efdata <- as.data.frame(expr)
efdata <- efdata[rownames(efdata) %in% symbolsCC[,1], ]
efdata[1:5,1:5]

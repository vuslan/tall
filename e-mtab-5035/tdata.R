library(data.table)

cat("\014")
efdata[1:5,1:5]

#transpose data frame
tfdata <- transpose(efdata)

#redefine row and column names
rownames(tfdata) <- colnames(efdata)
colnames(tfdata) <- rownames(efdata)

efdata <- tfdata

efdata[1:5,1:5]

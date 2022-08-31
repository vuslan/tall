sdrf <- read.table(file.choose(),sep="\t", header=TRUE)
colnames(sdrf)
efdata['Source.Name'] = sdrf['Source.Name']
efdata['Characteristics.disease.'] = sdrf['Characteristics.disease.']
# 1...44700 geneset 44701...44702 Characteristics
cat("\014")
efdata[1:5,44701:44702]

ist = grep("T-cell", as.character(efdata$Characteristics.disease.))
length(ist)
T = efdata[ist,]
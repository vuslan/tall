setClass = lapply(efdata$Characteristics.disease., function(t) (if ( ("T-cell acute lymphoblastic leukemia" %in% t) ) { return(1) } else { return(0) } ))
efClass = as.data.frame(setClass)
table(unlist(efClass))
efdata['class'] = t(efClass)
# 1...44700 geneset 44701...44703 Characteristics+Class
cat("\014")
efdata[1:5,44701:44703]

rm(list=ls())
cat("\014")

library(ALL)
data(ALL)

cat("\014")
ist = grep("T", as.character(ALL$BT))
length(ist)
T = ALL[, ist]

cat("\014")
samplesCR = c('01003',
              '01007',
              '04018',
              '09002',
              '10005',
              '11002',
              '15006',
              '16002',
              '16007',
              '19002',
              '19014',
              '19017',
              '20005',
              '24006',
              '26009',
              '28008',
              '28009',
              '37001',
              '43015',
              '44001',
              '49004',
              '56007',
              '65003',
              '83001')

CR = T[,samplesCR]
CR
pdata = pData(CR)
summary(pdata)

cat("\014")
T22 = CR[,-c(2,4)]
pdata = pData(T22)
summary(pdata)
dfdata = as.data.frame(pData(T22))

cat("\014")
setClass = lapply(T22$f.u, function(t) (if ( ("REL" %in% t) ) { return(1) } else { return(0) } ))
dfClass = as.data.frame(setClass)
names(dfClass)=sampleNames(T22)
table(unlist(dfClass))
dfdata['class'] = t(dfClass)

D = T22[c("676_g_at","572_at","37173_at","37171_at","39872_at","40726_at","38158_at","38834_at","41692_at","37027_at","39621_at","527_at","35614_at","37231_at","40738_at","31801_at","37577_at","37333_at","37111_g_at","41050_at","39785_at","38344_at","32767_at","32617_at","32263_at"),]
featureNames(D) = c("IFITM1","TTK","CENPE","KIF23","GTSE1","KIF11","ESPL1","TOPBP1","SYNJ1","AHNAK","OTUD3","CENPA","TCFL5","DLGAP5","CD2","TMF1","ARHGAP19","DNMT1","PFKFB3","TAF5","CEP57","ALMS1","STIL","FANCI","CCNB2")

D
exprs(D)

?heatmap
colormap = function(t) (if ( ("REL" %in% t) ) { return("#FF0000") } else { return("#0000FF") } )
patientColors = unlist(lapply(D$f.u, colormap))
heatmap(exprs(D), ColSideColors = patientColors, Rowv = NA)
sampleNames(D)[dfClass==1]
sampleNames(D)[dfClass==0]

pData(D)["83001",]

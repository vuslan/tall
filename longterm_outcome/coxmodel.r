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

cat("\014")
nGene = 12625L
genes <- featureNames(T22)
head(genes)
matT22exprs = t(exprs(T22)[,])

library(iterativeBMA)
genesel <- BssWssFast(matT22exprs, t(dfClass), 2)

geneidx = genesel$ix[1:nGene]
genes <- featureNames(T22)[geneidx]
head(genes)

generatio = genesel$x[1:nGene]
head(generatio)
tail(generatio)
round(generatio[1:25],2)

cat("\014")
library(hgu95av2.db)
?hgu95av2ENTREZID
D = T22[genes,]
symbolsD = select(hgu95av2.db,
                  featureNames(D),
                  c("SYMBOL","ENTREZID", "GENENAME"))

cat("\014")
genesize = 25
D25 = D[c(1:genesize),]
featureNames(D25)
symbolsD25 = select(hgu95av2.db,
                    featureNames(D25),
                    c("SYMBOL","ENTREZID", "GENENAME"))

cat("\014")
genesize = 3
D3 = D25[c(1,10,19),]
featureNames(D3)
symbolsD3 = select(hgu95av2.db,
                   featureNames(D3),
                   c("SYMBOL","ENTREZID", "GENENAME"))

expdata = cbind(as.data.frame(t(round(exprs(D3),digits=2))),as.numeric(dfClass))
colnames(expdata)[genesize+1] = "class"

cat("\014")
t21pdata <- pData(T22)[-c(12),] # 24006
t21pdata
library(lubridate)
datelastseen = mdy(t21pdata$`date last seen`)
datecr = mdy(t21pdata$date.cr)
datediff = datelastseen - datecr
datediff
imputedate = mean(datediff)
imputedate
imputecr = mdy(pData(T22)[12,21])-imputedate
imputecr
months(imputecr)
pData(T22)[12,8] =
  as.character.Date(imputecr, "%m/%d/%Y")
pData(T22)[12,]
datelastseen =
  mdy(pData(T22)$`date last seen`)
datecr = mdy(pData(T22)$date.cr)
datediff = as.numeric(datelastseen - datecr)
datediff
dfdata['fu_time'] = datediff

cat("\014")
idx=1
gene1 = round(exprs(D25)[idx,],2)
dfdata['gene1'] = gene1
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
IFITM1 <- unlist(dfdata[,"gene1"])
medIFITM1 = round(median(IFITM1),2)
overIFITM1 <- ifelse(gene1>=medIFITM1,1,0) # dichotomise ifitm1
overIFITM1
table(overIFITM1, exclude = NULL) # inspect the numbers
table(IFITM1,overIFITM1, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ IFITM1, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overIFITM1, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "IFITM1 over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overIFITM1)
summary(cox)

cat("\014")
idx=10
gene10 = round(exprs(D25)[idx,],2)
dfdata['gene10'] = gene10
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
AHNAK <- unlist(dfdata[,"gene10"])
medAHNAK = round(median(AHNAK),2)
overAHNAK <- ifelse(gene10>=medAHNAK,1,0) # dichotomise ahnak
overAHNAK
table(medAHNAK, exclude = NULL) # inspect the numbers
table(AHNAK, medAHNAK, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ AHNAK, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overAHNAK, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "AHNAK over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overAHNAK)
summary(cox)

cat("\014")
idx=19
gene19 = round(exprs(D25)[idx,],2)
dfdata['gene19'] = gene19
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
PFKFB3 <- unlist(dfdata[,"gene19"])
medPFKFB3 = round(median(PFKFB3),2)
overPFKFB3 <- ifelse(gene19>=medPFKFB3,1,0) # dichotomise pfkfb3
overPFKFB3
table(medPFKFB3, exclude = NULL) # inspect the numbers
table(PFKFB3, medPFKFB3, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ PFKFB3, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overPFKFB3, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "PFKFB3 over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overPFKFB3)
summary(cox)

cox <- coxph(Surv(fu_time, relapse) ~ overIFITM1 + overAHNAK + overPFKFB3)
summary(cox)

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
idx=2
gene2 = round(exprs(D25)[idx,],2)
dfdata['gene2'] = gene2
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
TTK <- unlist(dfdata[,"gene2"])
medTTK = round(median(TTK),2)
overTTK <- ifelse(gene2>=medTTK,1,0) # dichotomise ttk
overTTK
table(overTTK, exclude = NULL) # inspect the numbers
table(TTK,overTTK, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ TTK, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overTTK, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "TTK over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overTTK)
summary(cox)

cat("\014")
idx=3
gene3 = round(exprs(D25)[idx,],2)
dfdata['gene3'] = gene3
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
CENPE <- unlist(dfdata[,"gene3"])
medCENPE = round(median(CENPE),2)
overCENPE <- ifelse(gene3>=medCENPE,1,0) # dichotomise cenpe
overCENPE
table(overCENPE, exclude = NULL) # inspect the numbers
table(CENPE,overCENPE, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ CENPE, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overCENPE, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "CENPE over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overCENPE)
summary(cox)

cat("\014")
idx=4
gene4 = round(exprs(D25)[idx,],2)
dfdata['gene4'] = gene4
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
KIF23 <- unlist(dfdata[,"gene4"])
medKIF23 = round(median(KIF23),2)
overKIF23 <- ifelse(gene4>=medKIF23,1,0) # dichotomise kif23
overKIF23
table(overKIF23, exclude = NULL) # inspect the numbers
table(KIF23,overKIF23, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ KIF23, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overKIF23, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "KIF23 over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overKIF23)
summary(cox)

cat("\014")
idx=5
gene5 = round(exprs(D25)[idx,],2)
dfdata['gene5'] = gene5
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
GTSE1 <- unlist(dfdata[,"gene5"])
medGTSE1 = round(median(GTSE1),2)
overGTSE1 <- ifelse(gene5>=medGTSE1,1,0) # dichotomise gtse1
overGTSE1
table(overGTSE1, exclude = NULL) # inspect the numbers
table(GTSE1,overGTSE1, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ GTSE1, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overGTSE1, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "GTSE1 over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overGTSE1)
summary(cox)

cat("\014")
idx=6
gene6 = round(exprs(D25)[idx,],2)
dfdata['gene6'] = gene6
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
KIF11 <- unlist(dfdata[,"gene6"])
medKIF11 = round(median(KIF11),2)
overKIF11 <- ifelse(gene6>=medKIF11,1,0) # dichotomise kif11
overKIF11
table(overKIF11, exclude = NULL) # inspect the numbers
table(KIF11,overKIF11, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ KIF11, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overKIF11, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "KIF11 over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overKIF11)
summary(cox)

cat("\014")
idx=7
gene7 = round(exprs(D25)[idx,],2)
dfdata['gene7'] = gene7
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
ESPL1 <- unlist(dfdata[,"gene7"])
medESPL1 = round(median(ESPL1),2)
overESPL1 <- ifelse(gene7>=medESPL1,1,0) # dichotomise espl1
overESPL1
table(overESPL1, exclude = NULL) # inspect the numbers
table(ESPL1,overESPL1, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ ESPL1, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overESPL1, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "ESPL1 over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overESPL1)
summary(cox)

cat("\014")
idx=8
gene8 = round(exprs(D25)[idx,],2)
dfdata['gene8'] = gene8
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
TOPBP1 <- unlist(dfdata[,"gene8"])
medTOPBP1 = round(median(ESPL1),2)
overTOPBP1 <- ifelse(gene8>=medTOPBP1,1,0) # dichotomise topbp1
overTOPBP1
table(overTOPBP1, exclude = NULL) # inspect the numbers
table(TOPBP1,overTOPBP1, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ TOPBP1, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overTOPBP1, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "TOPBP1 over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overTOPBP1)
summary(cox)

cat("\014")
idx=9
gene9 = round(exprs(D25)[idx,],2)
dfdata['gene9'] = gene9
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
SYNJ1 <- unlist(dfdata[,"gene9"])
medSYNJ1 = round(median(SYNJ1),2)
overSYNJ1 <- ifelse(gene9>=medSYNJ1,1,0) # dichotomise synj1
overSYNJ1
table(overSYNJ1, exclude = NULL) # inspect the numbers
table(SYNJ1,overSYNJ1, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ SYNJ1, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overSYNJ1, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "SYNJ1 over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overSYNJ1)
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
idx=11
gene11 = round(exprs(D25)[idx,],2)
dfdata['gene11'] = gene11
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
OTUD3 <- unlist(dfdata[,"gene11"])
medOTUD3 = round(median(OTUD3),2)
overOTUD3 <- ifelse(gene11>=medOTUD3,1,0) # dichotomise otud3
overOTUD3
table(medOTUD3, exclude = NULL) # inspect the numbers
table(OTUD3, medOTUD3, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ OTUD3, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overOTUD3, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "OTUD3 over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overOTUD3)
summary(cox)

cat("\014")
idx=12
gene12 = round(exprs(D25)[idx,],2)
dfdata['gene12'] = gene12
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
CENPA <- unlist(dfdata[,"gene12"])
medCENPA = round(median(CENPA),2)
overCENPA <- ifelse(gene12>=medCENPA,1,0) # dichotomise cenpa
overCENPA
table(medCENPA, exclude = NULL) # inspect the numbers
table(CENPA, medCENPA, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ CENPA, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overCENPA, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "CENPA over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overCENPA)
summary(cox)

cat("\014")
idx=13
gene13 = round(exprs(D25)[idx,],2)
dfdata['gene13'] = gene13
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
TCFL5 <- unlist(dfdata[,"gene13"])
medTCFL5 = round(median(TCFL5),2)
overTCFL5 <- ifelse(gene13>=medTCFL5,1,0) # dichotomise tcfl5
overTCFL5
table(medTCFL5, exclude = NULL) # inspect the numbers
table(TCFL5, medTCFL5, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ TCFL5, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overTCFL5, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "TCFL5 over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overTCFL5)
summary(cox)

cat("\014")
idx=14
gene14 = round(exprs(D25)[idx,],2)
dfdata['gene14'] = gene14
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
DLGAP5 <- unlist(dfdata[,"gene14"])
medDLGAP5 = round(median(DLGAP5),2)
overDLGAP5 <- ifelse(gene14>=medDLGAP5,1,0) # dichotomise dlgap5
overDLGAP5
table(medDLGAP5, exclude = NULL) # inspect the numbers
table(DLGAP5, medDLGAP5, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ DLGAP5, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overDLGAP5, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "DLGAP5 over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overDLGAP5)
summary(cox)

cat("\014")
idx=15
gene15 = round(exprs(D25)[idx,],2)
dfdata['gene15'] = gene15
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
CD2 <- unlist(dfdata[,"gene15"])
medCD2 = round(median(CD2),2)
overCD2 <- ifelse(gene15>=medCD2,1,0) # dichotomise cd2
overCD2
table(medCD2, exclude = NULL) # inspect the numbers
table(CD2, medCD2, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ CD2, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overCD2, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "CD2 over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overCD2)
summary(cox)

cat("\014")
idx=16
gene16 = round(exprs(D25)[idx,],2)
dfdata['gene16'] = gene16
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
TMF1 <- unlist(dfdata[,"gene16"])
medTMF1 = round(median(TMF1),2)
overTMF1 <- ifelse(gene16>=medTMF1,1,0) # dichotomise tmf1
overTMF1
table(medTMF1, exclude = NULL) # inspect the numbers
table(TMF1, medTMF1, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ TMF1, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overTMF1, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "TMF1 over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overTMF1)
summary(cox)

cat("\014")
idx=17
gene17 = round(exprs(D25)[idx,],2)
dfdata['gene17'] = gene17
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
ARHGAP19 <- unlist(dfdata[,"gene17"])
medARHGAP19 = round(median(ARHGAP19),2)
overARHGAP19 <- ifelse(gene17>=medARHGAP19,1,0) # dichotomise arhgap19
overARHGAP19
table(medARHGAP19, exclude = NULL) # inspect the numbers
table(ARHGAP19, medARHGAP19, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ ARHGAP19, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overARHGAP19, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "ARHGAP19 over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overARHGAP19)
summary(cox)

cat("\014")
idx=18
gene18 = round(exprs(D25)[idx,],2)
dfdata['gene18'] = gene18
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
DNMT1 <- unlist(dfdata[,"gene18"])
medDNMT1 = round(median(DNMT1),2)
overDNMT1 <- ifelse(gene18>=medDNMT1,1,0) # dichotomise dnmt1
overDNMT1
table(medDNMT1, exclude = NULL) # inspect the numbers
table(DNMT1, medDNMT1, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ DNMT1, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overDNMT1, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "DNMT1 over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overDNMT1)
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

cat("\014")
idx=20
gene20 = round(exprs(D25)[idx,],2)
dfdata['gene20'] = gene20
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
TAF5 <- unlist(dfdata[,"gene20"])
medTAF5 = round(median(TAF5),2)
overTAF5 <- ifelse(gene20>=medTAF5,1,0) # dichotomise taf5
overTAF5
table(medTAF5, exclude = NULL) # inspect the numbers
table(TAF5, medTAF5, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ TAF5, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overTAF5, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "TAF5 over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overTAF5)
summary(cox)

cat("\014")
idx=21
gene21 = round(exprs(D25)[idx,],2)
dfdata['gene21'] = gene21
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
CEP57 <- unlist(dfdata[,"gene21"])
medCEP57 = round(median(CEP57),2)
overCEP57 <- ifelse(gene21>=medCEP57,1,0) # dichotomise cep57
overCEP57
table(medCEP57, exclude = NULL) # inspect the numbers
table(CEP57, medCEP57, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ CEP57, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overCEP57, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "CEP57 over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overCEP57)
summary(cox)

cat("\014")
idx=22
gene22 = round(exprs(D25)[idx,],2)
dfdata['gene22'] = gene22
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
ALMS1 <- unlist(dfdata[,"gene22"])
medALMS1 = round(median(ALMS1),2)
overALMS1 <- ifelse(gene22>=medALMS1,1,0) # dichotomise alms1
overALMS1
table(medALMS1, exclude = NULL) # inspect the numbers
table(ALMS1, medALMS1, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ ALMS1, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overALMS1, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "ALMS1 over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overALMS1)
summary(cox)

cat("\014")
idx=23
gene23 = round(exprs(D25)[idx,],2)
dfdata['gene23'] = gene23
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
STIL <- unlist(dfdata[,"gene23"])
medSTIL = round(median(STIL),2)
overSTIL <- ifelse(gene23>=medSTIL,1,0) # dichotomise stil
overSTIL
table(medSTIL, exclude = NULL) # inspect the numbers
table(STIL, medSTIL, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ STIL, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overSTIL, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "STIL over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overSTIL)
summary(cox)

cat("\014")
idx=24
gene24 = round(exprs(D25)[idx,],2)
dfdata['gene24'] = gene24
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
FANCI <- unlist(dfdata[,"gene24"])
medFANCI = round(median(FANCI),2)
overFANCI <- ifelse(gene24>=medFANCI,1,0) # dichotomise fanci
overFANCI
table(medFANCI, exclude = NULL) # inspect the numbers
table(FANCI, medFANCI, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ FANCI, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overFANCI, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "FANCI over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overFANCI)
summary(cox)

cat("\014")
idx=25
gene25 = round(exprs(D25)[idx,],2)
dfdata['gene25'] = gene25
View(dfdata)
fu_time <- dfdata[,"fu_time"] # continuous variable (numeric)
relapse <- dfdata[,"class"] # binary variable (numeric)
CCNB2 <- unlist(dfdata[,"gene25"])
medCCNB2 = round(median(CCNB2),2)
overCCNB2 <- ifelse(gene25>=medCCNB2,1,0) # dichotomise ccnb2
overCCNB2
table(medCCNB2, exclude = NULL) # inspect the numbers
table(CCNB2, medCCNB2, exclude = NULL) # check
survdiff(Surv(fu_time, relapse) ~ CCNB2, rho=0)
km_fit <- survfit(Surv(fu_time, relapse) ~ overCCNB2, data = dfdata)
library(survminer)
ggsurvplot(
  fit = km_fit,
  legend.title = "CCNB2 over-expression",
  palette = c("blue", "red"),
  pval = TRUE,
  xlab = "Time (Days)", 
  ylab = "CCR probability")
cox <- coxph(Surv(fu_time, relapse) ~ overCCNB2)
summary(cox)

cox <- coxph(Surv(fu_time, relapse) ~ overIFITM1 + overAHNAK + overPFKFB3)
summary(cox)

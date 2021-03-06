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

library(e1071)
cost = 1.0
epsilon = 0.1
gamma = 0.1
totalacc = 0; totalerr=0; samperr = c()
totalprd = c(); totalact = c();

for(i in 1:22) {
  
  # build svm classifier
  classifier = svm(formula = class ~ ., 
                   data = expdata[-i,], 
                   type = 'C-classification',
                   cost = cost,
                   epsilon = epsilon,
                   gamma = gamma,
                   kernel = 'linear')
  
  svm.predic = predict(classifier, newdata = expdata[i,])   
  predic = svm.predic
  actual = factor(dfClass[1,i],levels=c(0,1))
  
  cfmat = table(predic, actual)
  cfmat
  error <- sum(cfmat) - sum(diag(cfmat))
  error
  accuracy <- round(100- (error * 100 / length(actual)))
  accuracy
  
  totalacc = totalacc + accuracy
  totalerr = totalerr + error*100
  totalprd = c(totalprd, predic)
  totalact = c(totalact, actual)
  print(paste("(",i,")","accuracy = ", as.character(accuracy), "%"), quote=FALSE)
  
  if(accuracy == 0)
  {
    samperr = c(samperr, sampleNames(D3)[i])
  }
  
  rm(svm.predic)
  rm(predic)
  rm(actual)
  rm(cfmat)
  rm(error)
  rm(accuracy)
  
}

loocv.acc = round(totalacc/22,digits=2); loocv.acc
loocv.err = round(totalerr/22,digits=2); loocv.err
samperr
totalcf = table(totalact, totalprd)
totalcf

library(pROC)
exproc <- roc(as.numeric(totalact),
              as.numeric(totalprd),
              smoothed = TRUE,
              # arguments for ci
              ci=TRUE, ci.alpha=0.9, stratified=FALSE,
              # arguments for plot
              plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
              print.auc=TRUE, show.thres=TRUE)

library(precrec)
df <- data.frame(totalact,totalprd)
colnames(df) <- c("totalact", "totalprd")
precrec_obj <- evalmod(scores = df$totalprd, labels = df$totalact)
ggplot2::autoplot(precrec_obj)

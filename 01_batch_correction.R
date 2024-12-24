list.files()

load("matrixa.rda")
dfa<-ok
rm(ok)
load("matrixb.rda")
dfb<-ok
rm(ok)

dim(dfb)

all<-merge(dfa,dfb,by="row.names")

row.names(all)<-all$Row.names
all$Row.names<-NULL

load("hb.rda")
df<-merge(all,edata,by="row.names")
row.names(df)<-df$Row.names
df$Row.names<-NULL
## annotation cross
pheno<-read.csv("phenocross.csv",h=T,row.names=1)
all<-df[,row.names(pheno)]
all(colnames(all)==row.names(pheno))


## combat covariates
library(sva)
batch = pheno$batch
mod = model.matrix(~1+as.factor(group), data=pheno)
+as.factor(group)
combat_edata = ComBat(dat=all,  batch=batch, mod=mod,par.prior=TRUE, prior.plots=TRUE)



edata<-combat_edata

library(preprocessCore)

#edata = edata[rowMeans(edata) > 3, ]

colramp = colorRampPalette(c(3,"white",2))(50)
plot(density(edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.3))
for(i in 1:196){lines(density(edata[,i]),lwd=3,col=colramp[i])}


norm_edata = normalize.quantiles(as.matrix(edata))
plot(density(norm_edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.3))
for(i in 1:196){lines(density(norm_edata[,i]),lwd=3,col=colramp[i])}

colnames(norm_edata)<-colnames(combat_edata)
rownames(norm_edata)<-rownames(combat_edata)

data<-as.data.frame(norm_edata)

save(data,file="normalized.rda")

library(transpipe15)
pcatrans(data,pheno,group="tissue",pal="Set1",alpha=1,names=F)
library(dplyr)
pheno%>%mutate(tissue=case_when(tissue=="Normal"~"Normal cerebellum"))->pheno

##filtration
pheno%>%filter(tissue!="Pilocytic")->pheno
pheno%>%filter(tissue!="Hepatoblastoma cell line")->pheno

data<-data[,row.names(pheno)]

pcatrans(data,pheno,group="group",pal="Set1",alpha=1,names=F)
save(data,file="medulloblastoma.rda")
save(pheno,file="medullopheno.rda")
load(file="medulloblastoma.rda")
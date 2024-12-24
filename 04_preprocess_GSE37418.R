library(GEOquery)

gse <- getGEO("GSE37418", GSEMatrix = TRUE)

# dataset access
set <- gse[[1]]

# expression data and features annotations
data<-exprs(set)
features<-fData(set)
# phenotypes
annot <- pData(set)
tail(features,n=100)

library(dplyr)
library(janitor)
features%>%clean_names()->features
colnames(features)
features%>%select(id,gene_symbol)->features

library(stringr)
df<-as.data.frame(str_split_fixed(features$gene_symbol," /// ",2))
features$gene_symbol<-df$V1
features[features$gene_symbol=='',]<-NA



all<-merge(features,data,by="row.names")

all$Row.names<-NULL


write.csv(annot,file="annotation.csv",row.names=T)
all$id<-NULL
all%>%dplyr::rename(gene="gene_symbol")->all
all<-all[complete.cases(all$gene),]


library(transpipe15)
ok<-filtermatrix(all)
save(ok,file="matrix.rda")

save(features,file="featuresGPL570.rda")

pheno<-read.csv("annotation.csv",h=T,row.names=1)

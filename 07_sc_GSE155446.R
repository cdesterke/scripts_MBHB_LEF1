list.files()





library(data.table)
data<-fread("GSE155446_human_raw_counts.csv")

head(data[,1:5])

data<-as.data.frame(data)
row.names(data)<-data$gene
data$gene<-NULL

pheno<-read.csv("GSE155446_human_cell_metadata.csv",h=T,row.names=1)

all(colnames(data)==row.names(pheno))
library(Seurat)

mb  <- CreateSeuratObject(counts = data, min.cells = 3, project = "mb")
mb <- AddMetaData(object = mb,  metadata = pheno)

Idents(object = mb) <- 'geo_sample_id'
Idents(object = mb)

library(dplyr)

mb%>%NormalizeData()%>%FindVariableFeatures()%>%ScaleData()->mb

mb%>%RunPCA()->mb

ElbowPlot(mb,ndims=50)

mb<-RunUMAP(mb,dims=1:15)



DimPlot(mb,reduction="umap")
DimPlot(mb,reduction="umap",group.by="tumor_subpopulation")

wnt <- list(c('LEF1','AXIN2','FGFR2','DKK1','NXN','TCF7L1','NFATC4','NKD1','RNF43','DKK4','LGR6','STK3','YAP1'))
mb<- AddModuleScore(
  object = mb,
  features = wnt, name = 'wnt.lef',
ctrl=13
)

library(ggplot2)
colgex=c("grey90",brewer.pal(7,"Reds"))
plotTheme=theme_classic(base_size=16)
FeaturePlot(mb,reduction="umap",pt.size=0.1,features="wnt.lef1",min.cutoff = "q9",col=c("grey90","darkred"),split.by="coarse_cell_type")+
	plotTheme+coord_fixed()


x<-mb[[]]
head(x)

save(mb,file="mb_seurat.rda")

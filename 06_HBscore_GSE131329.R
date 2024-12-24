



##########
list.files()
pheno<-read.table("pheno_tumor.tsv",h=T,row.names=1)


wnt<-c("LEF1","AXIN2","FGFR2","DKK1","NXN","TCF7L1","NFATC4","NKD1","RNF43","DKK4","LGR6","STK3","YAP1")

data<-read.table("matrix.csv",sep=",",h=T)

library(transpipe15)
data<-data[complete.cases(data$gene),]

ok<-filtermatrix(data)

annot<-read.csv("input_data_meta.csv",h=T,row.names=1)

all(colnames(ok)==row.names(annot))

process<-ok[row.names(ok)%in%wnt,]


## perform PCA without ellipses
pcatrans(process,annot,group="tissue",pal="Set1",alpha=1,names=F)


library(dplyr)
annot%>%select(tissue)->annotation
## draw heatmap
bestheat(process,annotation,font=10,rownames=T,scale="none")

process2<-process[,row.names(pheno)]
## perform PCA without ellipses
pcatrans(process2,pheno,group="pretext_stage",pal="Set1",alpha=1,names=F)
pcatransellipses(process2,pheno,group="pretext_stage",pal="Set1",alpha=0.5,names=F,x=1,y=2,level=0.75)
pheno%>%select(cairo,pretext_stage,clinical_course)->phenotype
bestheat(process2,phenotype,font=10,rownames=T,scale="none")



best<-function (process, pheno, scale = "row", font = 6, rownames = TRUE) 
{
    if (!require(dplyr)) {
        install.packages("dplyr")
    }
    library(dplyr)
    if (!require(pheatmap)) {
        install.packages("pheatmap")
    }
    library(pheatmap)
    if (rownames == TRUE) {
        pheatmap(process, scale = scale, color = colorRampPalette(c("navy", 
            "white", "firebrick3"))(50), annotation = pheno, 
            fontsize = font, cutree_rows = 1, cutree_col = 1, 
            clustering_method = "ward.D2", clustering_distance_cols = "correlation", 
            clustering_distance_rows = "euclidean", show_colnames = F, 
            show_rownames = T, annotation_names_col = T, annotation_names_row = T)
    }
    else {
        pheatmap(process, scale = scale, color = colorRampPalette(c("navy", 
            "white", "firebrick3"))(50), annotation = pheno, 
            fontsize = font, cutree_rows = 1, cutree_col = 1, 
            clustering_method = "ward.D2", clustering_distance_cols = "correlation", 
            clustering_distance_rows = "euclidean", show_colnames = F, 
            show_rownames = F, annotation_names_col = T, annotation_names_row = T)
    }
}

bestheat(process2,phenotype,font=10,rownames=T,scale="none")


pheno<-read.csv("input_data_meta.csv",h=T,row.names=1)
library(Publish)
colnames(pheno)

u<-univariateTable(tissue~gender+age_months+pretext_stage+chic_risk_stratification+
ctnnbi_gene_alteration.ch1+clinical_course+histological_type,data=pheno)
res<-summary(u,show.missing =  "never")
res
write.table(res,file="gse131329baseline.tsv",row.names=F,sep="\t")

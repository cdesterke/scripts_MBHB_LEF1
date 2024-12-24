load("datamb.rda")
load("phenomb.rda")

data<-log(sel+1,2)


lef1<-read.table("promoter.tsv",h=T,sep="\t")


library(dplyr)

lef<-lef1[complete.cases(lef1$gene_symbol),]

library(dplyr)

lef%>%select(gene_symbol)%>%distinct()%>%pull(gene_symbol)->vector

small<-data[row.names(data)%in%vector,]

save(small,file="data_lef1.rda")

table(annot$MOLECULAR_SUBTYPE,annot$CANCER_TYPE_DETAILED)

library(Hmisc)

describe(annot)

pheno<-annot[complete.cases(annot$MOLECULAR_SUBTYPE),]



small<-small[,row.names(pheno)]

trans<-as.data.frame(t(small))

pheno%>%mutate(group=case_when(MOLECULAR_SUBTYPE=="MB, Group3"~"G3",MOLECULAR_SUBTYPE=="MB, Group4"~"G4",
						MOLECULAR_SUBTYPE=="MB, SHH"~"SHH",MOLECULAR_SUBTYPE=="MB, WNT"~"WNT"))->pheno

pheno%>%mutate(wnt=ifelse(group=="WNT","YES","no"))->pheno


library(transpipe15)

res<-deg(small,pheno$wnt,control="no")

vollimma(res,nb=500,fc=1,p=0.0161,size=4,alpha=1)

sig<-filtresig(res)

write.table(sig,file="limma_sig.tsv",row.names=T,sep="\t")

sig%>%filter(logFC>1)->pos

dim(pos)

pheno%>%select(group,wnt)->phenotype



process<-reducedf(pos,small,n=141)


## perform PCA without ellipses
pcatrans(process,phenotype,group="group",pal="Set1",alpha=1,names=F)

## draw heatmap
bestheat(process,phenotype,font=8,rownames=F)


pos$gene_symbol<-row.names(pos)

pos%>%left_join(lef,by="gene_symbol")->annotations

write.table(annotations,file="signature_integrate.tsv",row.names=F,sep="\t")

library(ggplot2)
library(ggrepel)

ggplot(annotations, aes(x = logFC, y = TSS_distance)) +
  geom_point(aes(size = macs_score), color = "blue") +  # Points avec taille variable
  geom_text_repel(aes(label = gene_symbol), size = 5, box.padding = 0.5) +  # Labels avec ggrepel
  theme_minimal() +  # Thème minimal


library(tfcombined)

vector<-row.names(pos)

res<-tfcrcalc(vector)
head(res,n=15)

res%>%filter(pvalues<=0.05)->res

tfcrnet(vector,res,layout=layout_nicely,cex=1,distance=1.5)
plotnlp(res)


library(clusterProfiler)
library(org.Hs.eg.db)


ego <- enrichGO(gene= pos$gene_symbol, OrgDb= org.Hs.eg.db,keyType= 'SYMBOL',ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
                
                
  ego <- enrichGO(gene= row.names(pos), OrgDb= org.Mm.eg.db,keyType= 'SYMBOL',ont= "BP",pAdjustMethod = "BH",qvalueCutoff  = 0.05)              
                


goplot(ego)

dotplot(ego, showCategory=10)
library(DOSE)
#BiocManager::install("enrichplot")
library(enrichplot)

genelist<-split(pos$logFC,pos$gene_symbol)

cnetplot(ego,foldChange=genelist,colorEdge=TRUE)
,showCategory=10)

ego2 <- pairwise_termsim(ego)
treeplot(ego2)

upsetplot(ego)


wnt<-c("LEF1","AXIN2","FGFR2","DKK1","NXN","TCF7L1","NFATC4","NKD1","RNF43","DKK4","LGR6","STK3","YAP1")

smallwnt<-process[row.names(process)%in%wnt,]

bestheat(smallwnt,phenotype,font=8,rownames=T)

pcatrans(smallwnt,phenotype,group="group",pal="Set1",alpha=1,names=F)

expression<-as.data.frame(t(smallwnt))
expression$wnt<-phenotype$wnt
library(logitloop)
library(dplyr)

expression$wnt<-as.factor(expression$wnt)
df<-logitloop(expression,outcome="wnt")
df
plotcoef(df,nb=13,title="")
df%>%filter(significance=="YES")->df
all<-cbind(pheno,expression)



id<-df$predictors
beta<-df$beta

equation<-paste(id,beta,sep="*",collapse=")+(")
equation<-paste("(",equation,")",sep="")
equation







expression%>%mutate(wnt_score=(LGR6*2.63142891330484)+(DKK1*0.955176811196931)+(STK3*1.65636190892634)+(RNF43*2.25764420724335)+
(YAP1*0.783227543742029)+(NFATC4*0.704123299746549)+(NXN*1.36092746269893)+(FGFR2*0.561794152596338)+
(TCF7L1*0.84193613738095)+(LEF1*3.21492490198806)+(DKK4*7.0586778644224)+(AXIN2*10.3595989356147)+
(NKD1*26.9323268141981))->expression

all$wnt_score<-expression$wnt_score

formule<-paste(wnt,collapse="+")
library(pROC)
library(Epi)

ROC(form = wnt~LEF1+FGFR2+DKK1+NXN+TCF7L1+NFATC4+RNF43+DKK4+LGR6+STK3+YAP1 , plot="ROC", data=expression)
all<-as.data.frame(all)
all%>%select(-51)%>%head()


library(cutpointr)

all<-all[,-51]
cp <- cutpointr(all, wnt_score,  wnt,method = maximize_metric, metric = sum_sens_spec)

plot(cp)
cp
all%>%mutate(wnt_cat=ifelse(wnt_score>=354.628,"HIGH","low"))->all

chisq.test(table(all$wnt_cat,all$group))



## mosaicplot

library(vcd)

struct <- structable(~ wnt_cat+CNS_REGION,data = all)
mosaic(struct, , direction = "h", pop = FALSE,colorize = T, shade = TRUE)
       #gp = gpar(fill = matrix(c("red","grey90" , "grey90","grey90" , "grey90", "green3"), 2, 3)))
labeling_cells(text = as.table(struct), margin = 0)(as.table(struct))
chisq.test(struct)


library(Publish)
u<-univariateTable(wnt_cat~SEX+TUMOR_TYPE+group+CNS_REGION+EFS_STATUS,data=all)
res<-summary(u)
res

write.table(res,file="univariate.tsv",row.names=F,sep="\t")

save(all,file="wntscore.rda")

table(all$wnt,all$wnt_cat)
load(file="wntscore.rda")

library(dplyr)
all%>%mutate(wnt_score=(LGR6*2.63142891330484)+(DKK1*0.955176811196931)+(STK3*1.65636190892634)+(RNF43*2.25764420724335)+
(YAP1*0.783227543742029)+(NFATC4*0.704123299746549)+(NXN*1.36092746269893)+(FGFR2*0.561794152596338)+
(TCF7L1*0.84193613738095)+(LEF1*3.21492490198806)+(DKK4*7.0586778644224)+(AXIN2*10.3595989356147)+
(NKD1*26.9323268141981))->all

tf<-read.csv("tf.csv",h=T)


tf%>%select(HGNC.symbol,DBD,Is.TF.)%>%filter(Is.TF.=="Yes")->tf

tf%>%dplyr::rename(gene="HGNC.symbol")->tf

library(transpipe)
res<-deg(data,pheno$group,control="control")
res<-read.table("limmacross.tsv",h=T,sep="\t")
## filter significant DEGs
sig<-filtresig(res)

## draw volcanoplot
vollimma(res,nb=500,fc=1,p=0.01,size=4,alpha=1)

process<-reducedf(sig,data,n=1259)

pcatrans(process,pheno,group="tissue",pal="Set1",alpha=1,names=F)

pheno%>%select(1:2,group)->annot

bestheat(process,annot,font=10,rownames=F)

write.table(res,file="limmacross.tsv",row.names=T,sep="\t")

write.table(sig,file="SIGcross.tsv",row.names=T,sep="\t")

save(pheno,file="phenofinal.rda")
save(data,file="datafinal.rda")

# Chargement des librairies nécessaires
library(ggplot2)
library(dplyr)
library(ggrepel)

# Transformation des p-values en -log10(pvalue)
res$gene <- row.names(res)
res <- res %>%
  mutate(logP = -log10(P.Value),
         is_inside = ifelse(gene %in% vector, "Yes", "No"),
         is_significant = ifelse(abs(logFC) > 1 & adj.P.Val <= 0.05, "Significant", "Not Significant"))

# Volcanoplot avec ggplot2
ggplot(res, aes(x = logFC, y = logP)) +
  # Points pour les gènes non significatifs ou non "inside"
  geom_point(data = res %>% filter(is_inside == "No" & is_significant == "Not Significant"),
             aes(x = logFC, y = logP),
             color = "grey80") +
  # Points pour les gènes "inside" mais non significatifs
  geom_point(data = res %>% filter(is_inside == "Yes" & is_significant == "Not Significant"),
             aes(x = logFC, y = logP),
             color = "orchid", size = 1) +
  # Points pour les gènes significatifs
  geom_point(data = res %>% filter(is_significant == "Significant"),
             aes(x = logFC, y = logP),
             color = "lightblue", size = 1.5) +
  # Ajouter les étiquettes pour les gènes "inside"
  ggrepel::geom_text_repel(
    data = res %>% filter(is_inside == "Yes"),  
    aes(label = gene, color = is_inside),  
    box.padding = 0.2, point.padding = 0.2,size=4.5) +
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(size = 20, face = "bold"), 
    axis.title = element_text(size = 16),  
    axis.text = element_text(size = 14),   
    strip.text = element_text(size = 16),  
    text = element_text(size = 14)) +
  labs(title = "",
       x = "Log2 Fold Change",
       y = "-log10(P-value)",
       color = "geneset")

load("harmony.rda")

x<-data[[]]
head(x)

wnt <- list(c('LEF1','AXIN2','FGFR2','DKK1','NXN','TCF7L1','NFATC4','NKD1','RNF43','DKK4','LGR6','STK3','YAP1'))
data<- AddModuleScore(
  object = data,
  features = wnt, name = 'wnt.lef',
ctrl=13
)

library(ggplot2)
colgex=c("grey90",brewer.pal(7,"Reds"))
plotTheme=theme_classic(base_size=12)
FeaturePlot(data,reduction="umap",pt.size=0.1,features="wnt.lef1",min.cutoff = "q9",col=c("grey90","darkred"),split.by="sample_group")+
	plotTheme+coord_fixed()

Cell.class

library(pals)
DimPlot(data,reduction="umap",group.by ="Cell.class",pt.size=1,cols=cols25(),label=F)

VlnPlot(data, features = c("wnt.lef1"), slot = "data", log = TRUE,pt.size=1,split.by="Cell.class_concise",col=cols25())





# Liste des groupes uniques
groupes <- unique(x$Cell.class)

# Créer un data frame pour stocker les résultats
resultats_df <- data.frame(
  Comparaison = character(),
  Moyenne_Groupe1 = numeric(),
  Moyenne_Groupe2 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Boucle pour faire les tests t pour chaque paire de groupes
for (i in 1:(length(groupes) - 1)) {
  for (j in (i + 1):length(groupes)) {
    # Sélection des deux groupes à comparer
    groupe1 <- x$wnt.lef1[x$Cell.class == groupes[i]]
    groupe2 <- x$wnt.lef1[x$Cell.class == groupes[j]]
    
    # Test t
    test <- t.test(groupe1, groupe2)
    
    # Ajouter les résultats dans le data frame
    resultats_df <- rbind(resultats_df, data.frame(
      Comparaison = paste(groupes[i], "vs", groupes[j]),
      Moyenne_Groupe1 = mean(groupe1),
      Moyenne_Groupe2 = mean(groupe2),
      p_value = test$p.value
    ))
  }
}

# Afficher les résultats sous forme de data frame
print(resultats_df)
write.csv(resultats_df,file="comparisons_cellsgroups_ferroptosisUP.csv",row.names=F)

tumor<-resultats_df[grepl("Tumor",resultats_df$Comparaison),]

write.table(tumor,file="lef1tumor.tsv",row.names=F,sep="\t")

write.table(x,file="lef1meta.tsv",row.names=T,sep="\t")

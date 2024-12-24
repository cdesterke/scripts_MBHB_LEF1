
wnt<-c("LEF1","AXIN2","FGFR2","DKK1","NXN","TCF7L1","NFATC4","NKD1","RNF43","DKK4","LGR6",
"STK3","YAP1")

pheno<-read.csv("annotationb.csv",h=T,row.names=1)

list.files()

load("matrixb.rda")

all(colnames(ok)==row.names(pheno))

small<-ok[row.names(ok)%in%wnt,]

library(transpipe15)

library(dplyr)
pheno%>%select(group,age_months,gender)->phenotype

bestheat(small,phenotype,font=10,rownames=T)

pcatrans(small,pheno,group="group",pal="Set1",alpha=1,names=F)


phenotype<-pheno[complete.cases(pheno$group),]
wnt<-small[,row.names(phenotype)]

phenotype%>%mutate(WNT_STATUS=ifelse(group=="WNT","YES","no"))->phenotype

library(logitloop)

expression<-as.data.frame(t(wnt))




library(glmnet)
y <- as.factor(phenotype$WNT_STATUS)
X<-as.matrix(expression)

library(caret)
trainIndex <- createDataPartition(y, p = 0.6, list = FALSE)
X_train <- X[trainIndex,]
y_train <- y[trainIndex]
X_test <- X[-trainIndex,]
y_test <- y[-trainIndex]


# sequence for tuning of lambda and alpha
lambda_grid <- 10^seq(3, -2, by = -0.1)
alpha_grid <- seq(0.1, 0.9, by = 0.1)

# initialisation of variables
results <- expand.grid(alpha = alpha_grid, lambda = lambda_grid)
results$auc <- NA
library(pROC)
# loop for lambda and alpha tuning
for (i in 1:nrow(results)) {
  alpha <- results$alpha[i]
  lambda <- results$lambda[i]
  
  model <- glmnet(X_train, y_train, alpha = alpha, lambda = lambda, family = "binomial")
  
  # Prédictions on test data
  probs <- predict(model, s = lambda, newx = X_test, type = "response")
  
  # compute AUC
  roc_obj <- roc(y_test, as.vector(probs))
  auc <- auc(roc_obj)
  
  # data results
  results$auc[i] <- auc
}


# Selection of best parameters
best_params <- results[which.max(results$auc),]
best_params
best_alpha <- best_params$alpha

best_lambda <- best_params$lambda
best_auc <- best_params$auc

library(ggplot2)
# tuning visualisation
results <- results %>%
  mutate(log_lambda = log10(lambda))
library(pals)
ggplot(results, aes(x = log_lambda, y = auc, color = as.factor(alpha))) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = log10(best_lambda), linetype = "dashed", color = "red",size=1) +
	scale_color_manual(values=glasbey())+
  labs(title = "AUC depending on lambda and alpha",
       x = "log10(lambda)",
       y = "AUC",
       color = "Alpha") +
  theme_minimal()+theme(text = element_text(size = 16))



### elastic net 0.1 alpha 

fit <- glmnet(X, y, family = "binomial",alpha=0.1)
plot(fit)

cvfit <- cv.glmnet(X, y, family = "binomial",alpha=0.1)
plot(cvfit)
coefficients<-coef(fit, s = cvfit$lambda.min)
coefficients<-as.matrix(coefficients)
coefficients<-as.data.frame(coefficients)

library(dplyr)
coefficients %>% dplyr::rename(coef = "s1") ->selcoef
selcoef$gene<-row.names(selcoef)
selcoef%>%arrange(desc(coef))->selcoef
selcoef%>%filter(gene != "(Intercept)")->selcoef
selcoef
selcoef%>%filter(coef>0)->pos
pos
dim(pos)

write.table(selcoef,file="selcoef.tsv",row.names=T,sep="\t")

library(ggplot2)
ggplot(data=selcoef,aes(x=reorder(gene,coef),y=coef))+geom_bar(stat="identity",fill="plum2")+
coord_flip()+theme_minimal()+xlab("Genesets") + ylab("Normalized enrichment score")+
geom_text(aes(label=round(coef,3)),hjust=0, vjust=0.5,color="darkblue",position= position_dodge(0),size=4,angle=0)+
xlab("Genes") + ylab("Elasticnet coefficients")+
ggtitle("") +theme(text = element_text(size = 14))+
theme(legend.position = "none")



id<-selcoef$gene
beta<-selcoef$coef

equation<-paste(id,beta,sep="*",collapse=")+(")
equation<-paste("(",equation,")",sep="")
equation



 
all<-cbind(phenotype,expression)

all%>%mutate(en_score=(AXIN2*1.2020701063996)+(RNF43*0.789193064889547)+(STK3*0.726759875595984)+(LEF1*0.707115524960767)+
(DKK4*0.675515112192303)+(NKD1*0.588331959358218)+(LGR6*0.489787531073418)+(DKK1*0.450485900900116)+
(NXN*0.290530670569544)+(FGFR2*0.24110284112139)+(YAP1*0.234049384669222)+(TCF7L1*0.0923628089758896)+
(NFATC4*0))->all

all%>%mutate(wnt_score=(LGR6*2.63142891330484)+(DKK1*0.955176811196931)+(STK3*1.65636190892634)+(RNF43*2.25764420724335)+
(YAP1*0.783227543742029)+(NFATC4*0.704123299746549)+(NXN*1.36092746269893)+(FGFR2*0.561794152596338)+
(TCF7L1*0.84193613738095)+(LEF1*3.21492490198806)+(DKK4*7.0586778644224)+(AXIN2*10.3595989356147)+(NKD1*26.9323268141981))->all

library(pROC)
library(Epi)

ROC(form = WNT_STATUS~LEF1+FGFR2+DKK1+NXN+TCF7L1+NFATC4+RNF43+DKK4+LGR6+STK3+YAP1 , plot="ROC", data=all)


library(cutpointr)


cp <- cutpointr(all, en_score,  WNT_STATUS,method = maximize_metric, metric = sum_sens_spec)

plot(cp)
cp
all%>%mutate(wnt_cat=ifelse(wnt_score>=485.982,"HIGH","low"))->all
all%>%mutate(en_cat=ifelse(en_score>=49.8293,"HIGH","low"))->all
chisq.test(table(all$wnt_cat,all$WNT_STATUS))

## mosaicplot

library(vcd)

struct <- structable(~ en_cat+WNT_STATUS,data = all)
mosaic(struct, , direction = "h", pop = FALSE,colorize = T, shade = TRUE)
       #gp = gpar(fill = matrix(c("red","grey90" , "grey90","grey90" , "grey90", "green3"), 2, 3)))
labeling_cells(text = as.table(struct), margin = 0)(as.table(struct))
chisq.test(struct)


library(Publish)
u<-univariateTable(wnt_cat~gender+age_months+anapath+ethnie+group,data=all)
res<-summary(u)
res

write.table(res,file="univariate.tsv",row.names=F,sep="\t")

save(all,file="wntscore.rda")
all$wnt_cat<-as.factor(all$wnt_cat)
all$wnt_cat<-relevel(all$wnt_cat,ref="low")
table(all$wnt,all$wnt_cat)
m<-glm(wnt_cat~gender+ethnie+age_months+anapath,family="binomial",data=all)
summary(m)




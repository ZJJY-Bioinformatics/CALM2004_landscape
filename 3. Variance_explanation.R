library(phyloseq)
library(microbiome)
library(metafor)
library(metap)
library(tidyverse)
library(foreach)
library(crayon)
library(vegan)
library(MASS)
library(lme4)
library(ggsci)
library(ggpubr)
library(ggVennDiagram)
library(car)
library(MatchIt)
library(glmnet)


## use lasso to do feature selection
lasso_fit=function(feature,outcome){
  if(length(colnames(feature))==0){
    result=data.frame(Protein=colnames(outcome),Feature="No",Lasso.beta=0)
  }else if(length(colnames(feature))==1){
    model=lm(outcome[,1]~ feature[,1])
    beta = summary(model)$coefficients[2]
    result=data.frame(Protein=colnames(outcome),Feature=colnames(feature),Lasso.beta=beta)
  }else{
    #x_train=model.matrix(~.-1,feature)
    cv=cv.glmnet(feature,as.matrix(outcome), alpha = 1, nfolds = 10, type.measure="mse",standardize=T)
    beta <- (coef(cv, s = "lambda.min"))
    beta=as.data.frame(beta[-1,1])
    beta$variable=rownames(beta)
    colnames(beta)[1]="beta"
    result=data.frame(Protein=colnames(outcome),Feature=beta$variable,Lasso.beta=beta$beta)
  }
}

### 1. run the feature selection part
res.all<-list()
for(mm in 1:100){
  
  set.seed(mm*233+1001)
  
  trainset= foreach(i=1:length(taxa.test),.combine = rbind) %do%  {
    
    tmp.protein=taxa.test[i]
    tmp.feature=as.character(colnames(factors))
    
    cat(yellow(tmp.protein,"lasso running","\n"))
    
    feature=na.omit(factors[,tmp.feature,drop=F])
    outcome=na.omit(deconvolution_sub_train[,tmp.protein,drop=F])
    
    feature=feature[rownames(feature) %in% rownames(outcome),,drop=F]
    outcome=outcome[rownames(outcome) %in% rownames(feature),,drop=F]
    
    feature=feature[order(rownames(feature)),,drop=F]
    outcome=outcome[order(rownames(outcome)),,drop=F]
    
    result=lasso_fit(feature,outcome)
    return.string = result
  }
  
  trainset=trainset[trainset$Lasso.beta!=0,]
  
  testset.variation=foreach(i=1:length(unique(trainset$Protein)),.combine = rbind) %do%  {
    
    pro=as.character(unique(trainset$Protein)[i])
    cov=as.character(trainset$Feature[trainset$Protein==pro])
    
    tmp.pro=deconvolution_sub_train[,pro,drop=F]
    tmp.cov=factors[,colnames(factors) %in% cov,drop=F]
    tmp.pro=na.omit(tmp.pro)
    tmp.cov=na.omit(tmp.cov)
    tmp.pro=(tmp.pro[rownames(tmp.pro) %in% rownames(tmp.cov),,drop=F])
    tmp.cov=(tmp.cov[rownames(tmp.cov) %in% rownames(tmp.pro),,drop=F])
    tmp.pro=tmp.pro[order(rownames(tmp.pro)),,drop=F]
    tmp.cov=tmp.cov[order(rownames(tmp.cov)),,drop=F]
    
    tmp.glm=glm(tmp.pro[,1]~.,data=as.data.frame(tmp.cov),family = gaussian)
    tmp.coef=as.data.frame(summary(tmp.glm)$coef)
    tmp.coef$FDR=p.adjust(tmp.coef$`Pr(>|t|)`)
    tmp.coef$CIup=tmp.coef$Estimate+1.96*tmp.coef$`Std. Error`
    tmp.coef$CIdown=tmp.coef$Estimate-1.96*tmp.coef$`Std. Error`
    tmp.av=anova(tmp.glm)
    tmp.av$Explain=NA
    for(j in 2:nrow(tmp.av)){
      tmp.av$Explain[j]=(tmp.av$`Resid. Dev`[j-1]-tmp.av$`Resid. Dev`[j])/tmp.av$`Resid. Dev`[1]
    }
    tmp.av=merge(tmp.av,tmp.coef,by="row.names",all=F)
    
    if(nrow(tmp.av)==0){
      return.string=data.frame(MultiVariate=NA,MultiVariate.FDR=NA,MultiVariate.explain=NA,Protein=pro,CIup=NA,CIdonw=NA,Estimate=NA,Run=mm)
    }else{
      return.string=data.frame(MultiVariate=tmp.av$Row.names,MultiVariate.FDR=tmp.av$FDR,
                               MultiVariate.explain=tmp.av$Explain,Protein=pro,CIup=tmp.av$CIup,CIdonw=tmp.av$CIdown,
                               Estimate=tmp.av$Estimate,Run=mm)
    }
  }
  
  res.all[[mm]]=testset.variation
  cat(green(mm,"fold is done","=====================","\n"))
}


### plot results
p1<-ggplot(df.plot, aes(fill=group, x=VE.mean, y=Bacteria)) + 
  geom_bar(position="stack", stat="identity",alpha=0.8)+
  theme_bw()+
  geom_text(aes(label=lab),position=position_stack(vjust=.5),size=2)+
  scale_y_dendrogram(hclust=hc,position = "right")+
  scale_fill_nejm()+ylab("")+xlab("Variance explained (R2)")+
  theme(legend.position = "bottom")


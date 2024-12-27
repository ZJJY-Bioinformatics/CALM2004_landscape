
### ===================================================
### 1. Single variable adonis

result.model1 = foreach(i=1:length(col.test),.combine = rbind) %do%  {
  
  pheno.col=col.test[i]
  tmp<-metadata[,c("NEW_ID",pheno.col)]
  tmp<-na.omit(tmp)
  
  dist.mat1<-dist_subset(dist.mat,tmp$NEW_ID)
  
  ad1 = adonis(dist.mat1 ~ tmp[,pheno.col],permutations = 1000)

  results = data.frame(ad1$aov.tab[1,])
  results$test = pheno.col
  results$N = dim(tmp)[1]
  results$model<-"Model1"
  
  cat(yellow("Model1 +++++",pheno.col,"\n"))
  
  return.string=results
}

### ===================================================
### 2. Adonis adj age, BMI and center

result.model3 = foreach(i=3:length(col.test),.combine = rbind) %do%  {
  
  pheno.col=col.test[i]
  cov=c("AGE","BMI","center")
  tmp<-metadata[,c("NEW_ID",pheno.col,cov)]
  tmp<-na.omit(tmp)
  colnames(tmp)<-c("ID","test","cov1","cov2","cov3")
  
  dist.mat1<-dist_subset(dist.mat,tmp$ID)
  
  ad1 = adonis2(dist.mat1 ~ test+cov1+cov2+cov3,data = tmp,permutations = 1000,by = "margin")
  
  results = data.frame(ad1[1,])
  results$test = pheno.col
  results$N = dim(tmp)[1]
  results$model<-"Model3"
  
  cat(yellow("Model3 +++++",pheno.col,"\n"))
  
  return.string=results
}

## =======================================================================
## 3. Multivariable linear model Shannon diversity ~ phenotypes

result.model2 = foreach(i=3:length(col.test),.combine = rbind) %do%  {
  pheno.col=col.test[i]
  cov= c("AGE","BMI","center")
  tmp<-df.pheno[,c("diversity_shannon",pheno.col,cov)]
  colnames(tmp)<-c("alpha","test","cov1","cov2","cov3")
  N.sample<-dim(na.omit(tmp))[1]
  
  mod<-lm(formula = alpha~test+cov1+cov2+cov3,data = tmp)
  coef<-summary(mod)$coefficients
  ci<-confint(mod, rownames(coef)[2], level=0.95)
  result=data.frame(Model="Model2",col=pheno.col,BETA=coef[2,1],SE=coef[2,2],Pval=coef[2,4],
                    CI.low=ci[1],CI.high=ci[2],N=N.sample)
  cat(yellow("Model2 +++++",pheno.col,"\n"))
  
  return.string=result
}

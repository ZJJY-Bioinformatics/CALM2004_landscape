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

### ============================================================================================
### Linear Model of associations between taxa (species and genus) with phenotypes
Assoc_taxa=function(mtx.taxa,dat.pheno,tmp.taxa,infection,cov.fix){
  
  ## merge data
  tmp.data=merge(mtx.taxa[,tmp.taxa,drop=F],dat.pheno[,c(infection,cov.fix),drop=F],by="row.names",all=F)
  tmp.data=na.omit(tmp.data)
  
  if(is.null(cov.fix)){
    colnames(tmp.data)<-c("ID","Bacteria","Infection")
    formula<-as.formula(paste("Bacteria~", "Infection"))
    
  }else{
    ## generate formular
    colid<-LETTERS[seq(from=1,to=length(cov.fix))]
    colnames(tmp.data)<-c("ID","Bacteria","Infection",colid)
    formula<-as.formula(paste("Bacteria~", "Infection+",paste(colid, collapse="+")))
  }
  
  ## run model
  tmp.mod1=lm(formula = formula,data = tmp.data)
  coef<-summary(tmp.mod1)
  ## summary
  tmp.result=data.frame(Bacteria=tmp.taxa,Var=infection,
                        estimate=coef$coefficients[2,1],se=coef$coefficients[2,2],
                        t=coef$coefficients[2,3],Pval=coef$coefficients[2,4],
                        N=dim(tmp.data)[1])
  return(tmp.result)
}

### ============================================================================================
### Plot all associatiosn (Figure 2)

### panel 1 adonis results
p1<-ggplot(res.beta,aes(x=R2,y=reorder(test,R2),fill=sig))+
  geom_bar(stat = "identity")+
  facet_grid(group~.,scales = "free_y",space = "free")+
  scale_fill_manual(values=c("gray","gold1","#D93907"))+
  theme_classic()+ylab("")+
  theme(legend.position = "top")

### panel 2 alpha diversity results
p2<-ggplot(res.alpha1,aes(x=BETA,y=col,color=sig))+
  geom_point(mapping = aes(x=BETA),size=3,shape=15)+
  facet_grid(group~.,scales = "free_y",space = "free")+
  geom_errorbar(aes(x=BETA,xmin=CI.low, xmax=CI.high),width=0)+
  scale_color_manual(values=c("gray","gold1","#D93907"))+
  theme_classic()+ylab("")+
  geom_vline(xintercept = 0, linetype="dotted",size=1)+
  theme(legend.position = "top")

### panel 3 taxa association results
p3<-ggplot(res.plot,aes(x=Bacteria,y=Var,fill=estimate))+
  geom_tile(colour="gray")+
  scale_fill_gradient2(low="#0052A5", high = "#d11141", mid = "white", midpoint = 0)+
  #scale_fill_manual(values=c("#000099","#0033CC","#0000FF","#6666FF","#9999FF","white","#FFF3B2","#FEB24C","#FC4E2A","#E31A1C","#990000"))+
  geom_text(aes(label=sig), size=4,color="black") +
  theme_classic()+
  scale_x_discrete(limits = taxa_order)+
  labs(fill="Beta")+
  facet_grid(group~.,scales = "free_y",space = "free")+
  #coord_fixed(ratio = 1)+
  theme(#axis.text.y = element_text(hjust = 0),
    axis.text.x = element_text(angle = 90, hjust = 1),
    #panel.border = element_rect(colour = "black", fill=NA, size=0.8),
    legend.position = "top")+ylab("")+xlab("")


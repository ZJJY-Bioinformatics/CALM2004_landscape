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

### calculate balance score
numerator_taxa_index<-res.sig$taxa[which(res.sig$estimate<0)]
denominator_taxa_index<-res.sig$taxa[which(res.sig$estimate>0)]
ftbl<-t(taxa1)
ftbl[ftbl==0]<-1e-10

balance_df <- metadata %>% mutate(balance_value=NA)

for(sample in rownames(ftbl)){
  balance_df[sample,]$balance_value <- log(exp(mean(log(ftbl[sample,numerator_taxa_index])))) - log(exp(mean(log(ftbl[sample,denominator_taxa_index]))))  
}

### plot balance score results
p11<-ggplot(na.omit(balance_df[,c("G_healthy","balance_value")]),aes(x=G_healthy,y=balance_value,fill=G_healthy))+
  geom_violin(width=0.5)+
  geom_boxplot(alpha=0.6,width=0.2,fill="white")+
  stat_compare_means()+
  theme_classic()+
  theme(legend.position = "none")+
  scale_fill_nejm()+xlab(" ")+ylab("Balance score")

p12<-ggplot(na.omit(balance_df[,c("STI","balance_value")]),aes(x=STI,y=balance_value,fill=STI))+
  geom_violin(width=0.5)+
  geom_boxplot(alpha=0.6,width=0.2,fill="white")+
  stat_compare_means()+
  theme_classic()+
  theme(legend.position = "none")+
  scale_fill_nejm()+xlab(" ")+ylab("Balance score")

p13<-ggplot(na.omit(balance_df[,c("EXA","balance_value")]),aes(x=EXA,y=balance_value,fill=EXA))+
  geom_violin(width=0.5)+
  geom_boxplot(alpha=0.6,width=0.2,fill="white")+
  stat_compare_means()+
  theme_classic()+
  theme(legend.position = "none")+
  scale_fill_nejm()+xlab(" ")+ylab("Balance score")

#Depth correction
#Goal: Adjust the number of observed variants at the sites with low sequencing coverage
#Step 1: Calculate the prop of obs/exp for each coverage: 1-39 and above 40
#Step 2: Plot the relationship between obs/exp and sequencing coverage at a normal scale
#Step 3: Plot the relationship of the two at a log scale
#Step 4: Get the scaling factor for each sequencing coverage

#Step 1:
#For sites with sequencing coverage between 1-39:
rm(list=ls())
#loading the data with expected number of variants for sites with various sequencing coverage
load("./low_coverage_exp_variants.RData")

#loading the data with obs number of variants for sites with various sequencing coverage
load("./gnomad_exome_observed_low_coverage.RData")

low_coverage<-merge(coverage_exp,coverage_observed,by="median")

low_coverage$ratio<-low_coverage$key/low_coverage$exp
low_coverage<-subset(low_coverage,select=c(median,ratio))

#For sites with sequencing coverage above 40
load("./high_coverage_exp.RData")
load("./high_coverage_obs_count.RData")

high_coverage<-merge(all_pos_count,all_obs_count,by="median")

#Calculuate the global scaling factor 
scale_factor<-sum(high_coverage$count)/sum(high_coverage$exp)

high_coverage$ratio<-high_coverage$count/high_coverage$exp
high_coverage<-subset(high_coverage,select=c(median,ratio))

coverage_ratio<-rbind(low_coverage,high_coverage)
colnames(coverage_ratio)[1]<-"coverage"
coverage_ratio$ratio<-coverage_ratio$ratio/scale_factor

##Step 2:
#Plot the relationship between the coverage and the ratio
library(ggplot2)

A=ggplot(data=coverage_ratio,aes(x=coverage,y=ratio))+geom_point(color="darkgrey", size=4)+
  labs(y="Observed/Expected")+theme_classic()+coord_fixed(ratio=100)


##Step 3:
#Fit a linear relationship between log10(coverage) with ratio
coverage_ratio$log10_coverage<-log10(coverage_ratio$coverage)

log10coverage_ratio_lm<-lm(ratio~log10_coverage,data=coverage_ratio)
summary(log10coverage_ratio_lm)

B=ggplot(data = coverage_ratio,aes(x=log10_coverage,y=ratio))+geom_point(color="grey",size=4)+
  geom_smooth(method = "lm",formula=y~x, se = FALSE,color="black",size=1)+labs(y="Observed/Expected")+
  labs(x="log10(coverage)")+
  theme_classic()+coord_fixed(ratio=1.8)

library(ggpubr)
figure<-ggarrange(A,B,ncol=2,nrow=1)
figure

correction_factor_exp<-predict(log10coverage_ratio_lm,newdata=coverage_ratio)
correction_factor_obs<-1/correction_factor_exp
correction_factor<-data.frame(coverage=seq(1,39,1),cor_factor_exp=correction_factor_exp[1:39],cor_factor_obs=correction_factor_obs[1:39])


save(correction_factor,file="./coverage_correction_factor.RData")


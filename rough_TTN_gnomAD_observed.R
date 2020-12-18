#Purpose: Get the observed number of synonmous mutations in TTN from gnomAD exome
# Step 1: Include the sites with low coverage. Multiply with correction factor
# Step 2: Combine the results from high coverage and low coverage

#Step 1: 
rm(list=ls())
TTN_low_coverage_gnomad<-read.table("./rough_TTN_gnomad_low_coverage.txt",sep="\t",header=T)

#(1)Get the mutation code
get_column_index<-function(col_name,col_list){
  return(which(col_list==col_name))
}

methyl_level<-function(x){
  methy_level = NA
  if(is.na(x)){methy_level=0;return(methy_level)}
  if(x<0.2){methy_level = 0}
  else if(x>0.6){methy_level = 2}
  else {methy_level=1}
  return(methy_level)
}

get_sub_mutation<-function(x,df){
  col_list=colnames(df)
  met_level=get_column_index("methyl_level",col_list)
  context = get_column_index("context",col_list)
  ref = get_column_index("REF",col_list)
  alt = get_column_index("ALT",col_list)
  
  sub=NA
  sub=paste(x[context],x[ref],x[alt],sep="_")
  
  return(sub)
}

#Issue: some variants are in CpG sites but are not methylation mutation. These sites should have zero methlyation level
#Solution: for all the variants in data frame, if they are not methylation mutation, then set the methylation level to 0
CpG_mutation<-c("ACG_C_T","CCG_C_T","GCG_C_T","TCG_C_T","CGA_G_A","CGC_G_A","CGG_G_A","CGT_G_A")

check_methylation_level<-function(x,df){
  col_list = colnames(df)
  sub_code<-get_column_index("sub",col_list)
  met_level<-get_column_index("methyl_level",col_list)
  if(!is.element(x[sub_code],CpG_mutation)){
    x[met_level]=0
  }
  return(x[met_level])
}

#Transform the methyl_mean into methyl_level
TTN_low_coverage_gnomad$methyl_level<-sapply(TTN_low_coverage_gnomad$methyl_mean,methyl_level)

#Get the substitution code (i.e. context_ref_alt)
TTN_low_coverage_gnomad$sub<-apply(TTN_low_coverage_gnomad,1,function(x){get_sub_mutation(x,TTN_low_coverage_gnomad)})

#Check the methylation level of the variants, if it's not methylation mutation, the methylation level should be set to 0
TTN_low_coverage_gnomad$methyl_level<-apply(TTN_low_coverage_gnomad,1,function(x){check_methylation_level(x,TTN_low_coverage_gnomad)})

#Update the substitution code with the right mehtylation level: met_context_ref_alt
TTN_low_coverage_gnomad$mutation<-paste(TTN_low_coverage_gnomad$methyl_level,TTN_low_coverage_gnomad$sub,sep="_")


##Step 2: combine the number of observed variants from both high coverage sites and low coverage sites
TTN_high_coverage_gnomad<-read.table("./rough_TTN_gnomad_filtered.txt",sep="\t",header=T)


#(2)Get the mutation code
get_column_index<-function(col_name,col_list){
  return(which(col_list==col_name))
}

methyl_level<-function(x){
  methy_level = NA
  if(is.na(x)){methy_level=0;return(methy_level)}
  if(x<0.2){methy_level = 0}
  else if(x>0.6){methy_level = 2}
  else {methy_level=1}
  return(methy_level)
}

get_sub_mutation<-function(x,df){
  col_list=colnames(df)
  met_level=get_column_index("methyl_level",col_list)
  context = get_column_index("context",col_list)
  ref = get_column_index("REF",col_list)
  alt = get_column_index("ALT",col_list)
  
  sub=NA
  sub=paste(x[context],x[ref],x[alt],sep="_")
  
  return(sub)
}

#Issue: some variants are in CpG sites but are not methylation mutation. These sites should have zero methlyation level
#Solution: for all the variants in data frame, if they are not methylation mutation, then set the methylation level to 0
CpG_mutation<-c("ACG_C_T","CCG_C_T","GCG_C_T","TCG_C_T","CGA_G_A","CGC_G_A","CGG_G_A","CGT_G_A")

check_methylation_level<-function(x,df){
  col_list = colnames(df)
  sub_code<-get_column_index("sub",col_list)
  met_level<-get_column_index("methyl_level",col_list)
  if(!is.element(x[sub_code],CpG_mutation)){
    x[met_level]=0
  }
  return(x[met_level])
}

#Transform the methyl_mean into methyl_level
TTN_high_coverage_gnomad$methyl_level<-sapply(TTN_high_coverage_gnomad$methyl_mean,methyl_level)

#Get the substitution code (i.e. context_ref_alt)
TTN_high_coverage_gnomad$sub<-apply(TTN_high_coverage_gnomad,1,function(x){get_sub_mutation(x,TTN_high_coverage_gnomad)})

#Check the methylation level of the variants, if it's not methylation mutation, the methylation level should be set to 0
TTN_high_coverage_gnomad$methyl_level<-apply(TTN_high_coverage_gnomad,1,function(x){check_methylation_level(x,TTN_low_coverage_gnomad)})

#Update the substitution code with the right mehtylation level: met_context_ref_alt
TTN_high_coverage_gnomad$mutation<-paste(TTN_high_coverage_gnomad$methyl_level,TTN_high_coverage_gnomad$sub,sep="_")


rough_TTN_obs<-rbind(TTN_low_coverage_gnomad,TTN_high_coverage_gnomad)
save(rough_TTN_obs,file="./rough_TTN_obs.RData")

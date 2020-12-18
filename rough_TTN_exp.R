#Purpose: Get the expected number of mutations in the rough TTN sites
#Step 1: get the total possible synonmous mutations for each substitution, each context and each methylation level
# I could combine the both the low coverage and high coverage sites into one to do that.
#Step 2: times the correponding expected subtition rate estimated from calibration on syn variants

#How to get the all possible syn sites:
#From the synthetic file of ExAC, extract the ones in the exon position with the corresponding conseq
setwd("~/Desktop/OneDrive/OneDrive - Imperial College London/TTN/circularRNA/Constraint_analysis")
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
  context = get_column_index("context_x",col_list)
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


high_coverage<-read.table("./rough_all_possible_snp_TTN_methy_filtered.txt",sep="\t",header = T)
low_coverage<-read.table("./rough_TTN_all_possible_low_coverage_methy_filtered.txt",sep="\t",header=T)

TTN_pos<-rbind(high_coverage,low_coverage)

#Transform the methyl_mean into methyl_level
TTN_pos$methyl_level<-sapply(TTN_pos$methyl_mean,methyl_level)

#Get the substitution code (i.e. context_ref_alt)
TTN_pos$sub<-apply(TTN_pos,1,function(x){get_sub_mutation(x,TTN_pos)})

#Check the methylation level of the variants, if it's not methylation mutation, the methylation level should be set to 0
TTN_pos$methyl_level<-apply(TTN_pos,1,function(x){check_methylation_level(x,TTN_pos)})

#Update the substitution code with the right mehtylation level: met_context_ref_alt
TTN_pos$mutation<-paste(TTN_pos$methyl_level,TTN_pos$sub,sep="_")

#add predicted probability of mutation
load("./all_predicted_prop.RData")
TTN_exp<-merge(TTN_pos,predicted_prop,by="mutation",all.x = TRUE)
save(TTN_exp,file="./rough_TTN_exp.RData")



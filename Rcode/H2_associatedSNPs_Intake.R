library(dplyr)
library(tidyr)
library(scales)
library(readxl)
library(Rmpfr)

ldchecking0.2 <- read.csv("/media/cailu/Expansion1/ukb/bfilemaf01/variants_generegionmaf0.1/FoodsIntakes/ldfoodintake/ld0.2.ld", sep="")
ldchecking0.5 <- read.csv("/media/cailu/Expansion1/ukb/bfilemaf01/variants_generegionmaf0.1/FoodsIntakes/ldfoodintake/ld0.5.ld", sep="")
ldchecking0.8 <- read.csv("/media/cailu/Expansion1/ukb/bfilemaf01/variants_generegionmaf0.1/FoodsIntakes/ldfoodintake/ld0.8.ld", sep="")

ass <- read.csv("/media/cailu/Expansion1/ukb/bfilemaf01/variants_generegionmaf0.1/FoodsIntakers20241028/ExtraAssociations_foodintake20241223.csv")

df <- read.csv("/media/cailu/Expansion1/ukb/bfilemaf01/variants_generegionmaf0.1/FoodsIntakers20241028/all_geno_pheon_data_20241223.csv")


ass$marker[ass$marker=="rs10743938"]<-"rs10743938...888"


names(df)[names(df)=="rs2074464...235"]<-"rs2074464...235"
names(df)[names(df)=="rs10743938...897"]<-"rs10743938...888"
names(df)[names(df)=="rs7130086...518"]<-"rs7130086"
names(df)[names(df)=="rs605734...717"]<-"rs605734"
names(df)[names(df)=="rs7124871...594"]<-"rs7124871"
names(df)[names(df)=="rs17496724...568"]<-"rs17496724"
names(df)[names(df)=="rs61890419...546"]<-"rs61890419"
names(df)[names(df)=="rs12224086...636"]<-"rs12224086"
names(df)[names(df)=="rs2173236...1045"]<-"rs2173236...1035"



##v2
ldchecking0.2 <- read.csv("/media/cailu/Expansion1/ukb/bfilemaf01/LD/ldchecking0.2.ld", sep="")%>%
  mutate(SNP_A = recode(SNP_A, "rs10743938.1" = "rs10743938...888", "rs2074464"= "rs2074464...235", 'rs2173236.1'= 'rs2173236...1035', "rs7130086.1"="rs7130086...515"))%>%
  mutate(SNP_B = recode(SNP_B, "rs10743938.1" = "rs10743938...888", "rs2074464"= "rs2074464...235", 'rs2173236.1'= 'rs2173236...1035', "rs7130086.1"="rs7130086...515"))

ldchecking0.5 <- read.csv("/media/cailu/Expansion1/ukb/bfilemaf01/LD/ldchecking0.5.ld", sep="")%>%
  mutate(SNP_A = recode(SNP_A, "rs10743938.1" = "rs10743938...888", "rs2074464"= "rs2074464...235", 'rs2173236.1'= 'rs2173236...1035', "rs7130086.1"="rs7130086...515"))%>%
  mutate(SNP_B = recode(SNP_B, "rs10743938.1" = "rs10743938...888", "rs2074464"= "rs2074464...235", 'rs2173236.1'= 'rs2173236...1035', "rs7130086.1"="rs7130086...515"))

ldchecking0.8 <- read.csv("/media/cailu/Expansion1/ukb/bfilemaf01/LD/ldchecking0.8.ld", sep="")%>%
  mutate(SNP_A = recode(SNP_A, "rs10743938.1" = "rs10743938...888", "rs2074464"= "rs2074464...235", 'rs2173236.1'= 'rs2173236...1035', "rs7130086.1"="rs7130086...515"))%>%
  mutate(SNP_B = recode(SNP_B, "rs10743938.1" = "rs10743938...888", "rs2074464"= "rs2074464...235", 'rs2173236.1'= 'rs2173236...1035', "rs7130086.1"="rs7130086...515"))



df$sex<-as.numeric(df$sex)

trait<-unique(ass$trait)
result<-data.frame(trait=character(0), r2=numeric(0), p=character(0), N=numeric(0), Sample.szie=numeric(0))
for (i in trait){
  temp<-data.frame(trait=NA, r2=NA, p=NA, N=NA,Sample.szie=NA)
  list_snps<-unique(ass[ass$trait==i,]$marker)
  ExVar <- toString(paste(list_snps, "+ ", collapse = ""))
  ExVar <- substr(ExVar, 1, nchar(ExVar)-3)  
  dat<-df[, c(i,list_snps, "sex", "age","PC1","PC2" ,"PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")]
  dat<-dat%>% drop_na()
  
  m0 <- lm(formula(paste(i,"~sex+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")),dat)
  m1 <- lm(formula(paste(i,"~",ExVar,"+sex+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")),dat)
  m0s <- summary(m0)
  m1s <- summary(m1)
  temp$trait<-i
  temp$r2 <- m1s$r.squared - m0s$r.squared
  m2ll <- as.numeric((-2)*(logLik(m0)-logLik(m1)))
  dfn<-length(coef(m1))-length(coef(m0))
  log_p_value<-pchisq(q=m2ll,df = dfn, lower.tail = F,log.p = TRUE) # calculate p-value from -2*loglikelihood
  log_p_value_mpfr <- mpfr(log_p_value, precBits = 256)
  temp$p<-formatMpfr(exp(log_p_value_mpfr),  digits = 3, scientific = TRUE)
  temp$N<-length(list_snps)
  temp$Sample.szie<-nrow(dat)
 result<-bind_rows(result, temp)
}


##clumping snp r2<0.2
result.2<-data.frame(trait=character(0), r2.2=numeric(0), p.2=character(0), N.2=numeric(0), Sample.szie.2=numeric(0))
for (i in trait){
  temp<-data.frame(trait=NA, r2.2=NA, p.2=NA, N.2=NA,Sample.szie.2=NA)
  subass<-ass[ass$trait==i,]
  list_snps<-unique(subass[order(-subass$p.value),]$marker)
  x<-ldchecking0.2[ldchecking0.2$SNP_A%in%list_snps&ldchecking0.2$SNP_B%in%list_snps,]
  list_snps_paired<-list_snps[list_snps%in%x$SNP_A|list_snps%in%x$SNP_B]
  snps_to_remove <- list_snps_paired[c(TRUE, FALSE)]
  list_snps_clumped<-list_snps[!list_snps%in%snps_to_remove]
  
  ExVar <- toString(paste(list_snps_clumped, "+ ", collapse = ""))
  ExVar <- substr(ExVar, 1, nchar(ExVar)-3)  
  dat<-df[, c(i,list_snps, "sex", "age","PC1","PC2" ,"PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")]
  dat<-dat%>%drop_na()
  
  m0 <- lm(formula(paste(i,"~sex+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")),dat)
  m1 <- lm(formula(paste(i,"~",ExVar,"+sex+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")),dat)
  m0s <- summary(m0)
  m1s <- summary(m1)
  temp$trait<-i
  temp$r2.2 <- m1s$r.squared - m0s$r.squared
  m2ll <- as.numeric((-2)*(logLik(m0)-logLik(m1)))
  dfn<-length(coef(m1))-length(coef(m0))
  log_p_value<- pchisq(q=m2ll,df = dfn, lower.tail = F,log.p = TRUE) # calculate p-value from -2*loglikelihood
  log_p_value_mpfr <- mpfr(log_p_value, precBits = 256)
  temp$p.2<-formatMpfr(exp(log_p_value_mpfr),  digits = 3, scientific = TRUE)
  temp$N.2<-length(list_snps_clumped)
  temp$Sample.szie.2<-nrow(dat)
  result.2<-bind_rows(result.2, temp)
}

##clumping snp r2<0.5
result.5<-data.frame(trait=character(0), r2.5=numeric(0), p.5=character(0), N.5=numeric(0), Sample.szie.5=numeric(0))
for (i in trait){
  temp<-data.frame(trait=NA, r2.5=NA, p.5=NA, N.5=NA,Sample.szie.5=NA)
  subass<-ass[ass$trait==i,]
  list_snps<-unique(subass[order(-subass$p.value),]$marker)
  x<-ldchecking0.5[ldchecking0.5$SNP_A%in%list_snps&ldchecking0.5$SNP_B%in%list_snps,]
  list_snps_paired<-list_snps[list_snps%in%x$SNP_A|list_snps%in%x$SNP_B]
  snps_to_remove <- list_snps_paired[c(TRUE, FALSE)]
  list_snps_clumped<-list_snps[!list_snps%in%snps_to_remove]
  
  ExVar <- toString(paste(list_snps_clumped, "+ ", collapse = ""))
  ExVar <- substr(ExVar, 1, nchar(ExVar)-3)  
  dat<-df[, c(i,list_snps, "sex", "age","PC1","PC2" ,"PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")]
  dat<-dat%>%
    drop_na()
  
  m0 <- lm(formula(paste(i,"~sex+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")),dat)
  m1 <- lm(formula(paste(i,"~",ExVar,"+sex+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")),dat)
  m0s <- summary(m0)
  m1s <- summary(m1)
  temp$trait<-i
  temp$r2.5 <- m1s$r.squared - m0s$r.squared
  m2ll <- as.numeric((-2)*(logLik(m0)-logLik(m1)))
  dfn<-length(coef(m1))-length(coef(m0))
  log_p_value<- pchisq(q=m2ll,df = dfn, lower.tail = F,log.p = TRUE) # calculate p-value from -2*loglikelihood
  log_p_value_mpfr <- mpfr(log_p_value, precBits = 256)
  temp$p.5<-formatMpfr(exp(log_p_value_mpfr),  digits = 3, scientific = TRUE)
  temp$N.5<-length(list_snps_clumped)
  temp$Sample.szie.5<-nrow(dat)
  result.5<-bind_rows(result.5, temp)
}

##clumping snp r2<0.8
result.8<-data.frame(trait=character(0), r2.8=numeric(0), p.8=character(0), N.8=numeric(0), Sample.szie.8=numeric(0))
for (i in trait){
  temp<-data.frame(trait=NA, r2.8=NA, p.8=NA, N.8=NA,Sample.szie.8=NA)
  subass<-ass[ass$trait==i,]
  list_snps<-unique(subass[order(-subass$p.value),]$marker)
  x<-ldchecking0.8[ldchecking0.8$SNP_A%in%list_snps&ldchecking0.8$SNP_B%in%list_snps,]
  list_snps_paired<-list_snps[list_snps%in%x$SNP_A|list_snps%in%x$SNP_B]
  snps_to_remove <- list_snps_paired[c(TRUE, FALSE)]
  list_snps_clumped<-list_snps[!list_snps%in%snps_to_remove]
  
  ExVar <- toString(paste(list_snps_clumped, "+ ", collapse = ""))
  ExVar <- substr(ExVar, 1, nchar(ExVar)-3)  
  dat<-df[, c(i,list_snps, "sex", "age","PC1","PC2" ,"PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")]
  dat<-dat%>%
    drop_na()
  
  m0 <- lm(formula(paste(i,"~sex+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")),dat)
  m1 <- lm(formula(paste(i,"~",ExVar,"+sex+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")),dat)
  m0s <- summary(m0)
  m1s <- summary(m1)
  temp$trait<-i
  temp$r2.8 <- m1s$r.squared - m0s$r.squared
  m2ll <- as.numeric((-2)*(logLik(m0)-logLik(m1)))
  dfn<-length(coef(m1))-length(coef(m0))
  log_p_value<- pchisq(q=m2ll,df = dfn, lower.tail = F,log.p = TRUE) # calculate p-value from -2*loglikelihood
  log_p_value_mpfr <- mpfr(log_p_value, precBits = 256)
  temp$p.8<-formatMpfr(exp(log_p_value_mpfr),  digits = 3, scientific = TRUE)
  temp$N.8<-length(list_snps_clumped)
  temp$Sample.szie.8<-nrow(dat)
  result.8<-bind_rows(result.8, temp)
}


results_final<-merge(result, merge(result.8, merge(result.5, result.2, by="trait"), by="trait"), by="trait")

write.csv(results_final, "/media/cailu/Expansion1/ukb/bfilemaf01/variants_generegionmaf0.1/FoodsIntakers20241028/Foodintake_SNP_h2_20241223.csv", row.names=F)





pacman::p_load(ukbtools,vcfR,readxl, broom, dplyr, tidyr, scales, caret, ggstatsplot, tibble,lm.beta)
df_liking<-read.csv("/media/cailu/Expansion1/ukb/bfilemaf01/variants_generegionmaf0.1/geno_pheno20240605.csv")
df_intake<-read.csv("/media/cailu/Seagate Backup Plus Drive1/data_20230311/UKB_geno/bfilemaf0.1/FoodsIntakes/likingSpecificfoodIntakes.csv")

##Onion intake
d_onion_intake<-df_intake[names(df_intake)%in%c("eid", "onion_intake_f104260")]
df<-merge(d_onion_intake,df_liking, by="eid" )
df1<-df[c("onion_intake_f104260","liking_for_onions_f20693_0_0", "sex","age",
         "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
##
model<-lm(onion_intake_f104260~liking_for_onions_f20693_0_0+sex+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data=df1)

summary(model)

cor(df1[,c("onion_intake_f104260", "liking_for_onions_f20693_0_0")], use = "complete.obs", method = "pearson")

cor.test(df$onion_intake_f104260, df$liking_for_onions_f20693_0_0, method = "pearson")

##
onion<-names(my_ukb_data)[grepl("onion", names(my_ukb_data))]
d<-my_ukb_data[c("eid", onion)]

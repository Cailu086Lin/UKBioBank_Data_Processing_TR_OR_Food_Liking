library(biomaRt)
library(dplyr)
library(readxl)
library(tidyr)

listMarts()

snps<-useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp") #GRCh38
snps<-useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_snp") #GRCh37

listDatasets(snps)

snps<-useDataset("hsapiens_snp",mart = snps)
listFilters(snps)


snps_info<-getBM(filters="chromosomal_region",attributes=c("refsnp_id","chr_name","allele","minor_allele","minor_allele_freq",
                                                           "clinical_significance","polyphen_prediction","polyphen_score","sift_prediction","sift_score"),
                 values=coordinates, mart=snps)


##sort genes
#gencode.v43lift37.metadata.HGNC <- read.delim("/media/cailu/Seagate Backup Plus Drive/data_20230311/results3/gencode.v43lift37.metadata.HGNC.gz", header=FALSE)
#genes<-unique(gencode.v43lift37.metadata.HGNC$V2[grepl("^TAS1R|^TAS2R|^OR",  gencode.v43lift37.metadata.HGNC$V2)])


gencode.v19.metadata.HGNC <- read.delim("/media/cailu/Seagate Backup Plus Drive/data_20230311/results3/gencode.v19.metadata.HGNC.gz", header=FALSE)
genes2<-unique(gencode.v19.metadata.HGNC $V2[grepl("^TAS1R|^TAS2R|^OR",  gencode.v19.metadata.HGNC $V2)])
drop<-genes2[grepl("*P", genes2)]
genes<-sort(genes2[!(genes2%in%drop)])


####get gene chromsomal regions.

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "grch37.ensembl.org",path="/biomart/martservice")


# Retrieve chromosome regions for the genes
gene_region <- data.frame()

for (gene in genes) {
  # Get the gene information
  tryCatch({
  gene_info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
                     filters = "hgnc_symbol", 
                     values = gene, 
                     mart = ensembl)
  .lock2(dbfile, exclusive = TRUE)
  }, error = function(e) {
  })
  if(nrow(gene_info)>0){
  # Get the missense variants for the gene
     gene_region1 <- getBM(attributes = c("chromosome_name", "start_position", "end_position"), 
                                 filters = "ensembl_gene_id",
                                 values = gene_info$ensembl_gene_id,
                                 mart = ensembl)
      
      # Append the variants to the list
      gene_region1$gene <- gene
      gene_region1<-gene_region1[!grepl("[^[:digit:]]", gene_region1$chromosome_name),]
      gene_region1$chromosome_name <-as.integer(gene_region1$chromosome_name)
      gene_region<-bind_rows(gene_region, gene_region1)
      
      }else{
          gene_region<-gene_region
      }
}

gene_region$key<-paste0(gene_region$chromosome_name,":", gene_region$start_position,"_", gene_region$end_position)

write.csv(gene_region,"/media/cailu/Seagate Backup Plus Drive/data_20230311/results3/gene_regions.csv", row.names=F)

# Print the extracted variants
for (gene in genes) {
  cat("Missense variants for", gene, ":\n")
  print(variants[[gene]])
  cat("\n")
}


###
gene_regions <- read_csv("/media/cailu/Seagate Backup Plus Drive/data_20230311/results3/gene_regions.csv")
gene_regions<-gene_regions[!gene_regions$gene%in%gene_regions$gene[grepl("^ORA|^ORC|^ORM", gene_regions$gene)],] ##425 genes

##Chr2
g <- read_excel("/media/cailu/Seagate Backup Plus Drive/data_20230311/results3/gene_regions.xlsx")
c2 <- read.delim("/media/cailu/Seagate Backup Plus Drive/data_20230311/UKB_geno/bfilemaf0.1/c2flt.bim", header=FALSE)
g1<-g[g$chromosome_name==2,]
snp<-data.frame()
for (i in 1: nrow(g1)){
   snpx<-c2[c2$V4>g1[i,]$start_position &c2$V4<g1[i,]$end_position,] 
   snpx$gene<-g1[i,]$gene
   snp<-bind_rows(snp, snpx)
}
v2<-merge(g1, snp, by="gene")

write.csv(v2, "/media/cailu/Seagate Backup Plus Drive1/data_20230311/UKB_geno/bfilemaf0.1/variants_generegionmaf0.1/variants_c2.csv", row.names=F)

v2_ana <- read_excel("ann_c2.xlsx")%>%
  select(variant,dbSNP, Ref, Alt, EUR)%>%
  drop_na(dbSNP)%>%
  distinct(variant, .keep_all=T)
names(v2_ana)[c(2, 5)]<-c("dbSNP.annotation", "EUR_freq")
v2a<-merge(v2,  v2_ana, by.x="V2", by.y="variant", all.x=T)%>%
  filter(dbSNP.annotation=="missense"|dbSNP.annotation=="frameshift"|dbSNP.annotation=="nonsense")

v2b<-bind_rows(v2a, g1[!(g1$gene%in%v2a$gene),])
write.csv(v2b, "chr2_missese.csv", row.names=F)

##chr1

c2 <- read.delim("/media/cailu/Seagate Backup Plus Drive/data_20230311/UKB_geno/bfilemaf0.1/c1flt.bim", header=FALSE)
g1<-g[g$chromosome_name==1,]
snp<-data.frame()
for (i in 1: nrow(g1)){
  snpx<-c2[c2$V4>g1[i,]$start_position &c2$V4<g1[i,]$end_position,] 
  if(nrow(snpx)>0){
  snpx$gene<-g1[i,]$gene
  snp<-bind_rows(snp, snpx)
  }else{snpx[1,]<-c(NA, NA, NA,NA, NA, NA)
  snpx$gene<-g1[i,]$gene
  snp<-bind_rows(snp, snpx)
  }
}
v2<-merge(g1, snp, by="gene")

write.csv(v2, "/media/cailu/Seagate Backup Plus Drive1/data_20230311/UKB_geno/bfilemaf0.1/variants_generegionmaf0.1/variants_c1.csv", row.names=F)

v1_ana <- read_excel("ann_c1A.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v1_ana1 <- read_excel("ann_c1A1.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v1_ana2 <- read_excel("ann_c1A2.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v1_ana3 <- read_excel("ann_c1A3.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v1_ana4 <- read_excel("ann_c1A4.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v1_ana<-bind_rows(v1_ana,v1_ana1,v1_ana2,v1_ana3,v1_ana4)
names(v1_ana)[c(2, 5)]<-c("dbSNP.annotation", "EUR_freq")
v2a<-merge(v2,  v1_ana, by.x="V2", by.y="variant", all.x=T)%>%
  filter(dbSNP.annotation=="missense"|dbSNP.annotation=="frameshift"|dbSNP.annotation=="nonsense")

v2b<-bind_rows(v2a, g1[!(g1$gene%in%v2a$gene),])
write.csv(v2b, "chr1_missese.csv", row.names=F)


##chr5

c2 <- read.delim("/media/cailu/Seagate Backup Plus Drive/data_20230311/UKB_geno/bfilemaf0.1/c5flt.bim", header=FALSE)
g1<-g[g$chromosome_name==5,]
snp<-data.frame()
for (i in 1: nrow(g1)){
  snpx<-c2[c2$V4>g1[i,]$start_position &c2$V4<g1[i,]$end_position,] 
  if(nrow(snpx)>0){
    snpx$gene<-g1[i,]$gene
    snp<-bind_rows(snp, snpx)
  }else{snpx[1,]<-c(NA, NA, NA,NA, NA, NA)
  snpx$gene<-g1[i,]$gene
  snp<-bind_rows(snp, snpx)
  }
}
v2<-merge(g1, snp, by="gene")

write.csv(v2, "/media/cailu/Seagate Backup Plus Drive1/data_20230311/UKB_geno/bfilemaf0.1/variants_generegionmaf0.1/variants_c5.csv", row.names=F)
v5_ana <- read_excel("ann_c5.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
names(v5_ana)[c(2, 5)]<-c("dbSNP.annotation", "EUR_freq")
v2a<-merge(v2,  v5_ana, by.x="V2", by.y="variant", all.x=T)%>%
  filter(dbSNP.annotation=="missense"|dbSNP.annotation=="frameshift"|dbSNP.annotation=="nonsense")

v2b<-bind_rows(v2a, g1[!(g1$gene%in%v2a$gene),])
write.csv(v2b, "chr5_missese.csv", row.names=F)




##chr7

c2 <- read.delim("/media/cailu/Seagate Backup Plus Drive/data_20230311/UKB_geno/bfilemaf0.1/c7flt.bim", header=FALSE)
g1<-g[g$chromosome_name==7,]
snp<-data.frame()
for (i in 1: nrow(g1)){
  snpx<-c2[c2$V4>g1[i,]$start_position &c2$V4<g1[i,]$end_position,] 
  if(nrow(snpx)>0){
    snpx$gene<-g1[i,]$gene
    snp<-bind_rows(snp, snpx)
  }else{snpx[1,]<-c(NA, NA, NA,NA, NA, NA)
  snpx$gene<-g1[i,]$gene
  snp<-bind_rows(snp, snpx)
  }
}
v2<-merge(g1, snp, by="gene")

write.csv(v2, "/media/cailu/Seagate Backup Plus Drive1/data_20230311/UKB_geno/bfilemaf0.1/variants_generegionmaf0.1/variants_c7.csv", row.names=F)
v7_ana <- read_excel("ann_c7.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
names(v7_ana)[c(2, 5)]<-c("dbSNP.annotation", "EUR_freq")
v2a<-merge(v2,  v7_ana, by.x="V2", by.y="variant", all.x=T)%>%
  filter(dbSNP.annotation=="missense"|dbSNP.annotation=="frameshift"|dbSNP.annotation=="nonsense")

v2b<-bind_rows(v2a, g1[!(g1$gene%in%v2a$gene),])
write.csv(v2b, "chr7_missese.csv", row.names=F)


##chr10

c2 <- read.delim("/media/cailu/Seagate Backup Plus Drive/data_20230311/UKB_geno/bfilemaf0.1/c10flt.bim", header=FALSE)
g1<-g[g$chromosome_name==10,]
snp<-data.frame()
for (i in 1: nrow(g1)){
  snpx<-c2[c2$V4>g1[i,]$start_position &c2$V4<g1[i,]$end_position,] 
  if(nrow(snpx)>0){
    snpx$gene<-g1[i,]$gene
    snp<-bind_rows(snp, snpx)
  }else{snpx[1,]<-c(NA, NA, NA,NA, NA, NA)
  snpx$gene<-g1[i,]$gene
  snp<-bind_rows(snp, snpx)
  }
}
v2<-merge(g1, snp, by="gene")

write.csv(v2, "/media/cailu/Seagate Backup Plus Drive1/data_20230311/UKB_geno/bfilemaf0.1/variants_generegionmaf0.1/variants_c10.csv", row.names=F)
v10_ana <- read_excel("ann_c10.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, EUR)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
names(v10_ana)[c(2, 5)]<-c("dbSNP.annotation", "EUR_freq")
v2a<-merge(v2,  v10_ana, by.x="V2", by.y="variant", all.x=T)%>%
  filter(dbSNP.annotation=="missense"|dbSNP.annotation=="frameshift"|dbSNP.annotation=="nonsense")

v2b<-bind_rows(v2a, g1[!(g1$gene%in%v2a$gene),])
write.csv(v2b, "chr10_missese.csv", row.names=F)


##chr11
c2 <- read.delim("/media/cailu/Seagate Backup Plus Drive/data_20230311/UKB_geno/bfilemaf0.1/c11flt.bim", header=FALSE)
g1<-g[g$chromosome_name==11,]
snp<-data.frame()
for (i in 1: nrow(g1)){
  snpx<-c2[c2$V4>g1[i,]$start_position &c2$V4<g1[i,]$end_position,] 
  if(nrow(snpx)>0){
    snpx$gene<-g1[i,]$gene
    snp<-bind_rows(snp, snpx)
  }else{snpx[1,]<-c(NA, NA, NA,NA, NA, NA)
  snpx$gene<-g1[i,]$gene
  snp<-bind_rows(snp, snpx)
  }
}
v2<-merge(g1, snp, by="gene")

write.csv(v2, "/media/cailu/Seagate Backup Plus Drive1/data_20230311/UKB_geno/bfilemaf0.1/variants_generegionmaf0.1/variants_c11.csv", row.names=F)
v11_ana1 <- read_excel("ann_c11A1.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v11_ana2 <- read_excel("ann_c11A2.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v11_ana3 <- read_excel("ann_c11A3.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v11_ana4 <- read_excel("ann_c11A4.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v11_ana5 <- read_excel("ann_c11A5.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v11_ana6 <- read_excel("ann_c11A6.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v11_ana7 <- read_excel("ann_c11A7.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v11_ana8 <- read_excel("ann_c11A8.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v11_ana9 <- read_excel("ann_c11A9.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)

v11_ana<-bind_rows(v11_ana1,v11_ana2,v11_ana3,v11_ana4,v11_ana5,v11_ana6,v11_ana7,v11_ana8,v11_ana9)

names(v11_ana)[c(2, 5)]<-c("dbSNP.annotation", "EUR_freq")
v2a<-merge(v2,  v11_ana, by.x="V2", by.y="variant", all.x=T)%>%
  filter(dbSNP.annotation=="missense"|dbSNP.annotation=="frameshift"|dbSNP.annotation=="nonsense")

v2b<-bind_rows(v2a, g1[!(g1$gene%in%v2a$gene),])
write.csv(v2b, "chr11_missese.csv", row.names=F)

##chr12
c2 <- read.delim("/media/cailu/Seagate Backup Plus Drive/data_20230311/UKB_geno/bfilemaf0.1/c12flt.bim", header=FALSE)
g1<-g[g$chromosome_name==12,]
snp<-data.frame()
for (i in 1: nrow(g1)){
  snpx<-c2[c2$V4>g1[i,]$start_position &c2$V4<g1[i,]$end_position,] 
  if(nrow(snpx)>0){
    snpx$gene<-g1[i,]$gene
    snp<-bind_rows(snp, snpx)
  }else{snpx[1,]<-c(NA, NA, NA,NA, NA, NA)
  snpx$gene<-g1[i,]$gene
  snp<-bind_rows(snp, snpx)
  }
}
v2<-merge(g1, snp, by="gene")

write.csv(v2, "/media/cailu/Seagate Backup Plus Drive1/data_20230311/UKB_geno/bfilemaf0.1/variants_generegionmaf0.1/variants_c12.csv", row.names=F)
v12_ana1 <- read_excel("ann_c12A1.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v12_ana2 <- read_excel("ann_c12A2.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v12_ana3 <- read_excel("ann_c12A3.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v12_ana4 <- read_excel("ann_c12A4.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v12_ana5 <- read_excel("ann_c12A5.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v12_ana6 <- read_excel("ann_c12A6.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v12_ana7 <- read_excel("ann_c12A7.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v12_ana8 <- read_excel("ann_c12A8.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v12_ana9 <- read_excel("ann_c12A9.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v12_ana10 <- read_excel("ann_c12A10.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v12_ana11 <- read_excel("ann_c12A11.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)

v12_ana<-bind_rows(v12_ana1, v12_ana2,v12_ana3,v12_ana4,v12_ana5,v12_ana6,v12_ana7,v12_ana8,v12_ana9,v12_ana10,v12_ana11)
names(v12_ana)[c(2, 5)]<-c("dbSNP.annotation", "EUR_freq")

v2a<-merge(v2,  v12_ana, by.x="V2", by.y="variant", all.x=T)%>%
  filter(dbSNP.annotation=="missense"|dbSNP.annotation=="frameshift"|dbSNP.annotation=="nonsense")

v2b<-bind_rows(v2a, g1[!(g1$gene%in%v2a$gene),])
write.csv(v2b, "chr12_missese.csv", row.names=F)




##chr17
c2 <- read.delim("/media/cailu/Seagate Backup Plus Drive/data_20230311/UKB_geno/bfilemaf0.1/c17flt.bim", header=FALSE)
g1<-g[g$chromosome_name==17,]
snp<-data.frame()
for (i in 1: nrow(g1)){
  snpx<-c2[c2$V4>g1[i,]$start_position &c2$V4<g1[i,]$end_position,] 
  if(nrow(snpx)>0){
    snpx$gene<-g1[i,]$gene
    snp<-bind_rows(snp, snpx)
  }else{snpx[1,]<-c(NA, NA, NA,NA, NA, NA)
  snpx$gene<-g1[i,]$gene
  snp<-bind_rows(snp, snpx)
  }
}
v2<-merge(g1, snp, by="gene")

write.csv(v2, "/media/cailu/Seagate Backup Plus Drive1/data_20230311/UKB_geno/bfilemaf0.1/variants_generegionmaf0.1/variants_c17.csv", row.names=F)
v17_ana <- read_excel("ann_c17.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
names(v17_ana)[c(2, 5)]<-c("dbSNP.annotation", "EUR_freq")
v2a<-merge(v2,  v17_ana, by.x="V2", by.y="variant", all.x=T)%>%
  filter(dbSNP.annotation=="missense"|dbSNP.annotation=="frameshift"|dbSNP.annotation=="nonsense")

v2b<-bind_rows(v2a, g1[!(g1$gene%in%v2a$gene),])
write.csv(v2b, "chr17_missese.csv", row.names=F)


##chr19
c2 <- read.delim("/media/cailu/Seagate Backup Plus Drive/data_20230311/UKB_geno/bfilemaf0.1/c19flt.bim", header=FALSE)
g1<-g[g$chromosome_name==19,]
snp<-data.frame()
for (i in 1: nrow(g1)){
  snpx<-c2[c2$V4>g1[i,]$start_position &c2$V4<g1[i,]$end_position,] 
  if(nrow(snpx)>0){
    snpx$gene<-g1[i,]$gene
    snp<-bind_rows(snp, snpx)
  }else{snpx[1,]<-c(NA, NA, NA,NA, NA, NA)
  snpx$gene<-g1[i,]$gene
  snp<-bind_rows(snp, snpx)
  }
}
v2<-merge(g1, snp, by="gene")

write.csv(v2, "/media/cailu/Seagate Backup Plus Drive1/data_20230311/UKB_geno/bfilemaf0.1/variants_generegionmaf0.1/variants_c19.csv", row.names=F)
v19_ana <- read_excel("ann_c19.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
names(v19_ana)[c(2, 5)]<-c("dbSNP.annotation", "EUR_freq")
v2a<-merge(v2,  v19_ana, by.x="V2", by.y="variant", all.x=T)%>%
  filter(dbSNP.annotation=="missense"|dbSNP.annotation=="frameshift"|dbSNP.annotation=="nonsense")

v2b<-bind_rows(v2a, g1[!(g1$gene%in%v2a$gene),])
write.csv(v2b, "chr19_missese.csv", row.names=F)

##chr6
c2 <- read.delim("/media/cailu/Seagate Backup Plus Drive/data_20230311/UKB_geno/bfilemaf0.1/c6flt.bim", header=FALSE)
g1<-g[g$chromosome_name==6,]
snp<-data.frame()
for (i in 1: nrow(g1)){
  snpx<-c2[c2$V4>g1[i,]$start_position &c2$V4<g1[i,]$end_position,] 
  if(nrow(snpx)>0){
    snpx$gene<-g1[i,]$gene
    snp<-bind_rows(snp, snpx)
  }else{snpx[1,]<-c(NA, NA, NA,NA, NA, NA)
  snpx$gene<-g1[i,]$gene
  snp<-bind_rows(snp, snpx)
  }
}
v2<-merge(g1, snp, by="gene")

write.csv(v2, "/media/cailu/Seagate Backup Plus Drive1/data_20230311/UKB_geno/bfilemaf0.1/variants_generegionmaf0.1/variants_c6.csv", row.names=F)

v6_ana1 <- read_excel("ann_c6A1.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v6_ana2 <- read_excel("ann_c6A2.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
v6_ana<-bind_rows(v6_ana1, v6_ana2)
names(v6_ana)[c(2, 5)]<-c("dbSNP.annotation", "EUR_freq")
v2a<-merge(v2,  v6_ana, by.x="V2", by.y="variant", all.x=T)%>%
  filter(dbSNP.annotation=="missense"|dbSNP.annotation=="frameshift"|dbSNP.annotation=="nonsense")

v2b<-bind_rows(v2a, g1[!(g1$gene%in%v2a$gene),])
write.csv(v2b, "chr6_missese.csv", row.names=F)


##chr8
c2 <- read.delim("/media/cailu/Seagate Backup Plus Drive/data_20230311/UKB_geno/bfilemaf0.1/c8flt.bim", header=FALSE)
g1<-g[g$chromosome_name==8,]
snp<-data.frame()
for (i in 1: nrow(g1)){
  snpx<-c2[c2$V4>g1[i,]$start_position &c2$V4<g1[i,]$end_position,] 
  if(nrow(snpx)>0){
    snpx$gene<-g1[i,]$gene
    snp<-bind_rows(snp, snpx)
  }else{snpx[1,]<-c(NA, NA, NA,NA, NA, NA)
  snpx$gene<-g1[i,]$gene
  snp<-bind_rows(snp, snpx)
  }
}
v2<-merge(g1, snp, by="gene")

write.csv(v2, "/media/cailu/Seagate Backup Plus Drive1/data_20230311/UKB_geno/bfilemaf0.1/variants_generegionmaf0.1/variants_c8.csv", row.names=F)

v8_ana <- read_excel("ann_c8.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
names(v8_ana)[c(2, 5)]<-c("dbSNP.annotation", "EUR_freq")
v8_ana<-v6_ana[1,]
v8_ana[1,]<-c(NA, NA, NA, NA,NA)
v2a<-merge(v2,  v8_ana, by.x="V2", by.y="variant", all.x=T)%>%
  filter(dbSNP.annotation=="missense"|dbSNP.annotation=="frameshift"|dbSNP.annotation=="nonsense")

v2b<-bind_rows(v2a, g1[!(g1$gene%in%v2a$gene),])
write.csv(v2b, "chr8_missese.csv", row.names=F)


##cX
c2 <- read.delim("/media/cailu/Seagate Backup Plus Drive/data_20230311/UKB_geno/bfilemaf0.1/cXflt.bim", header=FALSE)
g1<-g[g$chromosome_name=="X",]
snp<-data.frame()
for (i in 1: nrow(g1)){
  snpx<-c2[c2$V4>g1[i,]$start_position &c2$V4<g1[i,]$end_position,] 
  if(nrow(snpx)>0){
    snpx$gene<-g1[i,]$gene
    snp<-bind_rows(snp, snpx)
  }else{snpx[1,]<-c(NA, NA, NA,NA, NA, NA)
  snpx$gene<-g1[i,]$gene
  snp<-bind_rows(snp, snpx)
  }
}
v2<-merge(g1, snp, by="gene")

write.csv(v2, "/media/cailu/Seagate Backup Plus Drive1/data_20230311/UKB_geno/bfilemaf0.1/variants_generegionmaf0.1/variants_cX.csv", row.names=F)
vx_ana <- read_excel("ann_cX.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)
names(vx_ana)[c(2, 5)]<-c("dbSNP.annotation", "EUR_freq")

v2a<-merge(v2,  vx_ana, by.x="V2", by.y="variant", all.x=T)%>%
  filter(dbSNP.annotation=="missense"|dbSNP.annotation=="frameshift"|dbSNP.annotation=="nonsense")

v2b<-bind_rows(v2a, g1[!(g1$gene%in%v2a$gene),])
write.csv(v2b, "chrx_missese.csv", row.names=F)

##15
c2 <- read.delim("/media/cailu/Seagate Backup Plus Drive/data_20230311/UKB_geno/bfilemaf0.1/c15flt.bim", header=FALSE)
g1<-g[g$chromosome_name==15,]
snp<-data.frame()
for (i in 1: nrow(g1)){
  snpx<-c2[c2$V4>g1[i,]$start_position &c2$V4<g1[i,]$end_position,] 
  if(nrow(snpx)>0){
    snpx$gene<-g1[i,]$gene
    snp<-bind_rows(snp, snpx)
  }else{snpx[1,]<-c(NA, NA, NA,NA, NA, NA)
  snpx$gene<-g1[i,]$gene
  snp<-bind_rows(snp, snpx)
  }
}
v2<-merge(g1, snp, by="gene")

write.csv(v2, "/media/cailu/Seagate Backup Plus Drive1/data_20230311/UKB_geno/bfilemaf0.1/variants_generegionmaf0.1/variants_c15.csv", row.names=F)
v15_ana <- read_excel("ann_c15.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)
names(v15_ana)[c(2, 5)]<-c("dbSNP.annotation", "EUR_freq")
v2a<-merge(v2,  v15_ana, by.x="V2", by.y="variant", all.x=T)%>%
  filter(dbSNP.annotation=="missense"|dbSNP.annotation=="frameshift"|dbSNP.annotation=="nonsense")

v2b<-bind_rows(v2a, g1[!(g1$gene%in%v2a$gene),])
write.csv(v2b, "chr15_missese.csv", row.names=F)

##chr3
c2 <- read.delim("/media/cailu/Seagate Backup Plus Drive/data_20230311/UKB_geno/bfilemaf0.1/c3flt.bim", header=FALSE)
g1<-g[g$chromosome_name==3,]
snp<-data.frame()
for (i in 1: nrow(g1)){
  snpx<-c2[c2$V4>g1[i,]$start_position &c2$V4<g1[i,]$end_position,] 
  if(nrow(snpx)>0){
    snpx$gene<-g1[i,]$gene
    snp<-bind_rows(snp, snpx)
  }else{snpx[1,]<-c(NA, NA, NA,NA, NA, NA)
  snpx$gene<-g1[i,]$gene
  snp<-bind_rows(snp, snpx)
  }
}
v2<-merge(g1, snp, by="gene")

write.csv(v2, "/media/cailu/Seagate Backup Plus Drive1/data_20230311/UKB_geno/bfilemaf0.1/variants_generegionmaf0.1/variants_c3.csv", row.names=F)
v3_ana <- read_excel("ann_c3.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
names(v3_ana)[c(2, 5)]<-c("dbSNP.annotation", "EUR_freq")

v2a<-merge(v2,  v3_ana, by.x="V2", by.y="variant", all.x=T)%>%
  filter(dbSNP.annotation=="missense"|dbSNP.annotation=="frameshift"|dbSNP.annotation=="nonsense")

v2b<-bind_rows(v2a, g1[!(g1$gene%in%v2a$gene),])
write.csv(v2b, "chr3_missese.csv", row.names=F)

##
##chr9
c2 <- read.delim("/media/cailu/Seagate Backup Plus Drive/data_20230311/UKB_geno/bfilemaf0.1/c9flt.bim", header=FALSE)
g1<-g[g$chromosome_name==9,]
snp<-data.frame()
for (i in 1: nrow(g1)){
  snpx<-c2[c2$V4>g1[i,]$start_position &c2$V4<g1[i,]$end_position,] 
  if(nrow(snpx)>0){
    snpx$gene<-g1[i,]$gene
    snp<-bind_rows(snp, snpx)
  }else{snpx[1,]<-c(NA, NA, NA,NA, NA, NA)
  snpx$gene<-g1[i,]$gene
  snp<-bind_rows(snp, snpx)
  }
}
v2<-merge(g1, snp, by="gene")

write.csv(v2, "/media/cailu/Seagate Backup Plus Drive1/data_20230311/UKB_geno/bfilemaf0.1/variants_generegionmaf0.1/variants_c9.csv", row.names=F)
v9_ana <- read_excel("ann_c9.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
names(v9_ana)[c(2, 5)]<-c("dbSNP.annotation", "EUR_freq")

v2a<-merge(v2,  v9_ana, by.x="V2", by.y="variant", all.x=T)%>%
  filter(dbSNP.annotation=="missense"|dbSNP.annotation=="frameshift"|dbSNP.annotation=="nonsense")

v2b<-bind_rows(v2a, g1[!(g1$gene%in%v2a$gene),])
write.csv(v2b, "chr9_missese.csv", row.names=F)

##chr16
c2 <- read.delim("/media/cailu/Seagate Backup Plus Drive/data_20230311/UKB_geno/bfilemaf0.1/c16flt.bim", header=FALSE)
g1<-g[g$chromosome_name==16,]
snp<-data.frame()
for (i in 1: nrow(g1)){
  snpx<-c2[c2$V4>g1[i,]$start_position &c2$V4<g1[i,]$end_position,] 
  if(nrow(snpx)>0){
    snpx$gene<-g1[i,]$gene
    snp<-bind_rows(snp, snpx)
  }else{snpx[1,]<-c(NA, NA, NA,NA, NA, NA)
  snpx$gene<-g1[i,]$gene
  snp<-bind_rows(snp, snpx)
  }
}
v2<-merge(g1, snp, by="gene")

write.csv(v2, "/media/cailu/Seagate Backup Plus Drive1/data_20230311/UKB_geno/bfilemaf0.1/variants_generegionmaf0.1/variants_c16.csv", row.names=F)
v16_ana <- read_excel("ann_c16.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(dbSNP)%>%
  distinct(variant, .keep_all=T)
names(v16_ana)[c(2, 5)]<-c("dbSNP.annotation", "EUR_freq")

v2a<-merge(v2,  v16_ana, by.x="V2", by.y="variant", all.x=T)%>%
  filter(dbSNP.annotation=="missense"|dbSNP.annotation=="frameshift"|dbSNP.annotation=="nonsense")

v2b<-bind_rows(v2a, g1[!(g1$gene%in%v2a$gene),])
write.csv(v2b, "chr16_missese.csv", row.names=F)

##chr22
c2 <- read.delim("/media/cailu/Seagate Backup Plus Drive/data_20230311/UKB_geno/bfilemaf0.1/c22flt.bim", header=FALSE)
g1<-g[g$chromosome_name==22,]
snp<-data.frame()
for (i in 1: nrow(g1)){
  snpx<-c2[c2$V4>g1[i,]$start_position &c2$V4<g1[i,]$end_position,] 
  if(nrow(snpx)>0){
    snpx$gene<-g1[i,]$gene
    snp<-bind_rows(snp, snpx)
  }else{snpx[1,]<-c(NA, NA, NA,NA, NA, NA)
  snpx$gene<-g1[i,]$gene
  snp<-bind_rows(snp, snpx)
  }
}
v2<-merge(g1, snp, by="gene")

write.csv(v2, "/media/cailu/Seagate Backup Plus Drive1/data_20230311/UKB_geno/bfilemaf0.1/variants_generegionmaf0.1/variants_c22.csv", row.names=F)
v22_ana <- read_excel("ann_c22.xlsx")%>%
select(variant,`dbSNP`,Ref,Alt, `EUR\nfreq`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
names(v22_ana)[c(2, 5)]<-c("dbSNP.annotation", "EUR_freq")


v2a<-merge(v2,  v22_ana, by.x="V2", by.y="variant", all.x=T)%>%
  filter(dbSNP.annotation=="missense"|dbSNP.annotation=="frameshift"|dbSNP.annotation=="nonsense")

v2b<-bind_rows(v2a, g1[!(g1$gene%in%v2a$gene),])
write.csv(v2b, "chr22_missese.csv", row.names=F)





##chr14
c2 <- read.delim("/media/cailu/Seagate Backup Plus Drive/data_20230311/UKB_geno/bfilemaf0.1/c14flt.bim", header=FALSE)
g1<-g[g$chromosome_name==14,]
snp<-data.frame()
for (i in 1: nrow(g1)){
  snpx<-c2[c2$V4>g1[i,]$start_position &c2$V4<g1[i,]$end_position,] 
  if(nrow(snpx)>0){
    snpx$gene<-g1[i,]$gene
    snp<-bind_rows(snp, snpx)
  }else{snpx[1,]<-c(NA, NA, NA,NA, NA, NA)
  snpx$gene<-g1[i,]$gene
  snp<-bind_rows(snp, snpx)
  }
}
v2<-merge(g1, snp, by="gene")

write.csv(v2, "/media/cailu/Seagate Backup Plus Drive1/data_20230311/UKB_geno/bfilemaf0.1/variants_generegionmaf0.1/variants_c14.csv", row.names=F)

v14_ana <- read_excel("ann_c14.xlsx")%>%
  select(variant,`dbSNP`,Ref,Alt, `EUR`)%>%
  drop_na(`dbSNP`)%>%
  distinct(variant, .keep_all=T)
names(v14_ana)[c(2, 5)]<-c("dbSNP.annotation", "EUR_freq")
v2a<-merge(v2,  v14_ana, by.x="V2", by.y="variant", all.x=T)%>%
  filter(dbSNP.annotation=="missense"|dbSNP.annotation=="frameshift"|dbSNP.annotation=="nonsense")

v2b<-bind_rows(v2a, g1[!(g1$gene%in%v2a$gene),])
write.csv(v2b, "chr14_missese.csv", row.names=F)

##merge missense variants
c1 <- read.csv("chr1_missese.csv")
c2 <- read.csv("chr2_missese.csv")
c3 <- read.csv("chr3_missese.csv")
c5 <- read.csv("chr5_missese.csv")
c6 <- read.csv("chr6_missese.csv")
c7 <- read.csv("chr7_missese.csv")
c9 <- read.csv("chr9_missese.csv")
c10 <- read.csv("chr10_missese.csv")
c11 <- read.csv("chr11_missese.csv")
c12 <- read.csv("chr12_missese.csv")
c14 <- read.csv("chr14_missese.csv")
c15 <- read.csv("chr15_missese.csv")
c16 <- read.csv("chr16_missese.csv")
c17 <- read.csv("chr17_missese.csv")
c19 <- read.csv("chr19_missese.csv")
c22 <- read.csv("chr22_missese.csv")
cX <- read.csv("chrx_missese.csv")

missense_variants<-bind_rows(c1,c2,c3,c5,c6,c7,c9,c10,c11,c12,c14,c15,c16,c17,c19,c22)
missense_variants$chromosome_name<-as.character(missense_variants$chromosome_name)
missense_variants$V1<-as.character(missense_variants$V1)
missense_variants<-bind_rows(missense_variants, cX)


write.csv(missense_variants, "missense&frameshift_variants_with_maf0.01.csv", row.names = F)

##check 
c1 <- read.csv("/media/cailu/Seagate Backup Plus Drive/data_20230311/results3/variants_c1.csv")
c2 <- read.csv("/media/cailu/Seagate Backup Plus Drive/data_20230311/results3/variants_c2.csv")
c3 <- read.csv("/media/cailu/Seagate Backup Plus Drive/data_20230311/results3/variants_c3.csv")
c5 <- read.csv("/media/cailu/Seagate Backup Plus Drive/data_20230311/results3/variants_c5.csv")
c6 <- read.csv("/media/cailu/Seagate Backup Plus Drive/data_20230311/results3/variants_c6.csv")
c7 <- read.csv("/media/cailu/Seagate Backup Plus Drive/data_20230311/results3/variants_c7.csv")
#c8 <- read.csv("/media/cailu/Seagate Backup Plus Drive/data_20230311/results3/variants_c8.csv")
c9 <- read.csv("/media/cailu/Seagate Backup Plus Drive/data_20230311/results3/variants_c9.csv")
c10 <- read.csv("/media/cailu/Seagate Backup Plus Drive/data_20230311/results3/variants_c10.csv")
c11 <- read.csv("/media/cailu/Seagate Backup Plus Drive/data_20230311/results3/variants_c11.csv")
c12 <- read.csv("/media/cailu/Seagate Backup Plus Drive/data_20230311/results3/variants_c12.csv")
c14 <- read.csv("/media/cailu/Seagate Backup Plus Drive/data_20230311/results3/variants_c14.csv")
c15 <- read.csv("/media/cailu/Seagate Backup Plus Drive/data_20230311/results3/variants_c15.csv")
c16 <- read.csv("/media/cailu/Seagate Backup Plus Drive/data_20230311/results3/variants_c16.csv")
c17 <- read.csv("/media/cailu/Seagate Backup Plus Drive/data_20230311/results3/variants_c17.csv")
c19 <- read.csv("/media/cailu/Seagate Backup Plus Drive/data_20230311/results3/variants_c19.csv")
c22 <- read.csv("/media/cailu/Seagate Backup Plus Drive/data_20230311/results3/variants_c22.csv")
cX <- read.csv("/media/cailu/Seagate Backup Plus Drive/data_20230311/results3/variants_cX.csv")

var<-bind_rows(c1,c2,c3,c5,c6,c7,c9,c10,c11,c12,c14,c15,c16,c17,c19,c22,cX)



###20230727 find common variants number within the 425 target genes (#5278 variants with maf >0.05)
gene_regions <- read_excel("/media/cailu/Seagate Backup Plus Drive/data_20230311/results3/gene_regions.xlsx")
chr<-unique(gene_regions$chromosome_name)
variants<-data.frame()
for (i in chr){
  g<-gene_regions[gene_regions$chromosome_name==i,]
  for (j in 1:nrow(g)){
    m<-read.delim(paste0("/media/cailu/Seagate Backup Plus Drive/data_20230311/UKB_geno/plink/c",i,".bim"), header=FALSE)
    tmp <- m[m$V4 <= g[j,]$end_position & m$V4>= g[j,]$start_position,]
    tmp$V1<-as.character(tmp$V1)
    variants<-bind_rows(variants, tmp)
  }
}

##format upload file
for (i in c(1,2,3,5,6,7,8,9,10,11,12,14,15,16,17,19, 22, "X")){
  d<-read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data_20230311/UKB_geno/bfilemaf0.1/variants_generegionmaf0.1/variants_c",i,".csv"))
  write.table(d$V2, paste0("/media/cailu/Seagate Backup Plus Drive1/data_20230311/UKB_geno/bfilemaf0.1/variants_generegionmaf0.1/upload_c",i,".txt"), sep = "\t", row.names = F, quote = F, col.names=F)
}

##20240202_sort nonsenses and framshift variants


var_non_fram<-data.frame()
for (i in c(1,3,5,6,7,9,10,11,12,14,15,17,19)){
  v<-read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data_20230311/results3/variants_c",i,".csv"))
  
  
  an<-read_excel(paste0("/media/cailu/Seagate Backup Plus Drive1/data_20230311/results3/v",i,"_ana.xlsx"))%>%
  select(variant,`dbSNP\nfunc annot`,Ref,Alt, `EUR\nfreq`)%>%
    drop_na(`dbSNP\nfunc annot`)%>%
    distinct(variant, .keep_all=T)
  names(an)[c(2, 5)]<-c("dbSNP.annotation", "EUR_freq")
  
  v2a<-merge(v,  an, by.x="V2", by.y="variant", all.x=T)%>%
    filter(dbSNP.annotation=="nonsense"|dbSNP.annotation=="frameshift")
  if (nrow(v2a)>0){
  var_non_fram<-bind_rows(var_non_fram,v2a)}else{
  var_non_fram<-var_non_fram  
  }
}
write.csv(var_non_fram, "/media/cailu/Seagate Backup Plus Drive1/data_20230311/results3/variants_nonsense_frameshift.csv", row.names = F)


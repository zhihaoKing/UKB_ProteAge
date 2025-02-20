setwd("/Users/zhihaojin/Documents/neuroscience_PhD_in_FDU/PhD_project/bioinformatics/UKBiobank_version_2")
## http://www.sthda.com/english/wiki/cox-proportional-hazards-model
## https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html
rm(list=ls());gc()
library(tidyverse)
library(ggplot2)
library(cowplot)     
library(ggforestplot)
library(psych)
library(broom)        
library(gtsummary)
library(htmlTable)
library(tableone)   
library(survival)    
library(Hmisc)
library(rms)
library(splines) 
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(haven)
library(patchwork)
library(grid)
library(gridExtra)
library(splines)
library(ggeffects)
library(ggpubr)
library(Metrics)
library(Hmisc)
library(reshape2)
library(RColorBrewer)
library(ggsci)
library(ggview)
library(smplot2)
library(ggpmisc)
library(cowplot)
library(doMC)
library(ClusterGVis)
library(limma)
library(tidyverse)
library(readr)
library(magrittr)
library(dplyr)
library(BioAge)
library(haven)
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(grid)
library(gridExtra)
library(psych)
library(splines) # for natural spline model to calculate BA residuals
library(ggeffects)
library(ggpubr)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggview)
library(pROC)
library(ggvenn)
library(ggVennDiagram)
library(pheatmap)
library(mice)
library(enrichplot)


################################## differential expression protein
res_f <- read.csv("result/differential expression analysis/result_f_logfc_cutoff_0.csv")  %>% filter(!Sig == "None") %>% arrange(desc(logFC))
res_f$id_f <- mapIds(org.Hs.eg.db,res_f$X,'ENTREZID','SYMBOL') 
res_f <- na.omit(res_f)
gene_list_f <- res_f$logFC
names(gene_list_f) <- res_f$id_f
## GSEA KEGG
set.seed(123)
gsea_kegg_f <- gseKEGG(
  gene_list_f,    
  organism = "hsa",    
  pvalueCutoff = 0.9,
  pAdjustMethod = "BH"
)
write.csv(gsea_kegg_f@result,"result/GSEA/female_gsea_kegg_all.csv")

kegg_f_res <- gsea_kegg_f@result

########## visualization of female kegg result

gse <- list()

for(i in c(1:30)){
gse[[i]] <- gseaplot2(gsea_kegg_f,
           #title = gsea_kegg_f@result[i,"Description"],
           title = paste0(gsea_kegg_f@result[i,"Description"]," (Female Cohort)"),
           geneSetID = i,
           pvalue_table = T,
           base_size = 18,
           color = "black")
#pdf(paste0("result/GSEA/female/female_gsea_kegg_",i,".pdf"),w = 9,h = 7)
gse[[i]]
ggsave(paste0("result/GSEA/female/female_gsea_kegg_",i,".pdf"),gse[[i]],w = 9,h = 7)
#dev.off()
}


## GSEA GO
#gsea_go_f <- gseGO(
#  geneList = gene_list_f,
#  OrgDb = org.Hs.eg.db,
#  ont = "ALL",                
#  pvalueCutoff = 0.9,       
# verbose = TRUE,
#  keyType = "ENTREZID"       
#)
#write.csv(gsea_go_f@result,"result/GSEA/female_gsea_go_all.csv")
#go_f_res <- gsea_go_f@result
#gseaplot2(gsea_go_f, 
#          title = "Female Cohort GO", 
#          geneSetID = 1:6,
#          pvalue_table = T,
#          base_size = 25,
#          ggsci::pal_lancet()(6))
#ggview(w = 12,h = 7.5)
#ggsave("submit/figure/figure4/log2(1.05)/female_gsea_go_all.pdf",w = 12,h = 7.5)



res_m <- read.csv("result/differential expression analysis/result_m_logfc_cutoff_0.csv") %>% filter(!Sig == "None") %>% arrange(desc(logFC))
res_m$id_m <- mapIds(org.Hs.eg.db,res_m$X,'ENTREZID','SYMBOL') 
res_m <- na.omit(res_m)
gene_list_m <- res_m$logFC
names(gene_list_m) <- res_m$id_m
## GSEA KEGG
set.seed(456)
gsea_kegg_m <- gseKEGG(
  gene_list_m,    
  organism = "hsa",    
  pvalueCutoff = 0.9,
  pAdjustMethod = "BH"
)
write.csv(gsea_kegg_m@result,"result/GSEA/male_gsea_kegg_all.csv")
kegg_m_res <- gsea_kegg_m@result
## visualize all the pathway


gse_m <- list()
for(i in c(1:11)){
gse_m[[i]] <-  gseaplot2(gsea_kegg_m,
           #title = gsea_kegg_f@result[i,"Description"],
           title = paste0(gsea_kegg_m@result[i,"Description"]," (Male Cohort)"),
           geneSetID = i,
           pvalue_table = T,
           base_size = 18,
           color = "black")
gse_m[[i]]
ggsave(paste0("result/GSEA/male/male_gsea_kegg_",i,".pdf"),gse_m[[i]],w = 9,h = 7)
}

## GSEA GO
#gsea_go_m <- gseGO(
#  geneList = gene_list_m,
#  OrgDb = org.Hs.eg.db,
#  ont = "ALL",                
#  pvalueCutoff = 0.9,       
#  verbose = TRUE,
#  keyType = "ENTREZID"       
#)
#write.csv(gsea_go_m@result,"submit/figure/figure4/log2(1)/figure4_GSEA/male_gsea_go_all.csv")#

#gseaplot2(gsea_go_m, 
#          title = "Male Cohort GO", 
#          geneSetID = 1,
#          pvalue_table = T,
#          base_size = 25,
#          ggsci::pal_lancet()(1))
#ggview(w = 12,h = 7.5)#

#ggsave("submit/figure/figure4/log2(1.05)/male_gsea_go_all.pdf",w = 12,h = 7.5)




################################## coefficient
## female
coef_f <- read.csv("result/ProteAge/female_coef_alpha0.05.csv",row.names = 1) %>% dplyr::filter(!coefficient == 0) %>% arrange(desc(coefficient))

#coef_f1 <- coef_f[c(1:300,1538:1837),]
names_f <- mapIds(org.Hs.eg.db,coef_f$protein,'ENTREZID','SYMBOL') 
gene_list_f <-  coef_f$coefficient
names(gene_list_f) <- names_f
gsea_go_f <- gseGO(
  geneList = gene_list_f,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",                
  pvalueCutoff = 0.9,       
  verbose = TRUE,
  keyType = "ENTREZID"       
)
write.csv(gsea_go_f@result,"submission/figure/figure 5/female_gsea_go_all.csv")

method <- 'p.adjust'
plot_data <- gsea_go_f[order(gsea_go_f[,method]),] %>% head(.,n=6)
plot_data$Description <- stringr::str_wrap(plot_data$Description, width = 40)

go <- enrichplot::gseaplot2(gsea_go_f,
                              plot_data$ID,
                              title = "Female Cohort GSEA GO",
                              color = c("#FF0000","#BA55D3","#6495ED","#48D1CC","#FFA500"),
                              base_size = 20,
                              rel_heights = c(1.5, 0.5, 0.5),
                              subplots = 1:3, 
                              ES_geom = "line",
                              pvalue_table = F) 
  go[[1]] <- go[[1]]+
    theme(axis.text = element_text(color = 'black'),
          axis.text.y = element_text(size = 15),
          axis.title = element_text(color = 'black',size = 16),
          panel.border = element_rect(fill = NA),
          legend.text = element_text(size = 14),
          legend.key.spacing.y = unit(0.1,'cm'),
          legend.position = 'top')+
    guides(color=guide_legend(ncol = 2))
  go[[3]] <- go[[3]]+
    theme(axis.text.x = element_text(size = 18,color = 'black'),
          axis.text.y = element_text(color = 'black'),
          axis.title.y = element_text(color = 'black',size = 16))
go
ggview(w = 11,h = 8)

ggsave("1st submission NC/Figures JZH 20241226/Figure 5/Figure 5D.pdf",w = 11,h = 8)



## male
coef_m <- read.csv("result/ProteAge/male_coef_alpha1.csv",row.names = 1) %>% dplyr::filter(!coefficient == 0) %>% arrange(desc(coefficient))
#coef_m1 <- coef_m[c(1:300,505:804),]
names_m <- mapIds(org.Hs.eg.db,coef_m$protein,'ENTREZID','SYMBOL') 
gene_list_m <-  coef_m$coefficient
names(gene_list_m) <- names_m

gsea_go_m <- gseGO(
  geneList = gene_list_m,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",                
  pvalueCutoff = 0.9,       
  verbose = TRUE,
  keyType = "ENTREZID"       
)
write.csv(gsea_go_m@result,"submission/figure/figure 5/male_gsea_go_all.csv")


method <- 'p.adjust'
plot_data <- gsea_go_m[order(gsea_go_m[,method]),] %>% head(.,n=6)
plot_data$Description <- stringr::str_wrap(plot_data$Description, width = 40)

go <- enrichplot::gseaplot2(gsea_go_m,
                              plot_data$ID,
                              title = "Male Cohort GSEA GO",
                              color = c("#FF0000","#BA55D3","#6495ED","#48D1CC","#FFA500"),
                              base_size = 20,
                              rel_heights = c(1.5, 0.5, 0.5),
                              subplots = 1:3, 
                              ES_geom = "line",
                              pvalue_table = F) 
  go[[1]] <- go[[1]]+
    theme(axis.text = element_text(color = 'black'),
          axis.text.y = element_text(size = 15),
          axis.title = element_text(color = 'black',size = 16),
          panel.border = element_rect(fill = NA),
          legend.text = element_text(size = 14),
          legend.key.spacing.y = unit(0.1,'cm'),
          legend.position = 'top')+
    guides(color=guide_legend(ncol = 2))
  go[[3]] <- go[[3]]+
    theme(axis.text.x = element_text(size = 18,color = 'black'),
          axis.text.y = element_text(color = 'black'),
          axis.title.y = element_text(color = 'black',size = 16))
go
ggview(w = 11,h = 8)

ggsave("1st submission NC/Figures JZH 20241226/Figure 5/Figure 5E.pdf",w = 11,h = 8)






################################## calculate correlation with age
### load data
p <- read.csv("data/Proteomics_data.csv")
p$X <- NULL
f <- read.csv("data/whole_cohort_sex_age.csv") 
f$X <- NULL

d <- merge(f,p,by = "EID")
d1 <- d %>% filter(Sex == "Male")
d2 <- d %>% filter(Sex == "Female")


#save(d1,d2,vvisCluster,file = "result/select_protein/change_with_age/data.RData")
load("result/select_protein/change_with_age/data.RData")

########################### male
for(i in 4:ncol(d1)) {
  d1[,i][is.na(d1[,i])] <- mean(d1[,i], na.rm=TRUE)
}
y1 <- as.numeric(d1[,"Age"])
d11 <- d1[,-c(1:3)]
colnames1 <- colnames(d11)
cor_data_df1 <- data.frame(colnames1)
for (i in 1:1463){
test1 <- cor.test(as.numeric(d11[,i]),y1,type="pearson")
cor_data_df1[i,2] <- test1$estimate
cor_data_df1[i,3] <- test1$p.value
}
names(cor_data_df1) <- c("symbol","correlation","pvalue")
cor_data_df1$FDR <- p.adjust(cor_data_df1$pvalue,method = "BH")
#cor_data_df1 <- cor_data_df1 %>% filter(FDR < 0.05) %>% arrange(correlation)

f1 <- read.csv("submit/supplementary/supplementary_table/protein_change_with_age/male_cohort.csv") %>% filter(cluster == "1"|cluster == "4")
p1 <- merge(cor_data_df1,f1,by.x = "symbol",by.y = "gene") %>% arrange(desc(correlation))
p1$entrez_id <- mapIds(org.Hs.eg.db,p1$symbol,'ENTREZID','SYMBOL') 
p11 <- p1[,c(2,39)] %>% na.omit()

gene_list1 <- p11$correlation
names(gene_list1) <- p11$entrez_id

######################################### GSEA KEGG
gsea_kegg_1 <- gseKEGG(
  gene_list1,    
  organism = "hsa",    
  pvalueCutoff = 0.9,
  pAdjustMethod = "none")
write.csv(gsea_kegg_1@result,"submit/figure/figure6/figure6_GSEA/figure6_male_gsea_kegg_all_11_5.csv")





######################################### GSEA GO
gsea_go_1 <- gseGO(
  geneList = gene_list1,
  OrgDb = org.Hs.eg.db,
  ont = "BP",                
  pvalueCutoff = 0.9,       
  verbose = TRUE,
  keyType = "ENTREZID",
  pAdjustMethod = "none"       
)
write.csv(gsea_go_1@result,"submit/figure/figure6/figure6_GSEA/figure6_male_gsea_go_BP_11_5.csv")





########################### female
for(i in 4:ncol(d2)) {
  d2[,i][is.na(d2[,i])] <- mean(d2[,i], na.rm=TRUE)
}
y2 <- as.numeric(d2[,"Age"])
d22 <- d2[,-c(1:3)]
colnames2 <- colnames(d22)
cor_data_df2 <- data.frame(colnames2)
for (i in 1:1463){
test2 <- cor.test(as.numeric(d22[,i]),y2,type="pearson")
cor_data_df2[i,2] <- test2$estimate
cor_data_df2[i,3] <- test2$p.value
}
names(cor_data_df2) <- c("symbol","correlation","pvalue")
cor_data_df2$FDR <- p.adjust(cor_data_df2$pvalue,method = "BH")


f2 <- read.csv("submit/supplementary/supplementary_table/protein_change_with_age/female_cohort.csv") %>% filter(cluster == "1"|cluster == "5")
p2 <- merge(cor_data_df2,f2,by.x = "symbol",by.y = "gene") %>% arrange(desc(correlation))
p2$entrez_id <- mapIds(org.Hs.eg.db,p2$symbol,'ENTREZID','SYMBOL') 
p22 <- p2[,c(2,39)] %>% na.omit()

gene_list2 <- p22$correlation
names(gene_list2) <- p22$entrez_id

######################################### GSEA KEGG
gsea_kegg_2 <- gseKEGG(
  gene_list2,    
  organism = "hsa",    
  pvalueCutoff = 0.9,
  pAdjustMethod = "none"
)
write.csv(gsea_kegg_2@result,"submit/figure/figure6/figure6_GSEA/figure6_female_gsea_kegg_all_11_5.csv")
female_kegg <- gsea_kegg_2@result

i = 4
gseaplot2(gsea_kegg_2, 
          title = paste0(gsea_kegg_2@result[i,"Description"]," (Female Cohort)"), 
          geneSetID = i,
          pvalue_table = T,
          base_size = 18,
          color = "black")
ggview(w = 12,h = 7.5)
#ggsave("submit/figure/figure6/figure6_GSEA/female_gsea_kegg_Gastric cancer.pdf",w = 7.5,h = 7.5)
ggsave(paste0("submit/figure/figure6/figure6_GSEA/female_gsea_kegg_",gsea_kegg_2@result[i,"Description"],".pdf"),w = 12,h = 7.5)
######################################### GSEA GO
gsea_go_2 <- gseGO(
  geneList = gene_list2,
  OrgDb = org.Hs.eg.db,
  ont = "BP",                
  pvalueCutoff = 0.9,       
  verbose = TRUE,
  keyType = "ENTREZID",
  pAdjustMethod = "none"      
)
write.csv(gsea_go_2@result,"submit/figure/figure6/figure6_GSEA/figure6_female_gsea_go_BP_11_5.csv")

gseaplot2(gsea_go_2, 
          title = "Female Cohort GO", 
          geneSetID = 1:6,
          pvalue_table = T,
          base_size = 25,
          ggsci::pal_lancet()(6))
ggview(w = 12,h = 7.5)
ggsave("submit/figure/figure6/female_gsea_go_all.pdf",w = 12,h = 7.5)


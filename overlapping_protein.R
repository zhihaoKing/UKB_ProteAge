setwd("/Users/zhihaojin/Documents/neuroscience_PhD_in_FDU/PhD_project/bioinformatics/UKBiobank_version_2")
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
library(VennDiagram)
library(Vennerable)

##################### differential expression analysis result ##################### 
f1 <- read.csv("result/differential expression analysis/result_f_logfc_cutoff_0.csv") 
f1_up <- f1 %>% filter(Sig == "Up") %>% dplyr ::rename(protein = X)
f1_down <- f1 %>% filter(Sig == "Down") %>% dplyr ::rename(protein = X)

m1 <- read.csv("result/differential expression analysis/result_m_logfc_cutoff_0.csv")
m1_up <-  m1 %>% filter(Sig == "Up") %>% dplyr ::rename(protein = X)
m1_down <-  m1 %>% filter(Sig == "Down") %>% dplyr ::rename(protein = X)
rm(list = c("f1","m1"))




##################### protein changed during aging result #####################
f2 <- read.csv("result/proteins_changed_with_age/female_cohort.csv",row.names = 1) %>% dplyr ::rename(protein = gene)
f2_up <- f2 %>% filter(cluster == "2")
f2_down <- f2 %>% filter(cluster == "5")

f2_up$X <- NULL
f2_down$X <- NULL

m2 <- read.csv("result/proteins_changed_with_age/male_cohort.csv",row.names = 1) %>% dplyr ::rename(protein = gene)
m2_up <- m2 %>% filter(cluster == "2")
m2_down <- m2 %>% filter(cluster == "5")

m2_up$X <- NULL
m2_down$X <- NULL
rm(list = c("f2","m2"))

# female proteage #992224

# male proteage #3E4F94


########################### Elastic Net Regression coef result #####################
f3 <- read.csv("result/ProteAge/female_coef_alpha0.05.csv",row.names = 1)
#colnames(f3) <- c("protein","coef")
f3_up <- f3 %>% filter(coefficient > 0) %>% mutate(group = "positive")
f3_down <- f3 %>% filter(coefficient < 0) %>% mutate(group = "negative")

m3 <- read.csv("result/ProteAge/male_coef_alpha1.csv",row.names = 1)
#colnames(m3) <- c("protein","coef")
m3_up <- m3 %>% filter(coefficient > 0) %>% mutate(group = "positive")
m3_down <- m3 %>% filter(coefficient < 0) %>% mutate(group = "negative")
rm(list = c("f3","m3"))

####### top10 protein
f3_down_top10 <- arrange(f3_down,coefficient)[c(1:20),]
f3_up_top10 <- arrange(f3_up,desc(coefficient))[c(1:20),]
m3_down_top10 <- arrange(m3_down,coefficient)[c(1:20),]
m3_up_top10 <- arrange(m3_up,desc(coefficient))[c(1:20),]

d_f <- rbind(f3_up_top10,f3_down_top10) %>% arrange(desc(coefficient))
d_f$protein <- factor(d_f$protein,levels = (d_f$protein))
d_m <- rbind(m3_up_top10,m3_down_top10) %>% arrange(desc(coefficient))
d_m$protein <- factor(d_m$protein,levels = (d_m$protein))

ggplot(d_f,aes(x = coefficient, y = protein))+
geom_point(size = 6,alpha = 0.7,aes(color = group))+
theme_classic(base_size = 14) + 
geom_segment(aes(x = 0, xend = coefficient, y = protein, yend = protein),size = 1.2,alpha = 0.5)+
theme(legend.position = "none",
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 18,face = "bold"),
      plot.title = element_text(size = 30,hjust = 0.5,face = "bold")) +
labs(x = "Coef in Elastic Net Regression", y = "Protein",title = "Female Cohort") +
scale_color_manual(values = c(positive = "#E3625D",
                             negative = "#3E90BF"))
ggview(w = 5,h = 8.5)
ggsave("1st submission NC/Figures JZH 20241226/Figure 5/Figure 5B.pdf",w = 5,h = 8.5,dpi = 300)

ggplot(d_m,aes(x = coefficient, y = protein))+
geom_point(size = 6, alpha = 0.7,aes(color = group))+
theme_classic(base_size = 14) + 
geom_segment(aes(x = 0, xend = coefficient, y = protein, yend = protein),size = 1.2,alpha = 0.5)+
theme(legend.position="none",
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 18,face = "bold"),
      plot.title = element_text(size = 30,hjust = 0.5,face = "bold")) +
labs(x = "Coef in Elastic Net Regression", y = "Protein",title = "Male Cohort") +
scale_color_manual(values = c(positive = "#F0C284",
                             negative = "#58B6E9"))
ggview(w = 5,h = 8.5)
ggsave("1st submission NC/Figures JZH 20241226/Figure 5/Figure 5C.pdf",w = 5,h = 8.5,dpi = 300)

### all protein coef
#f3 <- rbind(f3_down,f3_up) %>% arrange(desc(coef))
#f3$protein <- factor(f3$protein,levels = (f3$protein))
#m3 <- rbind(m3_down,m3_up) %>% arrange(desc(coef))
#m3$protein <- factor(m3$protein,levels = (m3$protein))

#ggplot(f3,aes(x = coef, y = protein))+
#geom_point(size = 0.1,alpha = 0.7,aes(color = group))+
#theme_classic(base_size = 14) + 
#geom_segment(aes(x = 0, xend = coef, y = protein, yend = protein),size = 0.1,alpha = 0.5)+
#theme(legend.position = "none",
#      axis.text = element_text(size = 0.5),
#      axis.title = element_text(size = 2,face = "bold"),
#      plot.title = element_text(size = 6,hjust = 0.5,face = "bold")) +
#labs(x = "Coef in Elastic Net Regression", y = "Protein",title = "Female Cohort") +
#scale_color_manual(values = c(positive = "#E3625D",
#                             negative = "#3E90BF"))
#ggview(w = 5,h = 8.5)
#ggsave("submit/supplymentary/supplymentary_figure/figure2_female_model.pdf",w = 5,h = 20,dpi = 300)

#ggplot(m3,aes(x = coef, y = protein))+
#geom_point(size = 0.1, alpha = 0.7,aes(color = group))+
#theme_classic(base_size = 14) + 
#geom_segment(aes(x = 0, xend = coef, y = protein, yend = protein),size = 0.1,alpha = 0.5)+
#theme(legend.position="none",
#      axis.text = element_text(size = 0.5),
#      axis.title = element_text(size = 2,face = "bold"),
#      plot.title = element_text(size = 6,hjust = 0.5,face = "bold")) +
#labs(x = "Coef in Elastic Net Regression", y = "Protein",title = "Male Cohort") +
#scale_color_manual(values = c(positive = "#F0C284",
#                             negative = "#58B6E9"))
#ggview(w = 5,h = 8.5)
#ggsave("submit/supplymentary/supplymentary_figure/figure2_male_model.pdf",w = 5,h = 20,dpi = 300)

######### overlapping proteins
f_down <- Reduce(intersect,list(f1_down$protein,f2_down$protein,f3_down$protein))
write.csv(f_down,"submission/figure/figure 7/deaccelerating aging rate protein in female cohort.csv")
p_f_down <- list(`Decreased in fast ager` = f1_down$protein,
                 `Decreased during aging` = f2_down$protein,
                 `Negative coef in ProteAge model` = f3_down$protein)


f_up <- Reduce(intersect,list(f1_up$protein,f2_up$protein,f3_up$protein))
write.csv(f_up,"submission/figure/figure 7/accelerating aging rate protein in female cohort.csv")
p_f_up <- list(`Increased in fast ager` = f1_up$protein,
               `Increased during aging` = f2_up$protein,
               `Positive coef in ProteAge model` = f3_up$protein)

m_down <- Reduce(intersect,list(m1_down$protein,m2_down$protein,m3_down$protein))
write.csv(m_down,"submission/figure/figure 7/deaccelerating aging rate protein in male cohort.csv")
p_m_down <- list(`Decreased in fast ager` = m1_down$protein,
                 `Decreased during aging` = m2_down$protein,
                 `Negative coef in ProteAge model` = m3_down$protein)

m_up <- Reduce(intersect,list(m1_up$protein,m2_up$protein,m3_up$protein))
write.csv(m_up,"submission/figure/figure 7/accelerating aging rate protein in male cohort.csv")
p_m_up <- list(`Increased in fast ager` = m1_up$protein,
                 `Increased during aging` = m2_up$protein,
                 `Positive coef in ProteAge model` = m3_up$protein)

#################### heatmap 
f11_up <- f1_up %>% filter(protein %in% f_up) %>% arrange(desc(logFC))
#f11_up <- f11_up %>% filter (protein %in% c("CCL19","LHB","XG","NTPROBNP","TFF1","ANXA10",
#                                            "NPPB","REG3A","CXCL10","GFAP","IGFBP4","MLN","HAVCR1","CHI3L1","LTBP2"))

f11_down <- f1_down %>% filter(protein %in% f_down) %>% arrange(logFC)
#f11_down <- f11_down %>% filter (protein %in% c("PRL","CTSV","CST6","RAB37","IL5","IL18RAP","CA6","PNPT1","GH2"))

m11_up <- m1_up %>% filter(protein %in% m_up) %>% arrange(desc(logFC))
#m11_up <- m11_up %>% filter (protein %in% c("CDH15","CDHR2","GP2","TSPAN1","ADGRG1","CHI3L1",
#                                            "EDA2R","MLN","PPY","CCN5","GDF15","NEFL","MMP12","TFF1","CDCP1"))

m11_down <- m1_down %>% filter(protein %in% m_down) %>% arrange(desc(logFC))

p <- read.csv("data/protein_data.csv")
colnames(p)[-1] <- toupper(colnames(p)[-1])
f <- read.csv("data/protein_cohort_sex_age.csv") %>% dplyr ::rename(Sex = p31,
                                                                    Age = p21022)


d <- merge(f,p,by = "eid")
d1 <- d %>% filter(Sex == "Male")
d2 <- d %>% filter(Sex == "Female")
col_color = c("#3E4F94", "white", "#992224")
col_range = c(-1, 0, 1)
col_fun = circlize::colorRamp2(col_range, col_color)
##### male cohort 
d1$Sex <- NULL
rownames(d1) <- d1$EID
d1$EID <- NULL
g1 <- list()
age <- c(39:70)
for(i in seq(age)){
  g1[[i]] <- d1 %>% filter(Age == age[i]) %>% as.tibble()
}
mean1 <- list()
for(i in seq(age)){
  mean1[[i]] <- apply(g1[[i]][,-1],2,function(x){mean(na.omit(x))})
}
names(mean1) <- age
dat1 <- data.frame(mean1[[1]],mean1[[2]],mean1[[3]],mean1[[4]],mean1[[5]],mean1[[6]],mean1[[7]],
                  mean1[[8]],mean1[[9]],mean1[[10]],mean1[[11]],mean1[[12]],mean1[[13]],mean1[[14]],
                  mean1[[15]],mean1[[16]],mean1[[17]],mean1[[18]],mean1[[19]],mean1[[20]],mean1[[21]],
                  mean1[[22]],mean1[[23]],mean1[[24]],mean1[[25]],mean1[[26]],mean1[[27]],mean1[[28]],
                  mean1[[29]],mean1[[30]],mean1[[31]],mean1[[32]])
colnames(dat1) <- age
dat1[,1] <- NULL
dat1 <- na.omit(dat1)

m_up_heatmap <- dat1[m11_up$protein,]
#png("result/overlapping_protein/heatmap_male_up_overlapping_protein.png",w = 6,h = 9, units = 'in', res = 600)
pdf("1st submission NC/Supplemental figures JZH 20241226/supplemental figure 3 JZH 20241226/supplemental figure 3b JZH 20241226.pdf",w = 6.5,h = 11)
ComplexHeatmap::Heatmap(as.matrix(m_up_heatmap),
                        cluster_rows = TRUE,
                        cluster_columns = FALSE,
                        show_row_dend = F,
                        name = "Z-score",
                        col = col_fun,
                        column_title = "Expression of accelerating aging rate proteins in male cohort",
                        row_names_gp = gpar(fontsize = 2.5),
                        column_names_gp = gpar(fontsize = 7),
                        column_names_centered = T,
                        column_title_gp = gpar(fontsize = 12),
                        column_names_rot = 0)

dev.off()

m_down_heatmap <- dat1[m11_down$protein,]
#png("result/overlapping_protein/heatmap_male_down_overlapping_protein.png",w = 6,h = 11, units = 'in', res = 600)
pdf("1st submission NC/Supplemental figures JZH 20241226/supplemental figure 3 JZH 20241226/supplemental figure 3d JZH 20241226.pdf",w = 6.5,h = 9)
ComplexHeatmap::Heatmap(as.matrix(m_down_heatmap),
                        cluster_rows = T,
                        cluster_columns = FALSE,
                        show_row_dend = F,
                        name = "Z-score",
                        col = col_fun,
                        column_title = "Expression of deaccelerating aging rate proteins in male cohort",
                        row_names_gp = gpar(fontsize = 9),
                        column_names_gp = gpar(fontsize = 7),
                        column_names_centered = T,
                        column_title_gp = gpar(fontsize = 12),
                        column_names_rot = 0)
dev.off()


##### female cohort          
d2$Sex <- NULL
rownames(d2) <- d2$EID
d2$EID <- NULL
g2 <- list()
age <- c(40:70)
for(i in seq(age)){
  g2[[i]] <- d2 %>% filter(Age == age[i]) %>% as.tibble()
}
mean2 <- list()
for(i in seq(age)){
  mean2[[i]] <- apply(g2[[i]][,-1],2,function(x){mean(na.omit(x))})
}
names(mean2) <- age
dat2 <- data.frame(mean2[[1]],mean2[[2]],mean2[[3]],mean2[[4]],mean2[[5]],mean2[[6]],mean2[[7]],
                  mean2[[8]],mean2[[9]],mean2[[10]],mean2[[11]],mean2[[12]],mean2[[13]],mean2[[14]],
                  mean2[[15]],mean2[[16]],mean2[[17]],mean2[[18]],mean2[[19]],mean2[[20]],mean2[[21]],
                  mean2[[22]],mean2[[23]],mean2[[24]],mean2[[25]],mean2[[26]],mean2[[27]],mean2[[28]],
                  mean2[[29]],mean2[[30]],mean2[[31]])
colnames(dat2) <- age
dat2[,1] <- NULL
dat2 <- na.omit(dat2)


f_up_heatmap <- dat2[f11_up$protein,]
pdf("1st submission NC/Supplemental figures JZH 20241226/supplemental figure 3 JZH 20241226/supplemental figure 3a JZH 20241226.pdf",w = 6.5,h = 13)
ComplexHeatmap::Heatmap(as.matrix(f_up_heatmap),
                        cluster_rows = T,
                        cluster_columns = FALSE,
                        show_row_dend = F,
                        name = "Z-score",
                        col = col_fun,
                        column_title = "Expression of accelerating aging rate proteins in female cohort",
                        row_names_gp = gpar(fontsize = 1.8),
                        column_names_gp = gpar(fontsize = 7.5),
                        column_names_centered = T,
                        column_title_gp = gpar(fontsize = 12),
                        column_names_rot = 0)
dev.off()

f_down_heatmap <- dat2[f11_down$protein,]
pdf("1st submission NC/Supplemental figures JZH 20241226/supplemental figure 3 JZH 20241226/supplemental figure 3c JZH 20241226.pdf",w = 6.5,h = 9)
ComplexHeatmap::Heatmap(as.matrix(f_down_heatmap),
                        cluster_rows = T,
                        cluster_columns = FALSE,
                        show_row_dend = F,
                        name = "Z-score",
                        col = col_fun,
                        column_title = "Expression of deaccelerating aging rate proteins in female cohort",
                        row_names_gp = gpar(fontsize = 8),
                        column_names_gp = gpar(fontsize = 8),
                        column_names_centered = T,
                        column_title_gp = gpar(fontsize = 11),
                        column_names_rot = 0)
dev.off()


## venn plot
# female proteage #992224
# female KDM #E3625D
# female phenoage #F0C284
# male proteage #3E4F94
# male KDM #3E90BF
# male phenoage #58B6E9
#jpeg("result/select_protein/down_regulated_female_cohort.jpeg",width = 26, height = 30, units="cm", res = 300, quality=100)
#pdf("result/overlapping_protein/down_regulated_female_cohort.pdf",w = 10.5,h = 12.5)
#ggvenn(p_f_down,fill_color = c("#992224", "#E3625D","#F0C284"),text_size = 7.8,set_name_size = 10.5) +
#theme(plot.title = element_text(hjust = 0.5,vjust = 2,size = 55,face = "bold")) + labs(title = "Female Cohort")
#ggview(w = 10.5,h = 12.5)
#dev.off()

p_f_down_venn <- Venn(list(`Decreased in fast ager` = f1_down$protein,
                 `Decreased during aging` = f2_down$protein,
                 `Negative coef in ProteAge model` = f3_down$protein))
pdf("1st submission NC/Figures JZH 20241226/Figure 7/Figure 7C JZH 20241226.pdf",w = 10.5,h = 12.5)
plot(p_f_down_venn)
dev.off()

#jpeg("result/select_protein/up_regulated_female_cohort.jpeg",width = 26, height = 30, units="cm", res = 300, quality=100)
#pdf("result/overlapping_protein/up_regulated_female_cohort.pdf",w = 10.5,h = 12.5)
#ggvenn(p_f_up,fill_color = c("#CC625F", "#FB7E00","#FCC3B4"),text_size = 7.8,set_name_size = 10.5) +
#theme(plot.title = element_text(hjust = 0.5,vjust = 2,size = 55,face = "bold")) + labs(title = "Female Cohort")
#ggview(w = 10.5,h = 12.5)
#dev.off()

p_f_up_venn <- Venn(list(`Increased in fast ager` = f1_up$protein,
               `Increased during aging` = f2_up$protein,
               `Positive coef in ProteAge model` = f3_up$protein))
pdf("1st submission NC/Figures JZH 20241226/Figure 7/Figure 7A JZH 20241226.pdf",w = 10.5,h = 12.5)
plot(p_f_up_venn)
dev.off()

#jpeg("result/select_protein/down_regulated_male_cohort.jpeg",width = 26, height = 30, units="cm", res = 300, quality=100)
#pdf("result/overlapping_protein/down_regulated_male_cohort.pdf",w = 10.5,h = 12.5)
#ggvenn(p_m_down,fill_color = c("#748EBB", "#95D1D7","#CCE1DE"),text_size = 7.8,set_name_size = 10.5) +
#theme(plot.title = element_text(hjust = 0.5,vjust = 2,size = 55,face = "bold")) + labs(title = "Male Cohort")
#ggview(w = 10.5,h = 12.5)
#dev.off()

p_m_down_venn <- Venn(list(`Decreased in fast ager` = m1_down$protein,
                 `Decreased during aging` = m2_down$protein,
                 `Negative coef in ProteAge model` = m3_down$protein))
pdf("1st submission NC/Figures JZH 20241226/Figure 7/Figure 7G JZH 20241226.pdf",w = 10.5,h = 12.5)                 
plot(p_m_down_venn)
dev.off()

#jpeg("result/select_protein/up_regulated_male_cohort.jpeg",width = 26, height = 30, units="cm", res = 300, quality=100)
#pdf("result/overlapping_protein/up_regulated_male_cohort.pdf",w = 10.5,h = 12.5)
#ggvenn(p_m_up,fill_color = c("#3E4F94", "#3E90BF","#58B6E9"),text_size = 7.8,set_name_size = 10.5) +
#theme(plot.title = element_text(hjust = 0.5,vjust = 2,size = 55,face = "bold")) + labs(title = "Male Cohort")
#ggview(w = 10,h = 11.5)
#dev.off()

p_m_up_venn <- Venn(list(`Increased in fast ager` = m1_up$protein,
                 `Increased during aging` = m2_up$protein,
                 `Positive coef in ProteAge model` = m3_up$protein))
pdf("1st submission NC/Figures JZH 20241226/Figure 7/Figure 7E JZH 20241226.pdf",w = 10.5,h = 12.5)
plot(p_m_up_venn)
dev.off()
########### protein for protein-protein-interaction analysis
write.csv(f_down,"submission/figure/figure 7/PPI/female_down.csv")
write.csv(f_up,"submission/figure/figure 7/PPI/female_up.csv")
write.csv(m_down,"submission/figure/figure 7/PPI/male_down.csv")
write.csv(m_up,"submission/figure/figure 7/PPI/male_up.csv")

#### accelerating aging rate proteins in female and male
f_m_up <- intersect(f_up,m_up)
write.csv(f_m_up,"submission/figure/figure 7/overlapping accelerating aging rate proteins in female and male cohort.csv")
#pdf("result/overlapping_protein/up_regulated_male_female.pdf",w = 10.5,h = 10.5)
#ggvenn(list(`Female cohort` = f_up,`Male cohort` = m_up),fill_color = c("#992224", "#F0C284"),text_size = 10,set_name_size = 11) +
#      theme(plot.title = element_text(hjust = 0.5,vjust = 5.5,size = 55,face = "bold")) + labs(title = "Up-regulated proteins")
#ggview(w = 10.5,h = 12.5)
#dev.off()

f_m_up_venn <- Venn(list(`Female cohort` = f_up,
               `Male cohort` = m_up))
pdf("1st submission NC/Figures JZH 20241226/Figure 7/Figure 7I JZH 20241226.pdf",w = 10.5,h = 10.5)
plot(f_m_up_venn)
dev.off()

#### deaccelerating aging rate proteins in female and male
f_m_down <- intersect(f_down,m_down)
write.csv(f_m_down,"submission/figure/figure 7/overlapping deaccelerating aging rate proteins in female and male cohort.csv")
#pdf("submit/figure/figure7/log2(1)/down_regulated_male_female_11_4.pdf",w = 10.5,h = 10.5)
#ggvenn(list(`Female cohort` = f_down,`Male cohort` = m_down),fill_color = c("#3E4F94", "#58B6E9"),text_size = 10,set_name_size = 11) +
#      theme(plot.title = element_text(hjust = 0.5,vjust = 5.5,size = 55,face = "bold")) + labs(title = "Down-regulated proteins")
#ggview(w = 10.5,h = 12.5)
#dev.off()

f_m_down_venn <- Venn(list(`Female cohort` = f_down,
               `Male cohort` = m_down))
pdf("1st submission NC/Figures JZH 20241226/Figure 7/Figure 7K JZH 20241226.pdf",w = 10.5,h = 10.5)
plot(f_m_down_venn)
dev.off()

write.csv(f_m_up,"submission/figure/figure 7/PPI/female_male_up.csv")
write.csv(f_m_down,"submission/figure/figure 7/PPI/female_male_down.csv")

#### accelerating in female and deaccelerating in male
f_up_m_down <- intersect(f_up,m_down)
#write.csv(f_m_up,"submit/supplementary/supplementary_table/supplymentary table 11/overlapping accelerating aging rate proteins in female and male cohort.csv")
#pdf("submit/figure/figure7/log2(1)/up_regulated_male_female_11_4.pdf",w = 10.5,h = 10.5)
#ggvenn(list(`Female cohort` = f_up,`Male cohort` = m_up),fill_color = c("#992224", "#F0C284"),text_size = 10,set_name_size = 11) +
#      theme(plot.title = element_text(hjust = 0.5,vjust = 5.5,size = 55,face = "bold")) + labs(title = "Up-regulated proteins")
#ggview(w = 10.5,h = 12.5)
#dev.off()

f_up_m_down_venn <- Venn(list(`Accelerating aging rate proteins in female` = f_up,
               `Deaccelerating aging rate proteins in male` = m_down))
pdf("submission/figure/figure 7/accelerating in female and deaccelerating in male.pdf",w = 14.5,h = 10.5)
plot(f_up_m_down_venn)
dev.off()

#### deaccelerating in female and accelerating in male
f_down_m_up <- intersect(f_down,m_up)
#write.csv(f_m_up,"submit/supplementary/supplementary_table/supplymentary table 11/overlapping accelerating aging rate proteins in female and male cohort.csv")
#pdf("submit/figure/figure7/log2(1)/up_regulated_male_female_11_4.pdf",w = 10.5,h = 10.5)
#ggvenn(list(`Female cohort` = f_up,`Male cohort` = m_up),fill_color = c("#992224", "#F0C284"),text_size = 10,set_name_size = 11) +
#      theme(plot.title = element_text(hjust = 0.5,vjust = 5.5,size = 55,face = "bold")) + labs(title = "Up-regulated proteins")
#ggview(w = 10.5,h = 12.5)
#dev.off()

f_down_m_up_venn <- Venn(list(`Deaccelerating aging rate proteins in female` = f_down,
               `Accelerating aging rate proteins in male` = m_up))
pdf("submission/figure/figure 7/deaccelerating in female and accelerating in male.pdf",w = 14.5,h = 10.5)
plot(f_down_m_up_venn)
dev.off()


f_m_up0 <- c("NTPROBNP","TFF1","CXCL9","CDCP1","SPINK4","MLN","ANXA10","GDF15","CHI3L1","MMP12",
            "NEFL","HAVCR1","ADGRG1","EDA2R","CXCL14")

f_m_up_heatmap1 <- dat1[f_m_up,]
#f_m_up_heatmap1 <- dat1[f_m_up,]
pdf("1st submission NC/Supplemental figures JZH 20241226/supplemental figure 4 JZH 20241226/supplemental figure 4c JZH 20241226.pdf",w = 7,h = 10)
ComplexHeatmap::Heatmap(as.matrix(f_m_up_heatmap1),
                        cluster_rows = T,
                        cluster_columns = FALSE,
                        show_row_dend = F,
                        name = "Z-score",
                        col = col_fun,
                        column_title = "Accelerating aging rate proteins in both female and male cohort (expression in male cohort)",
                        row_names_gp = gpar(fontsize = 4),
                        column_names_gp = gpar(fontsize = 9),
                        column_names_centered = T,
                        column_title_gp = gpar(fontsize = 10),
                        column_names_rot = 0)
dev.off()

f_m_up_heatmap2 <- dat2[f_m_up,]
pdf("1st submission NC/Supplemental figures JZH 20241226/supplemental figure 4 JZH 20241226/supplemental figure 4a JZH 20241226.pdf",w = 7,h = 10)
ComplexHeatmap::Heatmap(as.matrix(f_m_up_heatmap2),
                        cluster_rows = T,
                        cluster_columns = FALSE,
                        show_row_dend = F,
                        name = "Z-score",
                        col = col_fun,
                        column_title = "Accelerating aging rate proteins in both female and male cohort (expression in female cohort)",
                        row_names_gp = gpar(fontsize = 4),
                        column_names_gp = gpar(fontsize = 9),
                        column_names_centered = T,
                        column_title_gp = gpar(fontsize = 10),
                        column_names_rot = 0)
dev.off()



f_m_down_heatmap1 <- dat1[f_m_down,]
pdf("1st submission NC/Supplemental figures JZH 20241226/supplemental figure 4 JZH 20241226/supplemental figure 4d JZH 20241226.pdf",w = 7,h = 8)
ComplexHeatmap::Heatmap(as.matrix(f_m_down_heatmap1),
                        cluster_rows = T,
                        cluster_columns = FALSE,
                        show_row_dend = F,
                        name = "Z-score",
                        col = col_fun,
                        column_title = "Deaccelerating aging rate proteins in both female and male cohort (expression in male cohort)",
                        row_names_gp = gpar(fontsize = 9),
                        column_names_gp = gpar(fontsize = 9),
                        column_names_centered = T,
                        column_title_gp = gpar(fontsize = 10),
                        column_names_rot = 0)
dev.off()

f_m_down_heatmap2 <- dat2[f_m_down,]
pdf("1st submission NC/Supplemental figures JZH 20241226/supplemental figure 4 JZH 20241226/supplemental figure 4b JZH 20241226.pdf",w = 7,h = 8)
ComplexHeatmap::Heatmap(as.matrix(f_m_down_heatmap2),
                        cluster_rows = T,
                        cluster_columns = FALSE,
                        show_row_dend = F,
                        name = "Z-score",
                        col = col_fun,
                        column_title = "Deaccelerating aging rate proteins in both female and male cohort (expression in female cohort)",
                        row_names_gp = gpar(fontsize = 9),
                        column_names_gp = gpar(fontsize = 9),
                        column_names_centered = T,
                        column_title_gp = gpar(fontsize = 10),
                        column_names_rot = 0)
dev.off()



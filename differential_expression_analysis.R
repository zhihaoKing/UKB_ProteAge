setwd("/Users/zhihaojin/Documents/neuroscience_PhD_in_FDU/PhD_project/bioinformatics/UKBiobank_version_2")

rm(list=ls());gc()
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
library(ggvenn)
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
library(ReactomePA)
library(Metrics)
####  function 
#get_BA_resids <- function(ukb_bioage,BA){
#  data = ukb_bioage %>% drop_na(BA)
#  # Basic model = regress on age alone
#  model <- parse(text = sprintf("lm(%s ~ ns(age, df = 3), data = data)", BA)) %>% eval()
#  model_predict <- ggpredict(model, terms = c("age"))
#  data[,"BA_res"] <- NA
#  data[!is.na(data[BA]),"BA_res"] <- resid(model)
#  return(residuals(model))
#}


go_analysis <- function(gene){
   go <- enrichGO(gene       = gene,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'ENTREZID',
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
    return(go)
}
kegg_analysis <- function(gene){
   kegg <- enrichKEGG(gene,
                      organism = "hsa",
                      keyType = "kegg",
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      #minGSSize = 10,
                      #maxGSSize = 500,
                      qvalueCutoff = 0.1)
   return(kegg)
}


########## upload data
p <- read.csv("data/protein_data.csv")
p$X <- NULL
colnames(p)[-1] <- toupper(colnames(p)[-1])
#ukb_bioage <- read.csv("result/bioage/UKB_BioAge_all_features.csv")
#ukb_bioage$X <- NULL

# female proteage #992224
# female KDM #E3625D
# female phenoage #F0C284
# male proteage #3E4F94
# male KDM #3E90BF
# male phenoage #58B6E9

########### female cohort
#proage_f <- read.csv("result/proteage/lasso_model/predicted_result/female_proteage_validation_cohort.csv")
#proage_f$X <- NULL
#ukb_bioage_f <- merge(proage_f,ukb_bioage,by = "eid") %>% dplyr::rename(age = age.x)
#ukb_bioage_f$age.y <- NULL

#for(BA in c("proteage","kdm","kdm_original","kdm_noglu","phenoage","phenoage_original","phenoage_noglu")){
#  BA_res <- paste0(BA, "_res")
#  ukb_bioage_f[,BA_res] = NA
#  ukb_bioage_f[!is.na(ukb_bioage_f[BA]),BA_res] <- get_BA_resids(ukb_bioage_f,BA)
#}
ukb_bioage_f <- read.csv("result/bioage/female_cohort.csv")
ukb_bioage_f$X <- NULL


mae_f <- mae(ukb_bioage_f$age,ukb_bioage_f$ProteAge) ## 2.237385
ukb_bioage_f$Ager <- ifelse(abs(ukb_bioage_f$ProteAge.res) < mae_f,"normal_ager","changed_ager")
ukb_bioage_f <- ukb_bioage_f %>% filter(! Ager == "normal_ager")
ukb_bioage_f$Ager <- ifelse(ukb_bioage_f$ProteAge.res > mae_f ,"faster_ager","slow_ager")
#ukb_bioage_f$ager <- ifelse(ukb_bioage_f$proteage - ukb_bioage_f$age < -mae_f,"slow_ager",ukb_bioage_f$ager)
table(ukb_bioage_f$Ager)



p_f <- merge(ukb_bioage_f[,c(1,10)],p,by = "eid")
p_f$Ager <- NULL
rownames(p_f) <- p_f$eid
p_f$eid <- NULL
ukb_bioage_f$Ager <- factor(ukb_bioage_f$Ager,levels = c("faster_ager","slow_ager"))
design <- model.matrix(~0 + Ager,data = ukb_bioage_f)
fit <- lmFit(t(p_f), design)
contr <- makeContrasts(Agerfaster_ager - Agerslow_ager, levels = design)
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
results_f <-  as.data.frame(top.table) %>% dplyr::rename(FDR = `adj.P.Val`,
                                                  Pvalue = `P.Value`)
range(results_f$logFC)
results_f$Sig = ifelse(results_f$FDR < 0.05 & 
                     abs(results_f$logFC) >= 0, 
                   ifelse(results_f$logFC > 0,'Up','Down'),'None')
table(results_f$Sig)

#results_f_up <- results_f %>% filter(Sig == "Up")
#results_f_down <- results_f %>% filter(Sig == "Down")
write.csv(results_f,"result/differential expression analysis/result_f_logfc_cutoff_0.csv")
#write.csv(results_f_up,"result/differential_expression_analysis/result_f_up.csv")
#write.csv(results_f_down,"result/differential_expression_analysis/result_f_down.csv")

##############
#proage_m <- read.csv("result/proteage/lasso_model/predicted_result/male_proteage_validation_cohort.csv")
#proage_m$X <- NULL
#ukb_bioage_m <- merge(proage_m,ukb_bioage,by = "eid") %>% dplyr::rename(age = age.x)
#ukb_bioage_m$age.y <- NULL

#for(BA in c("proteage","kdm","kdm_original","kdm_noglu","phenoage","phenoage_original","phenoage_noglu")){
#  BA_res <- paste0(BA, "_res")
#  ukb_bioage_m[,BA_res] = NA
#  ukb_bioage_m[!is.na(ukb_bioage_m[BA]),BA_res] <- get_BA_resids(ukb_bioage_m,BA)
#}
ukb_bioage_m <- read.csv("result/bioage/male_cohort.csv")
ukb_bioage_m$X <- NULL


mae_m <- mae(ukb_bioage_m$age,ukb_bioage_m$ProteAge) ## 2.445928
ukb_bioage_m$Ager <- ifelse(abs(ukb_bioage_m$ProteAge.res) < mae_m,"normal_ager","changed_ager")
ukb_bioage_m <- ukb_bioage_m %>% filter(! Ager == "normal_ager")
ukb_bioage_m$Ager <- ifelse(ukb_bioage_m$ProteAge.res > mae_m ,"faster_ager","slow_ager")
#ukb_bioage_f$ager <- ifelse(ukb_bioage_f$proteage - ukb_bioage_f$age < -mae_f,"slow_ager",ukb_bioage_f$ager)
table(ukb_bioage_m$Ager)



p_m <- merge(ukb_bioage_m[,c(1,10)],p,by = "eid")
p_m$Ager <- NULL
rownames(p_m) <- p_m$eid
p_m$eid <- NULL
ukb_bioage_m$Ager <- factor(ukb_bioage_m$Ager,levels = c("faster_ager","slow_ager"))
design <- model.matrix(~0 + Ager,data = ukb_bioage_m)
fit <- lmFit(t(p_m), design)
contr <- makeContrasts(Agerfaster_ager - Agerslow_ager, levels = design)
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
results_m <-  as.data.frame(top.table) %>% dplyr::rename(FDR = `adj.P.Val`,
                                                  Pvalue = `P.Value`)
range(results_m$logFC)
results_m$Sig = ifelse(results_m$FDR < 0.05 & 
                     abs(results_m$logFC) >= 0, 
                   ifelse(results_m$logFC > 0,'Up','Down'),'None')
table(results_m$Sig)


#results_m_up <- results_m %>% filter(Sig == "Up")
#results_m_down <- results_m %>% filter(Sig == "Down")
write.csv(results_m,"result/differential expression analysis/result_m_logfc_cutoff_0.csv")
#write.csv(results_m_up,"result/differential_expression_analysis/result_m_up.csv")
#write.csv(results_m_down,"result/differential_expression_analysis/result_m_down.csv")

results_f$lab <- rownames(results_f)
results_m$lab <- rownames(results_m)
save(results_f,results_m,file = "result/differential expression analysis/DEP_result.RData")
load("result/differential expression analysis/DEP_result.RData")

########## overlapping protein ###########
m_up <- results_m %>% filter(Sig == "Up") %>% arrange(desc(logFC))
m_down <- results_m %>% filter(Sig == "Down") %>% arrange(logFC)

f_up <- results_f %>% filter(Sig == "Up") %>% arrange(desc(logFC))
f_down <- results_f %>% filter(Sig == "Down") %>% arrange(logFC)

fe_ma_up <- intersect(f_up$lab,m_up$lab) #531
fe_up_speci <- setdiff(f_up$lab,m_up$lab) #613
ma_up_speci <- setdiff(m_up$lab,f_up$lab) #104
length(fe_ma_up) <- 613
length(ma_up_speci) <- 613
length(fe_up_speci) <- 613
up_protein <- tibble(fe_ma_up,fe_up_speci,ma_up_speci)

write.csv(up_protein,"result/differential expression analysis/overlapping upregulated DEPs in male and females.csv")
#png("result/differential expression analysis/up_regulated_male_female.png",w = 10.5,h = 10.5, units = 'in', res = 600)
pdf("1st submission NC/Figures JZH 20241226/Figure 4/Figure 4C.pdf",w = 10.5,h = 10.5)
ggvenn(list(`Female cohort` = f_up$lab,`Male cohort` = m_up$lab),fill_color = c("#992224", "#F0C284"),text_size = 10,set_name_size = 11) +
  theme(plot.title = element_text(hjust = 0.5,vjust = 5.5,size = 55,face = "bold")) + labs(title = "Up-regulated proteins")
ggview(w = 10.5,h = 10)
dev.off()

fe_ma_down <- intersect(f_down$lab,m_down$lab) # 23
fe_down_speci <- setdiff(f_down$lab,m_down$lab) # 48
ma_down_speci <- setdiff(m_down$lab,f_down$lab) # 93
length(fe_ma_down) <- 93
length(fe_down_speci) <- 93
length(ma_down_speci) <- 93
down_protein <- tibble(fe_ma_down,fe_down_speci,ma_down_speci)

write.csv(down_protein,"result/differential expression analysis/overlapping downregulated DEPs in male and females.csv")
#png("result/differential expression analysis/down_regulated_male_female.png",w = 10.5,h = 10.5, units = 'in', res = 600)
pdf("1st submission NC/Figures JZH 20241226/Figure 4/Figure 4D.pdf",w = 10.5,h = 10.5)
ggvenn(list(`Female cohort` = f_down$lab,`Male cohort` = m_down$lab),fill_color = c("#3E4F94", "#58B6E9"),text_size = 10,set_name_size = 11) +
  theme(plot.title = element_text(hjust = 0.5,vjust = 5.5,size = 55,face = "bold")) + labs(title = "Down-regulated proteins")
ggview(w = 10.5,h = 10)
dev.off()

########## kegg/go enrichment analysis ###########
bar_plot <- function(go,col,grou,tit,pos){
   d_f_up <- go@result
   d_f_up$logP <- -log10(d_f_up$p.adjust)
   d_f_up <- d_f_up %>% arrange(desc(Count))
   d_f_up$Description <- factor(d_f_up$Description, levels = d_f_up$Description)
   d_f_up <- d_f_up[c(1:10),]
   d_f_up$group <- grou
   p <- ggplot(d_f_up,aes(x = Count,y = Description)) +
                  geom_point(aes(size = logP*4),color = col) +
                  theme_classic(base_size = 14) +
                  theme(axis.text.y = element_text(size = 15,face = "bold"),
                        axis.text.x = element_text(size = 15,face = "bold"),
                        axis.title = element_text(size = 15,face = "bold"),
                        legend.title = element_text(face = "bold",size = 10),
                        plot.title = element_text(size = 16,face = "bold",hjust = 0.98),
                        legend.background = element_rect(linetype="solid",colour ="black"),
                        legend.key.size = unit(1, units = "pt"),
                        legend.position = pos) +
                  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
                  labs(y = "", x = "Count",
                       title = tit, 
                       size = "-logP.adj") +
                  scale_size("-logP.adj", range = c(5, 10))
   return(p)
}
## female male up regulated
g_fe_ma_up <- mapIds(org.Hs.eg.db, fe_ma_up, 'ENTREZID','SYMBOL') 

go_fe_ma_up <- go_analysis(g_fe_ma_up)
go_fe_ma_up@result$GeneRatio <- sapply(go_fe_ma_up@result$GeneRatio, function(x) eval(parse(text=x)))
go_fe_ma_up@result$GeneRatio <- round(go_fe_ma_up@result$GeneRatio,2)
write.csv(go_fe_ma_up@result %>% arrange(desc(Count)),"result/differential expression analysis/go analysis of up-regulated proteins in female and male.csv")
bar_plot(go_fe_ma_up,"#de946c","Up-regulated in Female and Male",
         "GO analysis of up-regulated proteins in female and male",c(0.8, 0.75))
ggview(w = 7.5,h = 6.4)
ggsave("1st submission NC/Figures JZH 20241226/Figure 4/Figure 4F.pdf",w = 7.5,h = 6.4,dpi = 300)

kegg_fe_ma_up <- kegg_analysis(g_fe_ma_up)
kegg_fe_ma_up@result$GeneRatio <- sapply(kegg_fe_ma_up@result$GeneRatio, function(x) eval(parse(text=x)))
kegg_fe_ma_up@result$GeneRatio <- round(kegg_fe_ma_up@result$GeneRatio,2)
write.csv(kegg_fe_ma_up@result %>% arrange(desc(Count)),"result/differential expression analysis/kegg_go/kegg analysis of up-regulated proteins in female and male.csv")
bar_plot(kegg_fe_ma_up,"#de946c","Up-regulated in Female and Male",
        "KEGG analysis of up-regulated proteins in female and male",c(0.75, 0.8))
ggview(w = 7.5,h = 6.4)
ggsave("1st submission NC/Figures JZH 20241226/Figure 4/Figure 4I.pdf",w = 7.5,h = 6.4,dpi = 300)

## female up regulated specific
g_fe_up_speci <- mapIds(org.Hs.eg.db, fe_up_speci, 'ENTREZID','SYMBOL') 
go_fe_up_speci <- go_analysis(g_fe_up_speci)
go_fe_up_speci@result$GeneRatio <- sapply(go_fe_up_speci@result$GeneRatio, function(x) eval(parse(text=x)))
go_fe_up_speci@result$GeneRatio <- round(go_fe_up_speci@result$GeneRatio,2)
write.csv(go_fe_up_speci@result %>% arrange(desc(Count)),"result/differential expression analysis/kegg_go/go analysis of up-regulated proteins specific in female.csv")
bar_plot(go_fe_up_speci,"#ce848d","Specificly Up-regulated in Female",
         "GO analysis of up-regulated proteins specific in female",c(0.75, 0.8))
ggview(w = 7.5,h = 6.4)
ggsave("1st submission NC/Figures JZH 20241226/Figure 4/Figure 4E.pdf",w = 7.5,h = 6.4,dpi = 300)

kegg_fe_up_speci <- kegg_analysis(g_fe_up_speci)
kegg_fe_up_speci@result$GeneRatio <- sapply(kegg_fe_up_speci@result$GeneRatio, function(x) eval(parse(text=x)))
kegg_fe_up_speci@result$GeneRatio <- round(kegg_fe_up_speci@result$GeneRatio,2)
write.csv(kegg_fe_up_speci@result %>% arrange(p.adjust),"result/differential expression analysis/kegg_go/kegg analysis of up-regulated proteins specific in female.csv")
bar_plot(kegg_fe_up_speci,"#ce848d","Specificly Up-regulated in Female",
         "KEGG analysis of up-regulated proteins specific in female",c(0.75, 0.8))
ggview(w = 7.5,h = 6.4)
ggsave("1st submission NC/Figures JZH 20241226/Figure 4/Figure 4H.pdf",w = 7.5,h = 6.4,dpi = 300)

## male up regulated specific
g_ma_up_speci <- mapIds(org.Hs.eg.db, ma_up_speci, 'ENTREZID','SYMBOL') 

go_ma_up_speci <- go_analysis(g_ma_up_speci)
go_ma_up_speci@result$GeneRatio <- sapply(go_ma_up_speci@result$GeneRatio, function(x) eval(parse(text=x)))
go_ma_up_speci@result$GeneRatio <- round(go_ma_up_speci@result$GeneRatio,2)

write.csv(go_ma_up_speci@result %>% arrange(desc(Count)),"result/differential expression analysis/kegg_go/go analysis of up-regulated proteins specific in male.csv")
bar_plot(go_ma_up_speci,"#fbd8b3","Specificly Up-regulated in Male",
         "GO analysis of up-regulated proteins specific in male",c(0.75, 0.8))
ggview(w = 7.5,h = 6.4)
ggsave("1st submission NC/Figures JZH 20241226/Figure 4/Figure 4G.pdf",w = 7.5,h = 6.4,dpi = 300)


kegg_ma_up_speci <- kegg_analysis(g_ma_up_speci)
kegg_ma_up_speci@result$GeneRatio <- sapply(kegg_ma_up_speci@result$GeneRatio, function(x) eval(parse(text=x)))
kegg_ma_up_speci@result$GeneRatio <- round(kegg_ma_up_speci@result$GeneRatio,2)
write.csv(kegg_ma_up_speci@result %>% arrange(desc(Count)),"result/differential expression analysis/kegg_go/kegg analysis of up-regulated proteins specific in male.csv")
bar_plot(kegg_ma_up_speci,"#fbd8b3","Specificly Up-regulated in Male",
         "KEGG analysis of up-regulated proteins specific in male",c(0.75, 0.8))
ggview(w = 7.5,h = 6.4)
ggsave("1st submission NC/Figures JZH 20241226/Figure 4/Figure 4J.pdf",w = 7.5,h = 6.4,dpi = 300)

##### select protein
f_m_up <- merge(m_up,f_up,by = "lab") 
f_m_up <- f_m_up %>% arrange(desc(logFC.x))

f_m_down <- merge(m_down,f_down,by = "lab")
f_m_down <- f_m_down %>% arrange(logFC.x)

overlap_down <- f_m_down[c(1:3),1]
overlap_up <- f_m_up[c(1:12),1]

female_up <- f_up[c(1:25),8]
male_up <- m_up[c(1:25),8]

female_down <- f_down[c(1:3),8]
male_down <- m_down[c(1:3),8]

female_lab <- c(overlap_up,overlap_down,female_up,female_down) %>% unique()
male_lab <- c(overlap_up,overlap_down,male_up,male_down) %>% unique()

##### volcano plot
sig <- as.data.frame(table(results_f$Sig))
total <- sig[1,2] + sig[3,2]
subti <- paste0("Down(",sig[1,2],") UP(",sig[3,2],") Total(",total,")")
results_f$lab <- rownames(results_f)
select_lab <- filter(results_f,!Sig == "None")
select_lab <- arrange(select_lab,logFC)
lab <- select_lab[c(1:5,815:834),8]
keyvals_f <- ifelse(
  results_f$Sig == "Down", '#3E90BF',
  ifelse(results_f$Sig == "Up", '#E3625D',
         'gray74'))
#keyvals[is.na(keyvals)] <- 'black'
names(keyvals_f)[keyvals_f == '#E3625D'] <- 'high'
names(keyvals_f)[keyvals_f == 'gray74'] <- 'mid'
names(keyvals_f)[keyvals_f == '#3E90BF'] <- 'low'
EnhancedVolcano(results_f,lab = rownames(results_f),
                x = 'logFC',
                y = 'FDR',
                xlab = bquote(~Log[2] ~ "FC"),
                ylab = bquote(~-Log[10] ~ italic(FDR)),
                title = "DEPs of Female Cohort",
                subtitle = subti, 
                #legendPosition = 'none',
                pCutoff = 0.05,
                FCcutoff = 0,
                labSize = 3.4,
                pointSize = 1.5,
                #labFace = 'italic',
                titleLabSize = 10,
                colCustom = keyvals_f,
                colAlpha = 1.2,
                drawConnectors = TRUE,
                selectLab = female_lab) +
  ggplot2::coord_cartesian(xlim=c(-1, 1)) + 
  theme_classic(base_size = 14) + 
  theme(legend.position="none",
        plot.title = element_text(size = 20,hjust = 0.5,face = "bold"))
ggview(w = 4.5,h = 6)
ggsave("1st submission NC/Figures JZH 20241226/Figure 4/Figure 4A JZH 20241228.pdf",w = 4.5,h = 6,dpi = 300)


sig <- as.data.frame(table(results_m$Sig))
total <- sig[1,2] + sig[3,2]
subti <- paste0("Down(",sig[1,2],") UP(",sig[3,2],") Total(",total,")")
results_m$lab <- rownames(results_m)
select_lab <- filter(results_m,!Sig == "None")
select_lab <- arrange(select_lab,logFC)
lab <- select_lab[c(1:2,463:477),8]
keyvals_m <- ifelse(
  results_m$Sig == "Down", '#58B6E9',
  ifelse(results_m$Sig == "Up", '#F0C284',
         'gray74'))
#keyvals[is.na(keyvals)] <- 'black'
names(keyvals_m)[keyvals_m == '#F0C284'] <- 'high'
names(keyvals_m)[keyvals_m == 'gray74'] <- 'mid'
names(keyvals_m)[keyvals_m == '#58B6E9'] <- 'low'
EnhancedVolcano(results_m,lab = rownames(results_m),
                x = 'logFC',
                y = 'FDR',
                xlab = bquote(~Log[2] ~ "FC"),
                ylab = bquote(~-Log[10] ~ italic(FDR)),
                title = "DEPs of Male Cohort",
                subtitle = subti, 
                #legendPosition = 'none',
                pCutoff = 0.05,
                FCcutoff = 0,
                labSize = 3.4,
                pointSize = 1.5,
                #labFace = 'italic',
                titleLabSize = 10,
                colCustom = keyvals_m,
                colAlpha = 1,
                drawConnectors = TRUE,
                selectLab = male_lab) +
  ggplot2::coord_cartesian(xlim=c(-1, 1)) + 
  theme_classic(base_size = 14) + 
  theme(legend.position="none",
        plot.title = element_text(size = 20,hjust = 0.5,face = "bold"))
ggview(w = 4.5,h = 6)
ggsave("1st submission NC/Figures JZH 20241226/Figure 4/Figure 4B JZH 20241228.pdf",w = 4.5,h = 6,dpi = 300)




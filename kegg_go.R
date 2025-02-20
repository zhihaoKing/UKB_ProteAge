setwd("/Users/zhihaojin/Documents/neuroscience_PhD_in_FDU/PhD_project/bioinformatics/UKBiobank_version_2")

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggview)
library(tidyverse)
library(cowplot)
rm(list=ls());gc()
### function
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



###### load data
female_up <- read.csv("result/kegg_go_PPI/PPI/log2(1)/female_up.csv")
g_f_up <- mapIds(org.Hs.eg.db,female_up$x,'ENTREZID','SYMBOL') 
go_f_up <- go_analysis(g_f_up)
kegg_f_up <- kegg_analysis(g_f_up)

female_down <- read.csv("result/kegg_go_PPI/PPI/log2(1)/female_down.csv")
g_f_down <- mapIds(org.Hs.eg.db,female_down$x,'ENTREZID','SYMBOL') 
go_f_down <- go_analysis(g_f_down)
kegg_f_down <-  kegg_analysis(g_f_down)

male_up <- read.csv("result/kegg_go_PPI/PPI/log2(1)/male_up.csv")
g_m_up <- mapIds(org.Hs.eg.db,male_up$x,'ENTREZID','SYMBOL') 
go_m_up <- go_analysis(g_m_up)
kegg_m_up <- kegg_analysis(g_m_up)

male_down <- read.csv("result/kegg_go_PPI/PPI/log2(1)/male_down.csv")
g_m_down <- mapIds(org.Hs.eg.db,male_down$x,'ENTREZID','SYMBOL') 
go_m_down <- go_analysis(g_m_down)
kegg_m_down <- kegg_analysis(g_m_down)

save(go_f_down,go_f_up,go_m_down,go_m_up,
     kegg_f_down,kegg_f_up,kegg_m_down,kegg_m_up,
     file = "result/kegg_go_PPI/kegg_go/log2(1)/enrichment_result_11_4.RData")

###########
rm(list=ls())
load("result/kegg_go_PPI/kegg_go/log2(1)/enrichment_result_11_4.RData")



#############################################    GO analysis    #############
####### female up
go_f_up@result$GeneRatio <- sapply(go_f_up@result$GeneRatio, function(x) eval(parse(text=x)))
go_f_up@result$GeneRatio <- round(go_f_up@result$GeneRatio,2)
d_f_up <- go_f_up@result
write.csv(d_f_up %>% arrange(desc(Count)),"submit/figure/figure7/log2(1)/figure7_GO_KEGG_ranked_count/go analysis of up-regulated proteins in female cohort 11 18.csv")
d_f_up$logP <- -log10(d_f_up$p.adjust)
d_f_up <- d_f_up %>% arrange(desc(Count))
d_f_up$Description <- factor(d_f_up$Description, levels = d_f_up$Description)
d_f_up <- d_f_up[c(1:10),]
d_f_up$group <- "Up-regulaed in Female Cohort"
ggplot(d_f_up,aes(x = Count,y = Description,size = logP)) +
                  geom_point(color = "#FA9474") +
                  theme_classic(base_size = 14) +
                  theme(axis.text.y = element_text(size = 11.2,face = "bold"),
                        axis.text.x = element_text(size = 10,face = "bold"),
                        axis.title = element_text(size = 12,face = "bold"),
                        legend.title = element_text(face = "bold",size = 10),
                        plot.title = element_text(size = 16,face = "bold",hjust = 0.98),
                        plot.subtitle = element_text(size = 16,face = "bold",hjust = -0.2),
                        legend.background = element_rect(linetype="solid",colour ="black"),
                        legend.key.size = unit(0.01, units = "cm"),
                        legend.position = c(0.7, 0.8)) +
                  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
                  labs(y = "", x = "Count",
                       title = "Accelerated aging proteins (female)", 
                       subtitle = "GO analysis", 
                       size = "-logP.adj") +
                       scale_size("-logP.adj", range = c(3, 6))
ggview(w = 5.2,h = 4.5)
ggsave("submit/figure/figure7/log2(1)/figure7_GO_KEGG_ranked_count/go analysis of up-regulated proteins in female cohort 11 18.pdf",w = 5.2,h = 4.5,dpi = 300)

####### male up
go_m_up@result$GeneRatio <- sapply(go_m_up@result$GeneRatio, function(x) eval(parse(text=x)))
go_m_up@result$GeneRatio <- round(go_m_up@result$GeneRatio,2)
d_m_up <- go_m_up@result
write.csv(d_m_up %>% arrange(desc(Count)),"submit/figure/figure7/log2(1)/figure7_GO_KEGG_ranked_count/go analysis of up-regulated proteins in male cohort 11 18.csv")
d_m_up$logP <- -log10(d_m_up$p.adjust)
d_m_up <- d_m_up %>% arrange(desc(Count))
d_m_up$Description <- factor(d_m_up$Description, levels = d_m_up$Description)
d_m_up <- d_m_up[c(1:10),]
d_m_up$group <- "Up-regulaed in Male Cohort"
ggplot(d_m_up,aes(x = Count,y = Description,size = logP)) +
                  geom_point(color = "#1F9BD0") +
                  theme_classic(base_size = 14) +
                  theme(axis.text.y = element_text(size = 11.2,face = "bold"),
                        axis.text.x = element_text(size = 10,face = "bold"),
                        axis.title = element_text(size = 12,face = "bold"),
                        legend.title = element_text(face = "bold",size = 10),
                        plot.title = element_text(size = 16,face = "bold",hjust = 0.98),
                        plot.subtitle = element_text(size = 16,face = "bold",hjust = -0.2),
                        legend.background = element_rect(linetype="solid",colour ="black"),
                        legend.key.size = unit(0.01, units = "cm"),
                        legend.position = c(0.79, 0.8)) +
                  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
                  labs(y = "", x = "Count",
                       title = "Accelerated aging proteins (male)", 
                       subtitle = "GO analysis",
                       size = "-logP.adj") +
                       scale_size("-logP.adj", range = c(3, 6))
ggview(w = 5.2,h = 4)
ggsave("submit/figure/figure7/log2(1)/figure7_GO_KEGG_ranked_count/go analysis of up-regulated proteins in male cohort 11 18.pdf",w = 5.2,h = 4.5,dpi = 300)



####### female down
#go_f_down@result$GeneRatio <- sapply(go_f_down@result$GeneRatio, function(x) eval(parse(text=x)))
#go_f_down@result$GeneRatio <- round(go_f_down@result$GeneRatio,2)
#d_f_down <- go_f_down@result
#write.csv(d_f_down %>% arrange(p.adjust),"submit/figure/figure7/log2(1)/figure7_GO_KEGG_ranked_by_pvalue/go analysis of down-regulated proteins in female cohort 11 4.csv")
#d_f_down$logP <- -log10(d_f_down$p.adjust)
#d_f_down <- d_f_down %>% arrange(p.adjust)
#d_f_down$Description <- factor(d_f_down$Description, levels = d_f_down$Description)
#d_f_down <- d_f_down[c(1:10),]
#d_f_down$group <- "Down-regulaed in Female Cohort"
#ggplot(d_f_down,aes(x = Count,y = Description,size = logP)) +
#                  geom_point(color = "#E7855F") +
#                  theme_classic(base_size = 14) +
#                 theme(axis.text.y = element_text(size = 11.2,face = "bold"),
#                        axis.text.x = element_text(size = 10,face = "bold"),
#                        axis.title = element_text(size = 12,face = "bold"),
#                        legend.title = element_text(face = "bold",size = 10),
#                        plot.title = element_text(size = 16,face = "bold",hjust = 1.1),
#                        plot.subtitle = element_text(size = 16,face = "bold",hjust = -0.5),
#                        legend.background = element_rect(linetype="solid",colour ="black"),
#                        legend.key.size = unit(0.01, units = "cm"),
#                        legend.position = c(0.8, 0.8)) +
#                  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
#                  labs(y = "", x = "Count",
#                       title = "Deaccelerated aging proteins (female)", 
#                       subtitle = "GO analysis",
#                       size = "-logP.adj") +
#                       scale_size("-logP.adj", range = c(3, 6))
#ggview(w = 5.2,h = 4.5)
#ggsave("submit/figure/figure7/log2(1)/figure7_GO_KEGG_ranked_by_pvalue/go analysis of down-regulated proteins in female cohort 11 4.pdf",w = 5.2,h = 4.5,dpi = 300)


#"#A4CCD4"
####### male down
go_m_down@result$GeneRatio <- sapply(go_m_down@result$GeneRatio, function(x) eval(parse(text=x)))
go_m_down@result$GeneRatio <- round(go_m_down@result$GeneRatio,2)
d_m_down <- go_m_down@result
write.csv(d_m_down %>% arrange(desc(Count)),"submit/figure/figure7/log2(1)/figure7_GO_KEGG_ranked_count/go analysis of down-regulated proteins in male cohort 11 18.csv")
d_m_down$logP <- -log10(d_m_down$p.adjust)
d_m_down <- d_m_down %>% arrange(desc(Count))
d_m_down$Description <- factor(d_m_down$Description, levels = d_m_down$Description)
d_m_down <- d_m_down[c(1:10),]
d_m_down$group <- "Down-regulaed in Male Cohort"
d_m_down$Description <- as.character(d_m_down$Description)
d_m_down[d_m_down == "phosphatidylinositol 3-kinase/protein kinase B signal transduction"] <- "PI3K/PKB signal transduction"
d_m_down[d_m_down == "regulation of phosphatidylinositol 3-kinase/protein kinase B signal transduction"] <- "regulation of PI3K/PKB signal transduction"
d_m_down[d_m_down == "positive regulation of phosphatidylinositol 3-kinase/protein kinase B signal transduction"] <- "positive regulation of PI3K/PKB signal transduction"
ggplot(d_m_down,aes(x = Count,y = Description,size = logP)) +
                  geom_point(color = "#A4CCD4") +
                  theme_classic(base_size = 14) +
                  theme(axis.text.y = element_text(size = 11.2,face = "bold"),
                        axis.text.x = element_text(size = 10,face = "bold"),
                        axis.title = element_text(size = 12,face = "bold"),
                        legend.title = element_text(face = "bold",size = 10),
                        plot.title = element_text(size = 16,face = "bold",hjust = 1.1),
                        plot.subtitle = element_text(size = 16,face = "bold",hjust = -0.5),
                        legend.background = element_rect(linetype="solid",colour ="black"),
                        legend.key.size = unit(0.01, units = "cm"),
                        legend.position = c(0.75, 0.83)) +
                  scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
                  labs(y = "", x = "Count",
                       title = "Deaccelerated aging proteins (male)", 
                       subtitle = "GO analysis",
                       size = "-logP.adj") +
                       scale_size("-logP.adj", range = c(3, 6))
ggview(w = 5.2,h = 4.5)
ggsave("submit/figure/figure7/log2(1)/figure7_GO_KEGG_ranked_count/go analysis of down-regulated proteins in male cohort 11 18.pdf",w = 5.2,h = 4.5,dpi = 300)


#############################################    KEGG analysis    #############
####### female up
kegg_f_up@result$GeneRatio <- sapply(kegg_f_up@result$GeneRatio, function(x) eval(parse(text=x)))
kegg_f_up@result$GeneRatio <- round(kegg_f_up@result$GeneRatio,2)
d_f_up <- kegg_f_up@result
write.csv(d_f_up %>% arrange(desc(Count)),"submit/figure/figure7/log2(1)/figure7_GO_KEGG_ranked_count/kegg analysis of up-regulated proteins in female cohort 11 18.csv")
d_f_up$logP <- -log10(d_f_up$p.adjust)
d_f_up <- d_f_up %>% arrange(desc(Count))
d_f_up$Description <- factor(d_f_up$Description, levels = d_f_up$Description)
d_f_up <- d_f_up[c(1:10),]
d_f_up$group <- "Up-regulaed in Female Cohort"
ggplot(d_f_up,aes(x = Count,y = Description,size = logP)) +
                  geom_point(color = "#FA9474") +
                  theme_classic(base_size = 14) +
                  theme(axis.text.y = element_text(size = 11.2,face = "bold"),
                        axis.text.x = element_text(size = 10,face = "bold"),
                        axis.title = element_text(size = 12,face = "bold"),
                        legend.title = element_text(face = "bold",size = 10),
                        plot.title = element_text(size = 16,face = "bold",hjust = 1.1),
                        plot.subtitle = element_text(size = 16,face = "bold",hjust = -0.7),
                        legend.background = element_rect(linetype="solid",colour ="black"),
                        legend.key.size = unit(0.01, units = "cm"),
                        legend.position = c(0.7, 0.8)) +
                  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
                  labs(y = "", x = "Count",
                       title = "Accelerated aging proteins (female)", 
                       subtitle = "KEGG analysis", 
                       size = "-logP.adj") +
                       scale_size("-logP.adj", range = c(3, 6))
ggview(w = 5.2,h = 4.5)
ggsave("submit/figure/figure7/log2(1)/figure7_GO_KEGG_ranked_count/kegg analysis of up-regulated proteins in female cohort 11 18.pdf",w = 5.2,h = 4.5,dpi = 300)

####### male up
kegg_m_up@result$GeneRatio <- sapply(kegg_m_up@result$GeneRatio, function(x) eval(parse(text=x)))
kegg_m_up@result$GeneRatio <- round(kegg_m_up@result$GeneRatio,2)
d_m_up <- kegg_m_up@result
write.csv(d_m_up %>% arrange(desc(Count)),"submit/figure/figure7/log2(1)/figure7_GO_KEGG_ranked_count/kegg analysis of up-regulated proteins in male cohort 11 18.csv")
d_m_up$logP <- -log10(d_m_up$p.adjust)
d_m_up <- d_m_up %>% arrange(desc(Count))
d_m_up$Description <- factor(d_m_up$Description, levels = d_m_up$Description)
d_m_up <- d_m_up[c(1:10),]
d_m_up$group <- "Up-regulaed in Male Cohort"
ggplot(d_m_up,aes(x = Count,y = Description,size = logP)) +
                  geom_point(color = "#1F9BD0") +
                  theme_classic(base_size = 14) +
                  theme(axis.text.y = element_text(size = 11.2,face = "bold"),
                        axis.text.x = element_text(size = 10,face = "bold"),
                        axis.title = element_text(size = 12,face = "bold"),
                        legend.title = element_text(face = "bold",size = 10),
                        plot.title = element_text(size = 16,face = "bold",hjust = 1.1),
                        plot.subtitle = element_text(size = 16,face = "bold",hjust = -0.7),
                        legend.background = element_rect(linetype="solid",colour ="black"),
                        legend.key.size = unit(0.01, units = "cm"),
                        legend.position = c(0.7, 0.8)) +
                  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
                  labs(y = "", x = "Count",
                       title = "Accelerated aging proteins (male)",
                       subtitle = "KEGG analysis", 
                       size = "-logP.adj") +
                       scale_size("-logP.adj", range = c(3, 6))
ggview(w = 5.2,h = 4.5)
ggsave("submit/figure/figure7/log2(1)/figure7_GO_KEGG_ranked_count/kegg analysis of up-regulated proteins in male cohort 11 18.pdf",w = 5.2,h = 4.5,dpi = 300)


####### female down
kegg_f_down@result$GeneRatio <- sapply(kegg_f_down@result$GeneRatio, function(x) eval(parse(text=x)))
kegg_f_down@result$GeneRatio <- round(kegg_f_down@result$GeneRatio,2)
d_f_down <- kegg_f_down@result
write.csv(d_f_down %>% arrange(desc(Count)),"submit/figure/figure7/log2(1)/figure7_GO_KEGG_ranked_count/kegg analysis of down-regulated proteins in female cohort 11 18.csv")
d_f_down$logP <- -log10(d_f_down$p.adjust)
d_f_down <- d_f_down %>% arrange(desc(Count))
d_f_down$Description <- factor(d_f_down$Description, levels = d_f_down$Description)
d_f_down <- d_f_down[c(1:10),]
d_f_down$group <- "Down-regulaed in Female Cohort"
ggplot(d_f_down,aes(x = Count,y = Description,size = logP)) +
                  geom_point(color = "#E7855F") +
                  theme_classic(base_size = 14) +
                  theme(axis.text.y = element_text(size = 11.2,face = "bold"),
                        axis.text.x = element_text(size = 10,face = "bold"),
                        axis.title = element_text(size = 12,face = "bold"),
                        legend.title = element_text(face = "bold",size = 10),
                        plot.title = element_text(size = 16,face = "bold",hjust = 1.1),
                        plot.subtitle = element_text(size = 16,face = "bold",hjust = -0.7),
                        legend.background = element_rect(linetype="solid",colour ="black"),
                        legend.key.size = unit(0.01, units = "cm"),
                        legend.position = c(0.7, 0.8)) +
                  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
                  labs(y = "", x = "Count",
                       title = "Deaccelerated aging proteins (female)",
                       subtitle = "KEGG analysis", 
                       size = "-logP.adj") +
                       scale_size("-logP.adj", range = c(3, 6))
ggview(w = 5.2,h = 4.5)
ggsave("submit/figure/figure7/log2(1)/figure7_GO_KEGG_ranked_count/kegg analysis of down-regulated proteins in female cohort 11 18.pdf",w = 5.2,h = 4.5,dpi = 300)

####### male down
kegg_m_down@result$GeneRatio <- sapply(kegg_m_down@result$GeneRatio, function(x) eval(parse(text=x)))
kegg_m_down@result$GeneRatio <- round(kegg_m_down@result$GeneRatio,2)
d_m_down <- kegg_m_down@result
write.csv(d_m_down,"submit/figure/figure7/log2(1)/figure7_GO_KEGG_ranked_count/kegg analysis of down-regulated proteins in male cohort 11 18.csv")
d_m_down$logP <- -log10(d_m_down$p.adjust)
#d_m_down <- d_m_down %>% arrange(desc(logP))
d_m_down <- d_m_down %>% arrange(desc(Count))
#d_m_down[1,3] <- "PI3K/PKB signal transduction"
#d_m_down[2,3] <- "regulation of PI3K/PKB signal transduction"
#d_m_down[3,3] <- "positive regulation of PI3K/PKB signal transduction"
d_m_down$Description <- factor(d_m_down$Description, levels = d_m_down$Description)
d_m_down <- d_m_down[c(1:10),]
d_m_down$group <- "Down-regulaed in Male Cohort"
ggplot(d_m_down,aes(x = Count,y = Description,size = logP)) +
                  geom_point(color = "#A4CCD4") +
                  theme_classic(base_size = 14) +
                  theme(axis.text.y = element_text(size = 11.2,face = "bold"),
                        axis.text.x = element_text(size = 10,face = "bold"),
                        axis.title = element_text(size = 12,face = "bold"),
                        legend.title = element_text(face = "bold",size = 10),
                        plot.title = element_text(size = 16,face = "bold",hjust = 1.1),
                        plot.subtitle = element_text(size = 16,face = "bold",hjust = -0.7),
                        legend.background = element_rect(linetype="solid",colour ="black"),
                        legend.key.size = unit(0.01, units = "cm"),
                        legend.position = c(0.8, 0.7)) +
                  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
                  labs(y = "", x = "Count",
                       title = "Deaccelerated aging proteins (male)",
                       subtitle = "KEGG analysis", 
                       size = "-logP.adj") +
                       scale_size("-logP.adj", range = c(3, 6))
ggview(w = 5.2,h = 4.5)
ggsave("submit/figure/figure7/log2(1)/figure7_GO_KEGG_ranked_count/kegg analysis of down-regulated proteins in male cohort 11 18.pdf",w = 5.2,h = 4.5,dpi = 300)






######### up-regulated in female and male
female_male_up <- read.csv("submission/figure/figure 7/PPI/female_male_up.csv")
female_male_up$X <- NULL
g_f_m_up <- mapIds(org.Hs.eg.db,female_male_up$x,'ENTREZID','SYMBOL') 
go_f_m_up <- enrichGO(gene       = g_f_m_up,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENTREZID',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
go_f_m_up@result$GeneRatio <- sapply(go_f_m_up@result$GeneRatio, function(x) eval(parse(text=x)))
go_f_m_up@result$GeneRatio <- round(go_f_m_up@result$GeneRatio,2)
d_f_m_up <- go_f_m_up@result
write.csv(d_f_m_up %>% arrange(desc(Count)),"submission/figure/figure 7/go analysis of up-regulated proteins in female and male cohort.csv")
d_f_m_up$logP <- -log10(d_f_m_up$p.adjust)
d_f_m_up <- d_f_m_up %>% arrange(desc(Count))
d_f_m_up$Description <- as.character(d_f_m_up$Description)
#d_f_m_up[9,3] <- "transmembrane RSTK signaling pathway"
d_f_m_up$Description <- factor(d_f_m_up$Description, levels = d_f_m_up$Description)
d_f_m_up <- d_f_m_up[c(1:10),]
d_f_m_up$group <- "Up-regulaed in Male and Female Cohort"
#ggplot(d_f_m_up,aes(x = Count,y = Description,size = logP)) +
#                  geom_point(color = "#E3625D") +
#                  theme_classic(base_size = 14) +
#                  theme(axis.text.y = element_text(size = 11.2,face = "bold"),
#                        axis.text.x = element_text(size = 10,face = "bold"),
#                        axis.title = element_text(size = 12,face = "bold"),
#                        legend.title = element_text(face = "bold",size = 10),
#                        plot.title = element_text(size = 14,face = "bold",hjust = 1),
#                        plot.subtitle = element_text(size = 16,face = "bold",hjust = -0.7),
#                        legend.background = element_rect(linetype="solid",colour ="black"),
#                        legend.key.size = unit(0.01, units = "cm"),
#                        legend.position = c(0.85, 0.7)) +
#                  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
#                  labs(y = "", x = "Count",
#                       title = "Accelerated aging proteins (female and male)",
#                       subtitle = "GO analysis", 
#                       size = "-logP.adj") +
#                       scale_size("-logP.adj", range = c(3, 6))
#ggview(w = 5.2,h = 4.5)
ggplot(d_f_m_up,aes(x = Description,y = Count)) +
                  geom_bar(fill = "#E3625D",stat = "identity") +
                  coord_flip() + 
                  theme_classic(base_size = 14) +
                  theme(axis.text.y = element_text(size = 13,face = "bold"),
                        axis.text.x = element_text(size = 13,face = "bold"),
                        axis.title = element_text(size = 15,face = "bold"),
                        plot.title = element_text(size = 18,face = "bold",hjust = 1.1)) +
                  scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
                  labs(y = "Count", x = "",
                       title = "GO analysis of accelerated aging proteins (female and male)",
                       size = "-logP.adj")
ggview(w = 8.5,h = 4)
ggsave("1st submission NC/Figures JZH 20241226/Figure 7/Figure 7M JZH 20241226.pdf",w = 8.5,h = 4,dpi = 300)

kegg_f_m_up <- enrichKEGG(g_f_m_up,
                        organism = "hsa",
                        keyType = "kegg",
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        #minGSSize = 10,
                        #maxGSSize = 500,
                        qvalueCutoff = 0.1)
kegg_f_m_up@result$GeneRatio <- sapply(kegg_f_m_up@result$GeneRatio, function(x) eval(parse(text=x)))
kegg_f_m_up@result$GeneRatio <- round(kegg_f_m_up@result$GeneRatio,2)
d_f_m_up <- kegg_f_m_up@result
#d_f_m_up <- read.csv("submission/figure/figure 7/kegg analysis of up-regulated proteins in female and male cohort.csv")
#rownames(d_f_m_up) <- d_f_m_up$X
#d_f_m_up$X <- NULL
write.csv(d_f_m_up %>% arrange(desc(Count)),"submission/figure/figure 7/kegg analysis of up-regulated proteins in female and male cohort.csv")
d_f_m_up$logP <- -log10(d_f_m_up$p.adjust)
d_f_m_up <- d_f_m_up %>% arrange(desc(Count))
#d_f_m_up$Description <- as.character(d_f_m_up$Description)
d_f_m_up$Description <- factor(d_f_m_up$Description, levels = d_f_m_up$Description)
d_f_m_up <- d_f_m_up[c(1:10),]
d_f_m_up$group <- "Up-regulaed in Male and Female Cohort"
#ggplot(d_f_m_up,aes(x = Count,y = Description,size = logP)) +
#                  geom_point(color = "#E3625D") +
#                  theme_classic(base_size = 14) +
#                  theme(axis.text.y = element_text(size = 11.2,face = "bold"),
#                        axis.text.x = element_text(size = 10,face = "bold"),
#                        axis.title = element_text(size = 12,face = "bold"),
#                        legend.title = element_text(face = "bold",size = 10),
#                        plot.title = element_text(size = 14,face = "bold",hjust = 1),
#                        plot.subtitle = element_text(size = 16,face = "bold",hjust = -0.7),
#                        legend.background = element_rect(linetype="solid",colour ="black"),
#                        legend.key.size = unit(0.01, units = "cm"),
#                        legend.position = c(0.8, 0.7)) +
#                  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
#                  labs(y = "", x = "Count",
#                       title = "Accelerated aging proteins (female and male)",
#                       subtitle = "KEGG analysis", 
#                       size = "-logP.adj") +
#                       scale_size("-logP.adj", range = c(3, 6))
#ggview(w = 5.2,h = 4.5)
ggplot(d_f_m_up,aes(x = Description,y = Count)) +
                  geom_bar(fill = "#E3625D",stat = "identity") +
                  coord_flip() + 
                  theme_classic(base_size = 14) +
                  theme(axis.text.y = element_text(size = 13,face = "bold"),
                        axis.text.x = element_text(size = 13,face = "bold"),
                        axis.title = element_text(size = 15,face = "bold"),
                        plot.title = element_text(size = 18,face = "bold",hjust = 1.1)) +
                  scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
                  labs(y = "Count", x = "",
                       title = "KEGG analysis of accelerated aging proteins (female and male)",
                       size = "-logP.adj")
#ggview(w = 8.5,h = 4)

ggsave("1st submission NC/Figures JZH 20241226/Figure 7/Figure 7N JZH 20241226.pdf",w = 8.5,h = 4,dpi = 300)


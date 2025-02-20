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
library(ezcox)
library(VennDiagram)
library(Vennerable)
### participant protein information
sex_age <- read.csv("data/protein_cohort_sex_age.csv") %>% dplyr::rename(Sex = p31,
                                                                        Age = p21022)
sex_age$X <- NULL
female <- sex_age %>% dplyr::filter(Sex == "Female")
male <- sex_age %>% dplyr::filter(Sex == "Male")

p <- read.csv("data/protein_data.csv")
colnames(p)[-1] <- toupper(colnames(p)[-1])
#p$X <- NULL

### participant survival information
f <- read.csv("data/protein_cohort_follow_up_information_whole.csv")
colnames(f) <- c("eid","Date_lost_to_follow_up","Date_of_attending_assessment_centre",
                 "Date_of_death","Age_when_attended_assessment_centre","Age_at_recruitment")

f <- 
  f %>%
  mutate(
    Date_lost_to_follow_up = ymd(Date_lost_to_follow_up),
    Date_of_attending_assessment_centre = ymd(Date_of_attending_assessment_centre),
    Date_of_death = ymd(Date_of_death)
    )
f <- 
  f %>%
  mutate(
    death_yrs = as.duration(Date_of_attending_assessment_centre %--% Date_of_death) / dyears(1)
  )

f <- 
  f %>%
  mutate(
    os_yrs = as.duration(Date_of_attending_assessment_centre %--% Date_lost_to_follow_up) / dyears(1)
  )

f$s <- ifelse(is.na(f$death_yrs) & is.na(f$os_yrs),NA,1)
f <- subset(f,s == 1)
f$s <- NULL
f$yrs <- ifelse(is.na(f$death_yrs),f$os_yrs,f$death_yrs)
f$status <- ifelse(is.na(f$death_yrs),0,1) 

### overlapping proteins
f_down <- read.csv("submission/figure/figure 7/PPI/female_down.csv") 
f_down_protein <- as.character(f_down$x)
f_up <- read.csv("submission/figure/figure 7/PPI/female_up.csv")
f_up_protein <- as.character(f_up$x)
m_down <- read.csv("submission/figure/figure 7/PPI/male_down.csv")
m_down_protein <- as.character(m_down$x)
m_up <- read.csv("submission/figure/figure 7/PPI/male_up.csv")
m_up_protein <- as.character(m_up$x)

f_m_up <- read.csv("submission/figure/figure 7/PPI/female_male_up.csv")
f_m_up_protein <- as.character(f_m_up$x)
f_m_down <- read.csv("submission/figure/figure 7/PPI/female_male_down.csv")
f_m_down_protein <- as.character(f_m_down$x)

f_m_protein <- c(f_m_up_protein,f_m_down_protein)

### https://shixiangwang.github.io/home/cn/post/2019-10-23-ezcox-for-batch-cox-models/

#### female up-regulated proteins 
#p_f_up <- p[,c("EID",f_up_protein)]
#p_f_up <- merge(female,p_f_up,by = "EID")
#p_f_up <- p_f_up[,-c(2,3)]

#p_f_up0 <- merge(f,p_f_up,by.x = "eid",by.y = "EID")
#p_f_up0 <- p_f_up0 %>% dplyr::rename(time = yrs)

#res_f_up = ezcox(p_f_up0, covariates = f_up_protein)
#res_f_up$Sig <- ifelse(res_f_up$p.value < 0.05,"yes","no")
#write.csv(res_f_up,"submit/supplementary/supplementary_table/protein_mortality/positive_aging_related_protein_female.csv")

#### female down-regulated proteins 
#p_f_down <- p[,c("EID",f_down_protein)]
#p_f_down <- merge(female,p_f_down,by = "EID")
#p_f_down <- p_f_down[,-c(2,3)]

#p_f_down0 <- merge(f,p_f_down,by.x = "eid",by.y = "EID")
#p_f_down0 <- p_f_down0 %>% dplyr::rename(time = yrs)

#res_f_down = ezcox(p_f_down0, covariates = f_down_protein)
#res_f_down$Sig <- ifelse(res_f_down$p.value < 0.05,"yes","no")
#write.csv(res_f_down,"submit/supplementary/supplementary_table/protein_mortality/negative_aging_related_protein_female.csv")

#### male up-regulated proteins 
#p_m_up <- p[,c("EID",m_up_protein)]
#p_m_up <- merge(male,p_m_up,by = "EID")
#p_m_up <- p_m_up[,-c(2,3)]

#p_m_up0 <- merge(f,p_m_up,by.x = "eid",by.y = "EID")
#p_m_up0 <- p_m_up0 %>% dplyr::rename(time = yrs)

#res_m_up = ezcox(p_m_up0, covariates = m_up_protein)
#res_m_up$Sig <- ifelse(res_m_up$p.value < 0.05,"yes","no")
#write.csv(res_m_up,"submit/supplementary/supplementary_table/protein_mortality/positive_aging_related_protein_male.csv")

#### male down-regulated proteins 
#p_m_down <- p[,c("EID",m_down_protein)]
#p_m_down <- merge(male,p_m_down,by = "EID")
#p_m_down <- p_m_down[,-c(2,3)]

#p_m_down0 <- merge(f,p_m_down,by.x = "eid",by.y = "EID")
#p_m_down0 <- p_m_down0 %>% dplyr::rename(time = yrs)

#res_m_down = ezcox(p_m_down0, covariates = m_down_protein)
#res_m_down$Sig <- ifelse(res_m_down$p.value < 0.05,"yes","no")
#write.csv(res_m_down,"submit/supplementary/supplementary_table/protein_mortality/negative_aging_related_protein_male.csv")

#### female up proteins
p_f <- p[,c("eid",f_m_protein)]
p_f <- merge(female,p_f,by = "eid")
p_f <- p_f[,-c(2,3)]

p_f0 <- merge(f,p_f,by = "eid")
p_f0 <- p_f0 %>% dplyr::rename(time = yrs)

res_f = ezcox(p_f0, covariates = f_m_protein)
res_f$Sig <- ifelse(res_f$p.value < 0.05,"yes","no")
write.csv(res_f,"submission/figure/figure 8/protein mortality prediction result in female cohort.csv")

#### male proteins
p_m <- p[,c("eid",f_m_protein)]
p_m <- merge(male,p_m,by = "eid")
p_m <- p_m[,-c(2,3)]

p_m0 <- merge(f,p_m,by = "eid")
p_m0 <- p_m0 %>% dplyr::rename(time = yrs)

res_m = ezcox(p_m0, covariates = f_m_protein)
res_m$Sig <- ifelse(res_m$p.value < 0.05,"yes","no")
write.csv(res_m,"submission/figure/figure 8/protein mortality prediction result in male cohort.csv")
##### 
m_sig <- res_m %>% filter(Sig == "yes")
f_sig <- res_f %>% filter(Sig == "yes")
mortality_associated_protein <- intersect(f_sig$Variable,m_sig$Variable)
protein_set <- c("TGFB1","PLAU","NRP1","MMP12","MMP1","ITGAV","IGF1R","GFAP","FGF5","ELN","DCN","CXCL9","CHI3L1","CHGA",
                 "CEACAM5","CDH2","CCL4","CCL11","ACTA2","TREM2",
                 "FASLG","KIT","PROK1","RET","CD1C","CR2","CTSV","KLK7")

######################### 
#res_f_sig <- res_f %>% dplyr::filter(Sig == "yes"| Variable %in% mortality_associated_protein) %>% dplyr::mutate(N = 2420) 
#res_f_sig <- res_f %>% dplyr::filter(Variable %in% mortality_associated_protein) %>% dplyr::mutate(N = 2420)
res_f_sig <- res_f %>% dplyr::filter(Sig == "yes"|Variable %in% protein_set) %>% dplyr::mutate(N = 2420) 
res_f_sig <- res_f_sig[,c(1,8:11,14)] %>% dplyr::rename(protein = Variable,
                                                              mean = HR,
                                                              lower = lower_95,
                                                              upper = upper_95,
                                                              p = p.value)
res_f_sig$HR2 <- paste0(res_f_sig$mean," (",res_f_sig$lower,",",res_f_sig$upper,")")
res_f_sig$trend <- ifelse(res_f_sig$protein %in% f_up_protein,"Accelerated aging protein","Deaccelerated aging protein")
res_f_sig <- res_f_sig %>% arrange(trend,desc(mean))
#res_f_up_sig <- res_f_up_sig[c(1:10),]

############ select certain proteins

#res_m_sig <- res_m %>% dplyr::filter(Sig == "yes") %>% arrange(desc(HR))
#select_protein <- c(mortality_associated_protein,"RET","ENPP5","CTSV","BCAN")
#res_m_sig <- res_m_sig %>% dplyr::filter( Variable %in% select_protein) %>% dplyr::mutate(N = 3476) ## 98 significant proteins in total
res_m_sig <- res_m %>% dplyr::filter(Variable %in% res_f_sig$protein) %>% dplyr::mutate(N = 3476)
#res_m_sig <- res_m %>% dplyr::filter(Sig == "yes"| Variable %in% mortality_associated_protein) %>% dplyr::mutate(N = 3476)
res_m_sig <- res_m_sig[,c(1,8:11,14)] %>% dplyr::rename(protein = Variable,
                                                              mean = HR,
                                                              lower = lower_95,
                                                              upper = upper_95,
                                                              p = p.value)
res_m_sig$HR2 <- paste0(res_m_sig$mean," (",res_m_sig$lower,",",res_m_sig$upper,")")
res_m_sig$trend <- ifelse(res_m_sig$protein %in% m_up_protein,"Accelerated aging protein","Deaccelerated aging protein")
res_m_sig <- res_m_sig %>% arrange(trend,desc(mean))

#top20_protein <- res_m_sig$protein[c(1:15,94:98)]
#select_protein <- c(top20_protein,hub_protein)

#res_m_sig <- res_m_sig %>% dplyr::filter(res_m_sig$protein %in% select_protein)
#res_m_sig <- res_m_sig %>% arrange(trend,desc(mean))

############ select all proteins
#res_m_sig <- res_m %>% dplyr::filter(Sig == "yes"| Variable %in% hub_protein) %>% dplyr::mutate(N = 2401)

#res_m_sig <- res_m_sig[,c(1,8:11,14)] %>% dplyr::rename(protein = Variable,
#                                                              mean = HR,
#                                                              lower = lower_95,
#                                                              upper = upper_95,
#                                                              p = p.value)
#res_m_sig$HR2 <- paste0(res_m_sig$mean," (",res_m_sig$lower,",",res_m_sig$upper,")")
#res_m_sig$trend <- ifelse(res_m_sig$protein %in% m_up_protein,"Accelerated aging protein","Deaccelerated aging protein")
#res_m_sig <- res_m_sig %>% arrange(trend,desc(mean))

######### female

dat_f <- rbind(c("Accelerated aging rate proteins",NA,NA,NA,"","",""),res_f_sig[c(1:35),],
             c("Deaccelerated aging rate proteins",NA,NA,NA,"","",""),res_f_sig[c(36:45),]) %>% as.data.frame

dat_f[,c(2:4)] <- lapply(dat_f[,c(2:4)],as.numeric)
dat_f[,c(5,6,7)] <- lapply(dat_f[,c(5,6,7)],as.character)
dat_f$` ` <- paste(rep(" ", 20), collapse = " ")
dat_f$protein <- ifelse(is.na(dat_f$mean),dat_f$protein,
                    paste0("  ", dat_f$protein))
dat_f <- dplyr::rename(dat_f,`HR (95% CI) `= HR2)
tm <- forestploter::forest_theme(base_size = 10,
                   ci_pch = 15,
                   ci_lty = 1)

g <- forestploter::forest(dat_f[,c(1,9,5:7)],
                          est = dat_f$mean,
                          lower = dat_f$lower,
                          upper = dat_f$upper,
                          sizes = 0.5,
                          ci_colum = 2,
                          xlim = c(0.5,1.75),
                          ticks_at = c(0.5, 1, 1.5),
                          ref_line = 1,
                          theme = tm)

g <- forestploter::edit_plot(g,row = c(1,37),gp = gpar(fontface = "bold"))

g <- forestploter::edit_plot(g,
               row = c(1,37),
               which = "background",
               gp = gpar(fill = "#CFCFCF"))

g <- forestploter::edit_plot(g,
               row = c(2:36,38:47),
               which = "background",
               gp = gpar(fill = "#FFFFFF"))

g <- forestploter::edit_plot(g,
               row = c(2:36),
               col = 2,
               which = "ci",
               gp = gpar(col = "#992224"))

g <- forestploter::edit_plot(g,
               row = c(38:47),
               col = 2,
               which = "ci",
               gp = gpar(col = "#3E4F94"))

#g <- forestploter::edit_plot(g,
#               row = c(2:25),
#               gp = gpar(col = "#992224",fill = "#FFFACD"))
#g <- forestploter::edit_plot(g,
#               row = c(27:28),
#               gp = gpar(col = "#3E4F94",fill = "#E0FFFF"))

g <- forestploter::insert_text(g,
                 text = "Female Cohort",
                 col = 1:6,
                 part = "header",
                 gp = gpar(fontface = "bold"))
pdf("1st submission NC/Figures JZH 20241226/Figure 8/Figure 8B JZH 20250106.pdf",width = 7, height = 15)   
g
dev.off()



######### male

dat_m <- rbind(c("Accelerated aging rate proteins",NA,NA,NA,"","",""),res_m_sig[c(1:35),],
               c("Deaccelerated aging rate proteins",NA,NA,NA,"","",""),res_m_sig[c(36:45),]) %>% as.data.frame

dat_m[,c(2:4)] <- lapply(dat_m[,c(2:4)],as.numeric)
dat_m[,c(5,6,7)] <- lapply(dat_m[,c(5,6,7)],as.character)
dat_m$` ` <- paste(rep(" ", 20), collapse = " ")
dat_m$protein <- ifelse(is.na(dat_m$mean),dat_m$protein,
                        paste0("   ", dat_m$protein))
dat_m <- dplyr::rename(dat_m,`HR (95% CI) `= HR2)
tm <- forestploter::forest_theme(base_size = 10,
                   ci_pch = 15,
                   ci_lty = 1)

g <- forestploter::forest(dat_m[,c(1,9,5:7)],
                          est = dat_m$mean,
                          lower = dat_m$lower,
                          upper = dat_m$upper,
                          sizes = 0.5,
                          ci_colum = 2,
                          xlim = c(0.5,1.75),
                          ticks_at = c(0.5, 1, 1.5),
                          ref_line = 1,
                          theme = tm)

g <- forestploter::edit_plot(g,row = c(1,37),gp = gpar(fontface = "bold"))

g <- forestploter::edit_plot(g,
               row = c(1,37),
               which = "background",
               gp = gpar(fill = "#CFCFCF"))

g <- forestploter::edit_plot(g,
               row = c(2:36,38:47),
               which = "background",
               gp = gpar(fill = "#FFFFFF"))

g <- forestploter::edit_plot(g,
               row = c(2:36),
               col = 2,
               which = "ci",
               gp = gpar(col = "#992224"))

g <- forestploter::edit_plot(g,
               row = c(38:47),
               col = 2,
               which = "ci",
               gp = gpar(col = "#3E4F94"))

#g <- forestploter::edit_plot(g,
#               row = c(2:25),
#               gp = gpar(col = "#992224",fill = "#FFFACD"))

#g <- forestploter::edit_plot(g,
#               row = c(27:28),
#               gp = gpar(col = "#3E4F94",fill = "#E0FFFF"))

g <- forestploter::insert_text(g,
                 text = "Male Cohort",
                 col = 1:6,
                 part = "header",
                 gp = gpar(fontface = "bold"))


pdf("1st submission NC/Figures JZH 20241226/Figure 8/Figure 8C JZH 20250106.pdf",width = 7, height = 15)   
g
dev.off()


## overlapping proteins
rm(list=ls());gc()
female <- read.csv("submission/figure/figure 8/protein mortality prediction result in female cohort.csv")
female_sig <- female %>% filter(Sig == "yes")

male <- read.csv("submission/figure/figure 8/protein mortality prediction result in male cohort.csv")
male_sig <- male %>% filter(Sig == "yes")

fe_ma_sig <- Venn(list(`Significant mortality-associated proteins in females` = female_sig$Variable,
               `Significant mortality-associated proteins in males` = male_sig$Variable))
pdf("1st submission NC/Figures JZH 20241226/Figure 8/Figure 8D JZH 20241226.pdf",w = 12.5,h = 8)
plot(fe_ma_sig)
dev.off()

mortality_associated_protein <- intersect(female_sig$Variable,male_sig$Variable)




































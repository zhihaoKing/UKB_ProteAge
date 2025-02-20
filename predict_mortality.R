setwd("/Users/zhihaojin/Documents/neuroscience_PhD_in_FDU/PhD_project/bioinformatics/UKBiobank_version_2")
## http://www.sthda.com/english/wiki/cox-proportional-hazards-model
## https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html

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
library(condsurv)
library(survival)
library(survminer)
library(BioAge)
library(haven)
library(readxl)
library(tidyr)
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
library(forestplot)
library(forestploter)
library(grid)
rm(list=ls());gc()
################################ Predict Mortality #################################

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


#ukb_bioage <- read.csv("result/bioage/UKB_BioAge_all_features.csv")

#ukb_bioage$X <- NULL 

#p_f <- read.csv("result/proteage/lasso_model/predicted_result/female_proteage_whole_cohort.csv")
#p_f$X <- NULL  
#p_m <- read.csv("result/proteage/lasso_model/predicted_result/male_proteage_whole_cohort.csv")
#p_m$X <- NULL  
#p <- rbind(p_m,p_f)
#ukb_bioage_whole <- merge(p,ukb_bioage,by = "eid") %>% dplyr::rename(age = age.x)
#ukb_bioage_whole$age.y <- NULL

#ukb_bioage_whole <- ukb_bioage_whole %>% dplyr::select(eid,age,gender,proteage,kdm_original,phenoage_original) %>% dplyr::rename(KDM = kdm_original,
                                                                                                                  # PhenoAge = phenoage_original)



#get_BA_resids <- function(BA){
#  data = ukb_bioage_whole %>% drop_na(BA)
#  # Basic model = regress on age alone
#  model <- parse(text = sprintf("lm(%s ~ ns(age, df = 3), data = data)", BA)) %>% eval()
#  model_predict <- ggpredict(model, terms = c("age"))
#  data[,"BA_res"] <- NA
#  data[!is.na(data[BA]),"BA_res"] <- resid(model)
#  return(residuals(model))
#}

#for(BA in c("proteage","KDM","PhenoAge")){
#  BA_res <- paste0(BA, "_res")
#  ukb_bioage_whole[,BA_res] = NA
#  ukb_bioage_whole[!is.na(ukb_bioage_whole[BA]),BA_res] <- get_BA_resids(BA)
#}


#ukb_bioage_whole <- ukb_bioage_whole %>% dplyr::rename(ProteAge = proteage,
#                                                `Chronological Age` = age,
#                                                Sex = gender,
#                                                `ProteAge Res` = proteage_res,
#                                                `KDM Res` = KDM_res,
#                                                `PhenoAge Res` = PhenoAge_res)  ## 35408 obs 37 variables  

## Survival analysis
#f0 <- merge(ukb_bioage_whole,f,by = "eid")

#proage_res_sd <- f0$`ProteAge Res` %>% na.omit() %>% sd()
#f0$`ProteAge Res/SD` <- as.numeric(f0$`ProteAge Res` / proage_res_sd)

#kdm_res_sd <- f0$`KDM Res` %>% na.omit() %>% sd()
#f0$`KDM Res/SD` <- as.numeric(f0$`KDM Res` / kdm_res_sd)

#pheno_res_sd <- f0$`PhenoAge Res` %>% na.omit() %>% sd()
#f0$`PhenoAge Res/SD` <- as.numeric(f0$`PhenoAge Res` / pheno_res_sd)


################################# female #################################
#########   40-70   #########
female <- read.csv("result/bioage/female_cohort_whole.csv",row.names = 1)
f1 <- merge(female,f, by = "eid") %>% dplyr::rename(`Chronological Age` = age)

proage_res_sd <- f1$`ProteAge.res` %>% na.omit() %>% sd()
f1$`ProteAge Res/SD` <- as.numeric(f1$`ProteAge.res` / proage_res_sd)

kdm_res_sd <- f1$`KDM.res` %>% na.omit() %>% sd()
f1$`KDM Res/SD` <- as.numeric(f1$`KDM.res` / kdm_res_sd)

pheno_res_sd <- f1$`PhenoAge.res` %>% na.omit() %>% sd()
f1$`PhenoAge Res/SD` <- as.numeric(f1$`PhenoAge.res` / pheno_res_sd)


cox1 <- coxph(Surv(yrs, status) ~ `Chronological Age` + `ProteAge Res/SD` + `KDM Res/SD` + `PhenoAge Res/SD`,
             data = f1)
tbl_regression(cox1,exp = TRUE)
c1 <- summary(cox1)
colnames(c1[["coefficients"]])
colnames(c1[["conf.int"]])

dat1 = as.data.frame(round(c1[["conf.int"]][, c(1, 3, 4)], 2))
dat1 = tibble::rownames_to_column(dat1, var = "Trait")
colnames(dat1)[2:4] = c("HR", "lower", "upper")

### https://cran.r-project.org/web/packages/forestplot/vignettes/forestplot.html
### https://cran.r-project.org/web/packages/forestploter/vignettes/forestploter-intro.html
### https://cloud.tencent.com/developer/article/1980956
### https://blog.csdn.net/dege857/article/details/127859291
dat1$HR2 <- paste0(dat1[, 2], " (", dat1[, 3], "-", dat1[, 4], ")")
dat1$p <- round(c1[["coefficients"]][,5],3)
dat1$Trait <- c("Chronological Age","ProteAge Res/SD","KDM Res/SD","PhenoAge Res/SD")
dat1 <- dplyr::rename(dat1,mean = HR)
dat1 <- as.tibble(dat1)
dat1$N <- 1510
#write.csv(dat1,"result/mortality_prediction/female/Female_cohort_aged_40_70.csv")
#jpeg("result/mortality_prediction/female/Female_cohort_aged_40_70.jpeg",width = 20, height = 7, units="cm", res = 300, quality=100)
#dat1 |>
#  forestplot::forestplot(labeltext = c(Trait, HR2, p,N),
             #clip = c(0.1, 2.5),
             #xlog = TRUE,
#             zero = 1,
#             lty.ci = "solid",
#             title = "HR per 1-SD increase of BA Residual of female cohort (aged 40-70)",
#             graph.pos = 2,
#             boxsize = 0.3) |>
#  fp_add_lines() |>
#  fp_set_style(box = "royalblue",
#               line = "darkblue",
#               summary = "royalblue") |> 
#  fp_add_header(Trait = c("", "Characters"),
#                HR2 = c("", "HR"),
#                p = c("","P-Value"),
#                N = c("","N")) |>
#  fp_set_zebra_style("#EFEFEF")
#dev.off()
#########   40-50   #########
f11 <- dplyr::filter(f1,`Chronological Age` >= 40 & `Chronological Age` < 50)
cox11 <- coxph(Surv(yrs, status) ~ `Chronological Age` + `ProteAge Res/SD` + `KDM Res/SD` + `PhenoAge Res/SD`,
             data = f11)
tbl_regression(cox11,exp = TRUE)
c11 <- summary(cox11)
colnames(c11[["coefficients"]])
colnames(c11[["conf.int"]])

dat11 = as.data.frame(round(c11[["conf.int"]][, c(1, 3, 4)], 2))
dat11 = tibble::rownames_to_column(dat11, var = "Trait")
colnames(dat11)[2:4] = c("HR", "lower", "upper")
dat11$HR2 <- paste0(dat11[, 2], " (", dat11[, 3], "-", dat11[, 4], ")")
dat11$p <- round(c11[["coefficients"]][,5],3)
dat11$Trait <- c("Chronological Age","ProteAge Res/SD","KDM Res/SD","PhenoAge Res/SD")
dat11 <- dplyr::rename(dat11,mean = HR)
dat11 <- as.tibble(dat11)
dat11$N <- 110
#write.csv(dat11,"result/mortality_prediction/female/Female_cohort_aged_40_50.csv")
#jpeg("result/mortality_prediction/female/Female_cohort_aged_40_50.jpeg",width = 20, height = 7, units="cm", res = 300, quality=100)
#dat11 |>
#  forestplot::forestplot(labeltext = c(Trait, HR2, p,N),
             #clip = c(0.1, 2.5),
             #xlog = TRUE,
#             zero = 1,
#             lty.ci = "solid",
#             title = "HR per 1-SD increase of BA Residual of female cohort (aged 40-50)",
#             graph.pos = 2,
#             boxsize = 0.3) |>
#  fp_add_lines() |>
#  fp_set_style(box = "royalblue",
#               line = "darkblue",
#               summary = "royalblue") |> 
#  fp_add_header(Trait = c("", "Characters"),
#                HR2 = c("", "HR"),
#                p = c("","P-Value"),
#                N = c("","N")) |>
#  fp_set_zebra_style("#EFEFEF")
#dev.off()
#########   50-60   #########
f12 <- dplyr::filter(f1,`Chronological Age` >= 50 & `Chronological Age` < 60)
cox12 <- coxph(Surv(yrs, status) ~ `Chronological Age` + `ProteAge Res/SD` + `KDM Res/SD` + `PhenoAge Res/SD`,
             data = f12)
tbl_regression(cox12,exp = TRUE)
c12 <- summary(cox12)
colnames(c12[["coefficients"]])
colnames(c12[["conf.int"]])

dat12 = as.data.frame(round(c12[["conf.int"]][, c(1, 3, 4)], 2))
dat12 = tibble::rownames_to_column(dat12, var = "Trait")
colnames(dat12)[2:4] = c("HR", "lower", "upper")
dat12$HR2 <- paste0(dat12[, 2], " (", dat12[, 3], "-", dat12[, 4], ")")
dat12$p <- round(c12[["coefficients"]][,5],3)
dat12$Trait <- c("Chronological Age","ProteAge Res/SD","KDM Res/SD","PhenoAge Res/SD")
dat12 <- dplyr::rename(dat12,mean = HR)
dat12 <- as.tibble(dat12)
dat12$N <- 332
#write.csv(dat12,"result/mortality_prediction/female/Female_cohort_aged_50_60.csv")
#jpeg("result/mortality_prediction/female/Female_cohort_aged_50_60.jpeg",width = 20, height = 7, units="cm", res = 300, quality=100)
#dat12 |>
#  forestplot::forestplot(labeltext = c(Trait, HR2, p,N),
             #clip = c(0.1, 2.5),
             #xlog = TRUE,
#             zero = 1,
#             lty.ci = "solid",
#             title = "HR per 1-SD increase of BA Residual of female cohort (aged 50-60)",
#             graph.pos = 2,
#             boxsize = 0.3) |>
#  fp_add_lines() |>
#  fp_set_style(box = "royalblue",
#               line = "darkblue",
#               summary = "royalblue") |> 
#  fp_add_header(Trait = c("", "Characters"),
#                HR2 = c("", "HR"),
#                p = c("","P-Value"),
#                N = c("","N")) |>
#  fp_set_zebra_style("#EFEFEF")
#dev.off()
#########   60-70   #########
f13 <- dplyr::filter(f1,`Chronological Age` >= 60 & `Chronological Age` <= 70)
cox13 <- coxph(Surv(yrs, status) ~ `Chronological Age` + `ProteAge Res/SD` + `KDM Res/SD` + `PhenoAge Res/SD`,
             data = f13)
tbl_regression(cox13,exp = TRUE)
c13 <- summary(cox13)
colnames(c13[["coefficients"]])
colnames(c13[["conf.int"]])

dat13 = as.data.frame(round(c13[["conf.int"]][, c(1, 3, 4)], 2))
dat13 = tibble::rownames_to_column(dat13, var = "Trait")
colnames(dat13)[2:4] = c("HR", "lower", "upper")
dat13$HR2 <- paste0(dat13[, 2], " (", dat13[, 3], "-", dat13[, 4], ")")
dat13$p <- round(c13[["coefficients"]][,5],3)
dat13$Trait <- c("Chronological Age","ProteAge Res/SD","KDM Res/SD","PhenoAge Res/SD")
dat13 <- dplyr::rename(dat13,mean = HR)
dat13 <- as.tibble(dat13)
dat13$N <- 1068
#write.csv(dat13,"result/mortality_prediction/female/Female_cohort_aged_60_70.csv")
#jpeg("result/mortality_prediction/female/Female_cohort_aged_60_70.jpeg",width = 20, height = 7, units="cm", res = 300, quality=100)
#dat13 |>
#  forestplot::forestplot(labeltext = c(Trait, HR2, p,N),
             #clip = c(0.1, 2.5),
             #xlog = TRUE,
#             zero = 1,
#             lty.ci = "solid",
#             title = "HR per 1-SD increase of BA Residual of female cohort (aged 60-70)",
#             graph.pos = 2,
#             boxsize = 0.3) |>
#  fp_add_lines() |>
#  fp_set_style(box = "royalblue",
#               line = "darkblue",
#               summary = "royalblue") |> 
#  fp_add_header(Trait = c("", "Characters"),
#                HR2 = c("", "HR"),
#                p = c("","P-Value"),
#                N = c("","N")) |>
#  fp_set_zebra_style("#EFEFEF")
#dev.off()

#################   put it together
dat_f <- rbind(c("40-70 years old",NA,NA,NA,"","",""),dat1,
             c("40-50 years old",NA,NA,NA,"","",""),dat11,
             c("50-60 years old",NA,NA,NA,"","",""),dat12,
             c("60-70 years old",NA,NA,NA,"","",""),dat13) %>% as.data.frame
#dat[is.na(dat)] <- ""
dat_f[,c(2:4)] <- lapply(dat_f[,c(2:4)],as.numeric)
dat_f[,c(5,6,7)] <- lapply(dat_f[,c(5,6,7)],as.character)
dat_f$group <- "Female"
write.csv(dat_f,"result/aging_clock_mortality/female_cohort.csv")

#jpeg("result/mortality_prediction/female/Female_cohort_aged.jpeg",width = 20, height = 20, units="cm", res = 300, quality=100)
#dat_f |>
#  forestplot::forestplot(labeltext = c(Trait, HR2, p,N),
#             #clip = c(0.5, 2.5),
#             xlog = TRUE,
#             zero = 1,
#             lty.ci = "solid",
#             title = "HR per 1-SD increase of BA Residual (Female Cohort)",
#             graph.pos = 2,
#             boxsize = 0.3,
#             is.summary = c(T,F,F,F,F,T,F,F,F,F,T,F,F,F,F,T,F,F,F,F)) |>
#  fp_add_lines() |>
#  fp_set_style(box = "royalblue",
##               line = "darkblue",
#              summary = "royalblue") |> 
# fp_add_header(Trait = c("", "Characters"),
#                HR2 = c("", "HR"),
#                p = c("","P-Value"),
#                N = c("","N")) |>
#  fp_set_zebra_style("#EFEFEF")
#dev.off()  
  
  
################################# male #################################
male <- read.csv("result/bioage/male_cohort_whole.csv",row.names = 1)
f2 <- merge(male,f, by = "eid") %>% dplyr::rename(`Chronological Age` = age)

proage_res_sd <- f2$`ProteAge.res` %>% na.omit() %>% sd()
f2$`ProteAge Res/SD` <- as.numeric(f2$`ProteAge.res` / proage_res_sd)

kdm_res_sd <- f2$`KDM.res` %>% na.omit() %>% sd()
f2$`KDM Res/SD` <- as.numeric(f2$`KDM.res` / kdm_res_sd)

pheno_res_sd <- f2$`PhenoAge.res` %>% na.omit() %>% sd()
f2$`PhenoAge Res/SD` <- as.numeric(f2$`PhenoAge.res` / pheno_res_sd)


cox2 <- coxph(Surv(yrs, status) ~ `Chronological Age` + `ProteAge Res/SD` + `KDM Res/SD` + `PhenoAge Res/SD`,
             data = f2)
tbl_regression(cox2,exp = TRUE)
c2 <- summary(cox2)
colnames(c2[["coefficients"]])
colnames(c2[["conf.int"]])

dat2 = as.data.frame(round(c2[["conf.int"]][, c(1, 3, 4)], 2))
dat2 = tibble::rownames_to_column(dat2, var = "Trait")
colnames(dat2)[2:4] = c("HR", "lower", "upper")

###https://cran.r-project.org/web/packages/forestplot/vignettes/forestplot.html
dat2$HR2 <- paste0(dat2[, 2], " (", dat2[, 3], "-", dat2[, 4], ")")
dat2$p <- round(c2[["coefficients"]][,5],3)
dat2$Trait <- c("Chronological Age","ProteAge Res/SD","KDM Res/SD","PhenoAge Res/SD")
dat2 <- dplyr::rename(dat2,mean = HR)
dat2 <- as.tibble(dat2)
dat2$N <- 2068
#write.csv(dat2,"result/mortality_prediction/male/Male_cohort_aged_40_70.csv")
#jpeg("result/mortality_prediction/male/Male_cohort_aged_40_70.jpeg",width = 20, height = 7, units="cm", res = 300, quality=100)
#dat2 |>
#  forestplot::forestplot(labeltext = c(Trait, HR2, p,N),
             #clip = c(0.1, 2.5),
             #xlog = TRUE,
#             zero = 1,
#             lty.ci = "solid",
#             title = "HR per 1-SD increase of BA Residual of male cohort (aged 40-70)",
#             graph.pos = 2,
#             boxsize = 0.3) |>
#  fp_add_lines() |>
#  fp_set_style(box = "royalblue",
#               line = "darkblue",
#               summary = "royalblue") |> 
#  fp_add_header(Trait = c("", "Characters"),
#                HR2 = c("", "HR"),
#                p = c("","P-Value"),
#                N = c("","N")) |>
#  fp_set_zebra_style("#EFEFEF")
#dev.off()

#########   40-50   #########
f21 <- dplyr::filter(f2,`Chronological Age` >= 40 & `Chronological Age` < 50)
cox21 <- coxph(Surv(yrs, status) ~ `Chronological Age` + `ProteAge Res/SD` + `KDM Res/SD` + `PhenoAge Res/SD`,
             data = f21)
tbl_regression(cox21,exp = TRUE)
c21 <- summary(cox21)
colnames(c21[["coefficients"]])
colnames(c21[["conf.int"]])

dat21 = as.data.frame(round(c21[["conf.int"]][, c(1, 3, 4)], 2))
dat21 = tibble::rownames_to_column(dat21, var = "Trait")
colnames(dat21)[2:4] = c("HR", "lower", "upper")
dat21$HR2 <- paste0(dat21[, 2], " (", dat21[, 3], "-", dat21[, 4], ")")
dat21$p <- round(c21[["coefficients"]][,5],3)
dat21$Trait <- c("Chronological Age","ProteAge Res/SD","KDM Res/SD","PhenoAge Res/SD")
dat21 <- dplyr::rename(dat21,mean = HR)
dat21 <- as.tibble(dat21)
dat21$N <- 120
#write.csv(dat21,"result/mortality_prediction/male/Male_cohort_aged_40_50.csv")
#jpeg("result/mortality_prediction/male/Male_cohort_aged_40_50.jpeg",width = 20, height = 7, units="cm", res = 300, quality=100)
#dat21 |>
#  forestplot::forestplot(labeltext = c(Trait, HR2, p,N),
             #clip = c(0.1, 2.5),
             #xlog = TRUE,
#             zero = 1,
#             lty.ci = "solid",
#             title = "HR per 1-SD increase of BA Residual of male cohort (aged 40-50)",
#             graph.pos = 2,
#             boxsize = 0.3) |>
#  fp_add_lines() |>
#  fp_set_style(box = "royalblue",
#               line = "darkblue",
#               summary = "royalblue") |> 
#  fp_add_header(Trait = c("", "Characters"),
#                HR2 = c("", "HR"),
#                p = c("","P-Value"),
#                N = c("","N")) |>
#  fp_set_zebra_style("#EFEFEF")
#dev.off()
#########   50-60   #########
f22 <- dplyr::filter(f2,`Chronological Age` >= 50 & `Chronological Age` < 60)
cox22 <- coxph(Surv(yrs, status) ~ `Chronological Age` + `ProteAge Res/SD` + `KDM Res/SD` + `PhenoAge Res/SD`,
             data = f22)
tbl_regression(cox22,exp = TRUE)
c22 <- summary(cox22)
colnames(c22[["coefficients"]])
colnames(c22[["conf.int"]])

dat22 = as.data.frame(round(c2[["conf.int"]][, c(1, 3, 4)], 2))
dat22 = tibble::rownames_to_column(dat22, var = "Trait")
colnames(dat22)[2:4] = c("HR", "lower", "upper")
dat22$HR2 <- paste0(dat22[, 2], " (", dat22[, 3], "-", dat22[, 4], ")")
dat22$p <- round(c22[["coefficients"]][,5],3)
dat22$Trait <- c("Chronological Age","ProteAge Res/SD","KDM Res/SD","PhenoAge Res/SD")
dat22 <- dplyr::rename(dat22,mean = HR)
dat22 <- as.tibble(dat22)
dat22$N <- 391
#write.csv(dat22,"result/mortality_prediction/male/Male_cohort_aged_50_60.csv")
#jpeg("result/mortality_prediction/male/Male_cohort_aged_50_60.jpeg",width = 20, height = 7, units="cm", res = 300, quality=100)
#dat22 |>
##  forestplot::forestplot(labeltext = c(Trait, HR2, p,N),
           #clip = c(0.1, 2.5),
             #xlog = TRUE,
#             zero = 1,
#             lty.ci = "solid",
#             title = "HR per 1-SD increase of BA Residual of male cohort (aged 50-60)",
#             graph.pos = 2,
#             boxsize = 0.3) |>
#  fp_add_lines() |>
#  fp_set_style(box = "royalblue",
#               line = "darkblue",
#               summary = "royalblue") |> 
#  fp_add_header(Trait = c("", "Characters"),
#                HR2 = c("", "HR"),
#                p = c("","P-Value"),
#                N = c("","N")) |>
#  fp_set_zebra_style("#EFEFEF")
#dev.off()
#########   60-70   #########
f23 <- dplyr::filter(f2,`Chronological Age` >= 60 & `Chronological Age` <= 70)
cox23 <- coxph(Surv(yrs, status) ~ `Chronological Age` + `ProteAge Res/SD` + `KDM Res/SD` + `PhenoAge Res/SD`,
             data = f23)
tbl_regression(cox23,exp = TRUE)
c23 <- summary(cox23)
colnames(c23[["coefficients"]])
colnames(c23[["conf.int"]])

dat23 = as.data.frame(round(c23[["conf.int"]][, c(1, 3, 4)], 2))
dat23 = tibble::rownames_to_column(dat23, var = "Trait")
colnames(dat23)[2:4] = c("HR", "lower", "upper")
dat23$HR2 <- paste0(dat23[, 2], " (", dat23[, 3], "-", dat23[, 4], ")")
dat23$p <- round(c23[["coefficients"]][,5],3)
dat23$Trait <- c("Chronological Age","ProteAge Res/SD","KDM Res/SD","PhenoAge Res/SD")
dat23 <- dplyr::rename(dat23,mean = HR)
dat23 <- as.tibble(dat23)
dat23$N <- 1575
#write.csv(dat23,"result/mortality_prediction/male/Male_cohort_aged_60_70.csv")
#jpeg("result/mortality_prediction/male/Male_cohort_aged_60_70.jpeg",width = 20, height = 7, units="cm", res = 300, quality=100)
#dat23 |>
#  forestplot::forestplot(labeltext = c(Trait, HR2, p,N),
#             #clip = c(0.1, 2.5),
#             #xlog = TRUE,
#             zero = 1,
#             lty.ci = "solid",
#             title = "HR per 1-SD increase of BA Residual of male cohort (aged 60-70)",
#             graph.pos = 2,
#             boxsize = 0.3) |>
#  fp_add_lines() |>
#  fp_set_style(box = "royalblue",
#               line = "darkblue",
#               summary = "royalblue") |> 
#  fp_add_header(Trait = c("", "Characters"),
#                HR2 = c("", "HR"),
#                p = c("","P-Value"),
#                N = c("","N")) |>
#  fp_set_zebra_style("#EFEFEF")
#dev.off()

#################   put it together
dat_m <- rbind(c("40-70 years old",NA,NA,NA,"","",""),dat2,
             c("40-50 years old",NA,NA,NA,"","",""),dat21,
             c("50-60 years old",NA,NA,NA,"","",""),dat22,
             c("60-70 years old",NA,NA,NA,"","",""),dat23) %>% as.data.frame
#dat[is.na(dat)] <- ""
dat_m[,c(2:4)] <- lapply(dat_m[,c(2:4)],as.numeric)
dat_m[,c(5,6,7)] <- lapply(dat_m[,c(5,6,7)],as.character)
dat_m$group <- "Male"
write.csv(dat_m,"result/aging_clock_mortality/male_cohort.csv")
#jpeg("result/mortality_prediction/male/male_cohort_aged.jpeg",width = 20, height = 20, units="cm", res = 300, quality=100)
#dat_m |>
#  forestplot::forestplot(labeltext = c(Trait, HR2, p,N),
#             #clip = c(0.1, 2.5),
#             xlog = TRUE,
#             zero = 1,
#             lty.ci = "solid",
#            title = "HR per 1-SD increase of BA Residual (Male Cohort)",
#             graph.pos = 2,
#             boxsize = 0.3,
#             is.summary = c(T,F,F,F,F,T,F,F,F,F,T,F,F,F,F,T,F,F,F,F)) |>
#  fp_add_lines() |>
#  fp_set_style(box = "royalblue",
#               line = "darkblue",
#               summary = "royalblue") |> 
#  fp_add_header(Trait = c("", "Characters"),
#                HR2 = c("", "HR"),
#                p = c("","P-Value"),
#                N = c("","N")) |>
#  fp_set_zebra_style("#EFEFEF")
#dev.off()


##################          
#dat_f$group <- ifelse(is.na(dat_f$mean),"",dat_f$group)
#dat_m$group <- ifelse(is.na(dat_m$mean),"",dat_m$group)
#dat_f$Trait <- paste0(dat_f$Trait," (",dat_f$group,")")
#dat_m$Trait <- paste0(dat_m$Trait," (",dat_m$group,")")


dat <- rbind(c("40-70 years old",NA,NA,NA,"","","",""),
             dat_f[2,],dat_m[2,],dat_f[3,],dat_m[3,],dat_f[4,],dat_m[4,],dat_f[5,],dat_m[5,],
             c("40-50 years old",NA,NA,NA,"","","",""),
             dat_f[7,],dat_m[7,],dat_f[8,],dat_m[8,],dat_f[9,],dat_m[9,],dat_f[10,],dat_m[10,],
             c("50-60 years old",NA,NA,NA,"","","",""),
             dat_f[12,],dat_m[12,],dat_f[13,],dat_m[13,],dat_f[14,],dat_m[14,],dat_f[15,],dat_m[15,],
             c("60-70 years old",NA,NA,NA,"","","",""),
             dat_f[17,],dat_m[17,],dat_f[18,],dat_m[18,],dat_f[19,],dat_m[19,],dat_f[20,],dat_m[20,])
#dat$group <- NULL
dat$` ` <- paste(rep(" ", 20), collapse = " ")
dat[,c(2:4)] <- lapply(dat[,c(2:4)],as.numeric)
dat[,c(5,6,7)] <- lapply(dat[,c(5,6,7)],as.character)
dat$Trait <- ifelse(is.na(dat$mean),dat$Trait,
                    paste0("   ", dat$Trait))
rownames(dat) <- c(1:36)
dat <- dplyr::rename(dat,`HR (95% CI) `= HR2)
tm <- forest_theme(base_size = 10,
                   ci_pch = 15,
                   ci_lty = 1)

                  
g <- forestploter::forest(dat[,c(1,8,9,5:7)],
       est = dat$mean,
       lower = dat$lower,
       upper = dat$upper,
       sizes = 0.5,
       ci_colum = 3,
       xlim = c(0.5,1.75),
       ticks_at = c(0.5, 1, 1.5),
       ref_line = 1,
       theme = tm)

g <- edit_plot(g,row = c(1,10,19,28),gp = gpar(fontface = "bold"))

g <- edit_plot(g,
               row = c(2,4,6,8,11,13,15,17,20,22,24,26,29,31,33,35),
               col = 3,
               which = "ci",
               gp = gpar(col = "#992224"))
g <- edit_plot(g,
               row = c(3,5,7,9,12,14,16,18,21,23,25,27,30,32,34,36),
               col = 3,
               which = "ci",
               gp = gpar(col = "#3E4F94"))
               
g <- edit_plot(g,
               row = c(2,4,6,8,11,13,15,17,20,22,24,26,29,31,33,35),
               gp = gpar(col = "#992224",fill = "#FFFACD"))
g <- edit_plot(g,
               row = c(3,5,7,9,12,14,16,18,21,23,25,27,30,32,34,36),
               gp = gpar(col = "#3E4F94",fill = "#E0FFFF"))

g <- edit_plot(g,
               row = c(2,4,6,8,11,13,15,17,20,22,24,26,29,31,33,35),
               which = "background",
               gp = gpar(fill = "#FFFFFF"))

g <- edit_plot(g,
               row = c(3,5,7,9,12,14,16,18,21,23,25,27,30,32,34,36),
               which = "background",
               gp = gpar(fill = "#FFFFFF"))
g <- edit_plot(g,
               row = c(1,10,19,28),
               which = "background",
               gp = gpar(fill = "#CFCFCF"))
g <- insert_text(g,
                 text = "HR per 1-SD increase of BA Residual",
                 col = 1:6,
                 part = "header",
                 gp = gpar(fontface = "bold"))
g
pdf("1st submission NC/Figures JZH 20241226/Figure 3/Figure 3.pdf",width = 7, height = 9.5)   
g
dev.off()
######## K-M plot
#ukb_bioage_whole$ager <- ifelse(ukb_bioage_whole$`ProteAge Residual` > 0, "Fast Ager","Slow Ager")
#k0 <- merge(ukb_bioage_whole,f,by = "eid")

#km <- survfit(Surv(yrs, status) ~ ager,data = k0,type = "kaplan-meier",conf.type = "log")
#ggsurvplot(km,
#           main = "Survival curve",
#           pval = TRUE)
#ggview(w = 5,h = 5)
#ggsave("submit/figure/figure3/KM_plot.jpeg",w = 4, h = 3)





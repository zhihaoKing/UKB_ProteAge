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
library(scales)



################################  select candidate biomarkers ###########################
biomarker_labels <- read_excel("data/BioAge_variables.xlsx", sheet = "UKB_variables") %>% as.data.frame() 
# select biomarkers
describe(BioAge::NHANES3)
data.frame(missing_percent=sapply(NHANES3, function(x) round(mean(is.na(x))*100, 1))) # proportion of missing data
names(which(colMeans(is.na(NHANES3)) > 0.2)) # several biomarkers have high missing, some with even 100% missing -> those missing >20% were not selected in the algorithm

### Select biomarkers to use in NHANES3 (select those with absolute correlation >0.1)
candidate_biomarkers_nhanes3 <- c("bmi","waist","fev_1000",
                                  "albumin_gL","alp","bun","creat_umol","glucose_mmol",
                                  "ttbl","uap","lymph","mcv","monopa",
                                  "rbc","rdw","wbc","crp","dbp","sbp",
                                  "pulse","hba1c","hdl","trig","totchol"
                                  #"grip_r","grip_l","basopa","eosnpa","neut","cyst","ggt","ldl" # Exclude the ones with high missingness
)


NHANES3_m <- dplyr::filter(NHANES3, gender == 1)
NHANES3_f <- dplyr::filter(NHANES3, gender == 2)

# Calculate Pearson correlations for all, male, and female datasets
pearson_m_nhanes <- data.frame(BM=candidate_biomarkers_nhanes3, r="", Selected=FALSE)
pearson_f_nhanes <- data.frame(BM=candidate_biomarkers_nhanes3, r="", Selected=FALSE)
pearson_all_nhanes <- data.frame(BM=candidate_biomarkers_nhanes3, r="", Selected=FALSE)

for (i in 1:length(candidate_biomarkers_nhanes3)){
  BM <- candidate_biomarkers_nhanes3[i]
  
  # All
  r_all <- cor(NHANES3[,BM], NHANES3$age, use="pairwise.complete.obs")
  pearson_all_nhanes[pearson_all_nhanes$BM==BM, "r"] <- r_all %>% round(digits=2)
  if(abs(r_all)>0.1 & !is.na(r_all)){pearson_all_nhanes[pearson_all_nhanes$BM==BM, "Selected"] <- TRUE}
  
  # Male
  r_m <- cor(NHANES3_m[,BM], NHANES3_m$age, use="pairwise.complete.obs")
  pearson_m_nhanes[pearson_m_nhanes$BM==BM, "r"] <- r_m %>% round(digits=2)
  if(abs(r_m)>0.1 & !is.na(r_m)){pearson_m_nhanes[pearson_m_nhanes$BM==BM, "Selected"] <- TRUE}
  
  # Female  
  r_f <- cor(NHANES3_f[,BM], NHANES3_f$age, use="pairwise.complete.obs")
  pearson_f_nhanes[pearson_f_nhanes$BM==BM, "r"] <- r_f %>% round(digits=2)
  if(abs(r_f)>0.1 & !is.na(r_f)){pearson_f_nhanes[pearson_f_nhanes$BM==BM, "Selected"] <- TRUE}
}
#rm(list=c("BM","r_all","r_m","r_f","i","NHANES3_m","NHANES3_f"))

biomarkers_all_nhanes <- data.frame(BM=candidate_biomarkers_nhanes3[pearson_all_nhanes$Selected],
                                    label=biomarker_labels$Description[match(candidate_biomarkers_nhanes3[pearson_all_nhanes$Selected],biomarker_labels$`Variable Name`)])
pearson_all_nhanes[pearson_all_nhanes$BM %in% biomarkers_all_nhanes[,1],] #
round(cor(NHANES3[,biomarkers_all_nhanes[,1]],use="pairwise.complete.obs"), 2)
biomarkers_all_nhanes <- subset(biomarkers_all_nhanes, biomarkers_all_nhanes$BM!="pulse")

################################  Training BioAge in NHANES   ###########################

### Train KDM bioage in NHANES III (separate training for men and women)
# Descriptive statistics of the dataset used for training of KDM
NHANES3 %>% 
  filter(age >= 30 & age <= 75 & pregnant == 0) %>% # Reference sample is NHANES III nonpregnant participants aged 30-75 years
  select(sampleID, year, wave, gender, age, biomarkers_all_nhanes[,1]) %>%
  na.omit() %>%
  describe()
# KDM algorithm based on 18 newly selected clinical biomarkers
kdm_nhanes_trained <- kdm_nhanes(biomarkers=biomarkers_all_nhanes[,1])
# KDM algorithm using original biomarkers in Levine et al (PMID: 23213031)
kdm_nhanes_trained_original <- kdm_nhanes(biomarkers=c("fev_1000","sbp","bun","hba1c","totchol","creat_umol","albumin_gL","alp","crp"))
# For sensitivity analysis for KDM, exclude HbA1c and serum glucose
kdm_nhanes_trained_noglu <- kdm_nhanes(biomarkers=biomarkers_all_nhanes[-c(7,16),1])



### Train PhenoAge in NHANES III
# Descriptive statistics of the dataset used for training of PhenoAge
NHANES3 %>% 
  filter(age >= 20 & age <= 84) %>%
  select(sampleID, year, wave, gender, age, biomarkers_all_nhanes[,1]) %>%
  na.omit() %>%
  describe()
# PhenoAge based on 18 newly selected clinical biomarkers
phenoage_nhanes_trained <- phenoage_nhanes(biomarkers=biomarkers_all_nhanes[,1])
# PhenoAge using original biomarkers in Levine et al
phenoage_nhanes_trained_original <- phenoage_nhanes(biomarkers=c("creat_umol","glucose_mmol","rdw","albumin_gL","alp","mcv","lymph","crp","wbc"))
# Sensitivity analysis for PhenoAge, exclude HbA1c and serum glucose
phenoage_nhanes_trained_noglu <- phenoage_nhanes(biomarkers=biomarkers_all_nhanes[-c(7,16),1])



######### Step 2: Testing BA algorithms in NHANES IV #########

### Assemble NHANES IV dataset with projected biological aging measures for analysis
BioAge_nhanes_trained_data <- merge(phenoage_nhanes_trained$data, kdm_nhanes_trained$data) %>% 
  merge(., kdm_nhanes_trained_original$data %>% dplyr::rename(kdm_original=kdm,kdm_advance_original=kdm_advance)) %>% 
  merge(., phenoage_nhanes_trained_original$data %>% dplyr::rename(phenoage_original=phenoage,phenoage_advance_original=phenoage_advance)) %>%
  merge(., kdm_nhanes_trained_noglu$data %>% dplyr::rename(kdm_noglu=kdm,kdm_advance_noglu=kdm_advance)) %>% 
  merge(., phenoage_nhanes_trained_noglu$data %>% dplyr::rename(phenoage_noglu=phenoage,phenoage_advance_noglu=phenoage_advance))

### Create BA residuals using regression model of a natural spline of CA with 3 degrees of freedom
get_BA_resids <- function(BA){
  data = BioAge_nhanes_trained_data %>% drop_na(BA)
  model <- parse(text = sprintf("lm(%s ~ ns(age, df = 3), data = data)", BA)) %>% eval()
  model_predict <- ggpredict(model, terms = c("age"))
  data[,"BA_res"] <- NA
  data[!is.na(data[BA]),"BA_res"] <- resid(model)
  return(residuals(model))
}
for(BA in c("kdm","kdm_original","kdm_noglu","phenoage","phenoage_original","phenoage_noglu")){
  BA_res <- paste0(BA, "_res")
  BioAge_nhanes_trained_data[,BA_res] = NA
  BioAge_nhanes_trained_data[!is.na(BioAge_nhanes_trained_data[BA]),BA_res] <- get_BA_resids(BA)
}
#rm(list=c("BA","BA_res"))

summary(BioAge_nhanes_trained_data %>% 
  select("kdm_original","kdm_original_res","kdm","kdm_res","kdm_noglu","kdm_noglu_res",
         "phenoage_original","phenoage_original_res","phenoage","phenoage_res","phenoage_noglu","phenoage_noglu_res"))


### Checking of the modified BioAge measures in NHANES4

# Correlations between BioAge measures and CA
agevar_1 <- c(
  "kdm"="KDM",
  "phenoage"="PhenoAge",
  #"hd_log"="HD (log)",
  "kdm_original"="Levine original\nKDM",
  "phenoage_original"="Levine original\nPhenoAge")
#jpeg(file="result/BioAge/NHANES4_BioAge_CA_correlations.jpg", width = 22, height = 15, units="cm", res = 300, quality=100)
#plot_ba(data=BioAge_nhanes_trained_data, agevar=names(agevar_1), label=as.vector(agevar_1))
#dev.off()

# Correlations between age residuals and CA
agevar_2 <- c(
  "kdm_res"="KDM residual",
  "phenoage_res"="PhenoAge residual",
  #"hd_log" = "HD (log)",
  "kdm_original_res"="Levine original\nKDM residual",
  "phenoage_original_res"="Levine original\nPhenoAge residual",
  "age"="Chronological age")
get_axis_type <- function(labels){
  return(rep("float", length(labels)) %>% setNames(names(labels))) # Create function to generate axis_type variables for BAA plots
}

################################ Upload UK Biobank data ###########################
ukb <- read.csv("data/protein_cohort_phenoage.csv")
#ukb$Sex <- ifelse(ukb$Sex == "1","Male","Female") ## 502367 obs,35 variables

colnames(ukb) <- c("eid","age","sex","bmi","waist","fev_1000","grip_l",
                  "albumin_gL","alp","bun_mmol","creat_umol","glucose_mmol","ttbl_umol",
                  "uap_umol","basopa","eosnpa","lymph","mcv","monopa","neut","rbc_L",
                  "rdw","wbc_L","crp_mgL","cyst","dbp","sbp","pulse","ggt","hba1c_mmol",
                  "hdl_mmol","ldl_mmol","trig_mmol","totchol_mmol","grip_r")

ukb$lnalp <- log(ukb$alp)
ukb$waist <- ukb$waist
ukb$fev <- ukb$fev_1000*1000
ukb$bun <- ukb$bun_mmol*2.8
ukb$lnbun <- log(ukb$bun)
ukb$lncreat_umol <- log(ukb$creat_umol)
ukb$ttbl <- ukb$ttbl_umol*.0585
ukb$uap <- ukb$uap_umol*.0168
ukb$rbc <- ukb$rbc_L
ukb$wbc <- ukb$wbc_L
ukb$crp <- ukb$crp_mgL/10
ukb$lncrp <- log(ukb$crp_mgL)
ukb$hba1c <- ukb$hba1c_mmol*0.0915 + 2.15
ukb$lnhba1c <- log(ukb$hba1c)
ukb$hdl <- ukb$hdl_mmol*38.67
ukb$ldl <- ukb$ldl_mmol*38.67
ukb$trig <- ukb$trig_mmol*88.57
ukb$totchol <- ukb$totchol_mmol*38.67


######### Step 3: Projecting BA algorithms onto UKB #########

#### Missing data pattern in UKB 
data.frame(missing_percent=sapply(ukb[,biomarkers_all_nhanes[,1]], function(x) round(mean(is.na(x))*100, 1))) # proportion of missing data
table(rowSums(is.na(ukb[,biomarkers_all_nhanes[,1]]))) # Number of missing biomarkers per individual; n=331,699 with all the 18 biomarkers available

### Projecting trained BA measures onto UKB (here we excluded UKB participants missing any of the 18 included biomarkers)

# KDM using 18 newly selected biomarkers (separate training for gender)
ukb_kdm_m <- kdm_calc(data = ukb[complete.cases( ukb[, biomarkers_all_nhanes[,1] ]), ] %>% 
                        filter(sex == "Male"), 
                      biomarkers=biomarkers_all_nhanes[,1],
                      fit = kdm_nhanes_trained$fit$male, 
                      s_ba2 = kdm_nhanes_trained$fit$male$s_b2)$data
ukb_kdm_f <- kdm_calc(data = ukb[complete.cases( ukb[, biomarkers_all_nhanes[,1] ]), ] %>% 
                        filter(sex == "Female"), 
                      biomarkers=biomarkers_all_nhanes[,1],
                      fit = kdm_nhanes_trained$fit$female, 
                      s_ba2 = kdm_nhanes_trained$fit$female$s_b2)$data
ukb_kdm_all <- rbind(ukb_kdm_m, ukb_kdm_f) # Combine the KDM datasets ## 331671 obs, 54 variables
rm(list=c("ukb_kdm_m","ukb_kdm_f"))

# Levine original KDM
ukb_kdm_original_m <- kdm_calc(data = ukb[complete.cases( ukb[, biomarkers_all_nhanes[,1] ]), ] %>%
                                 filter(sex == "Male"), 
                               biomarkers=c("fev_1000","sbp","bun","hba1c","totchol","creat_umol","albumin_gL","alp","crp"),
                               fit = kdm_nhanes_trained_original$fit$male, 
                               s_ba2 = kdm_nhanes_trained_original$fit$male$s_b2)$data
ukb_kdm_original_f <- kdm_calc(data = ukb[complete.cases( ukb[, biomarkers_all_nhanes[,1] ]), ] %>% 
                                 filter(sex == "Female"), 
                               biomarkers=c("fev_1000","sbp","bun","hba1c","totchol","creat_umol","albumin_gL","alp","crp"),
                               fit = kdm_nhanes_trained_original$fit$female, 
                               s_ba2 = kdm_nhanes_trained_original$fit$female$s_b2)$data
ukb_kdm_original_all <- rbind(ukb_kdm_original_m, ukb_kdm_original_f) %>% # Combine the KDM datasets 
  dplyr::rename(kdm_original=kdm, kdm_advance_original=kdm_advance) ## 331671 obs, 54 variables
rm(list=c("ukb_kdm_original_m","ukb_kdm_original_f"))

# Sensitivity analysis for KDM, remove HbA1c and serum glucose
ukb_kdm_noglu_m <- kdm_calc(data = ukb[complete.cases( ukb[, biomarkers_all_nhanes[,1] ]), ] %>%
                              filter(sex == "Male"), 
                            biomarkers=biomarkers_all_nhanes[-c(7,16),1],
                            fit = kdm_nhanes_trained_noglu$fit$male, 
                            s_ba2 = kdm_nhanes_trained_noglu$fit$male$s_b2)$data
ukb_kdm_noglu_f <- kdm_calc(data = ukb[complete.cases( ukb[, biomarkers_all_nhanes[,1] ]), ] %>% 
                              filter(sex == "Female"), 
                            biomarkers=biomarkers_all_nhanes[-c(7,16),1],
                            fit = kdm_nhanes_trained_noglu$fit$female, 
                            s_ba2 = kdm_nhanes_trained_noglu$fit$female$s_b2)$data
ukb_kdm_noglu_all <- rbind(ukb_kdm_noglu_m, ukb_kdm_noglu_f) %>% # Combine the KDM datasets
  dplyr::rename(kdm_noglu=kdm, kdm_advance_noglu=kdm_advance) ## 331671 obs, 54 variables
rm(list=c("ukb_kdm_noglu_m","ukb_kdm_noglu_f"))


# PhenoAge using 18 newly selected biomarkers
ukb_pheno_all <- phenoage_calc(data = ukb[complete.cases( ukb[, biomarkers_all_nhanes[,1] ]), ], 
                               biomarkers=biomarkers_all_nhanes[,1], 
                               fit = phenoage_nhanes_trained$fit)$data
dim(subset(ukb_pheno_all, ukb_pheno_all$phenoage==Inf))[1] # 34 individuals with phenoage=inf
ukb_pheno_all[which(ukb_pheno_all$phenoage==Inf),c("phenoage", "phenoage_advance")] <- NA  # Exclude individuals with infinite PhenoAge

# Levine original PhenoAge
ukb_pheno_original_all <- phenoage_calc(data = ukb[complete.cases( ukb[, biomarkers_all_nhanes[,1] ]), ], 
                                        biomarkers=c("creat_umol","glucose_mmol","rdw","albumin_gL","alp","mcv","lymph","crp","wbc"), 
                                        fit = phenoage_nhanes_trained_original$fit)$data %>%
  dplyr::rename(phenoage_original=phenoage, phenoage_advance_original=phenoage_advance)
dim(subset(ukb_pheno_original_all, ukb_pheno_original_all$phenoage_original==Inf))[1] # 40 individuals with phenoage=inf
ukb_pheno_original_all[which(ukb_pheno_original_all$phenoage_original==Inf),c("phenoage_original", "phenoage_advance_original")] <- NA  # Exclude individuals with infinite PhenoAge

# Sensitivity analysis for PhenoAge, remove HbA1c and serum glucose
ukb_pheno_noglu_all <- phenoage_calc(data = ukb[complete.cases( ukb[, biomarkers_all_nhanes[,1] ]), ], 
                                     biomarkers=biomarkers_all_nhanes[-c(7,16),1], 
                                     fit = phenoage_nhanes_trained_noglu$fit)$data %>%
  dplyr::rename(phenoage_noglu=phenoage, phenoage_advance_noglu=phenoage_advance)
dim(subset(ukb_pheno_noglu_all, ukb_pheno_noglu_all$phenoage_noglu==Inf))[1] # 37 individuals with phenoage=inf
ukb_pheno_noglu_all[which(ukb_pheno_noglu_all$phenoage_noglu==Inf),c("phenoage_noglu", "phenoage_advance_noglu")] <- NA  # Exclude individuals with infinite PhenoAge

### Merge all biological aging measures
ukb_bioage <- left_join(ukb[complete.cases( ukb[, biomarkers_all_nhanes[,1] ]), ] %>% 
                          select("eid","sex","age",biomarkers_all_nhanes[,1]), 
                        ukb_kdm_all[, c("eid", "kdm")], by = "eid") %>%
  left_join(., ukb_kdm_original_all[, c("eid", "kdm_original")], by = "eid") %>%
  left_join(., ukb_kdm_noglu_all[, c("eid", "kdm_noglu")], by = "eid") %>%
  left_join(., ukb_pheno_all[, c("eid","phenoage")], by = "eid") %>%
  left_join(., ukb_pheno_original_all[, c("eid", "phenoage_original")], by = "eid") %>%
  left_join(., ukb_pheno_noglu_all[, c("eid", "phenoage_noglu")], by = "eid")
summary(ukb_bioage %>% select(kdm, kdm_original, kdm_noglu, phenoage, phenoage_original, phenoage_noglu)) 


colnames(ukb_bioage)[2] <- "gender"
### Exclude extreme outliers (those outside 5SD from mean)
table(abs(as.vector(scale(ukb_bioage$kdm)))>5) # n=55 with KDM outside 5 SD
table(abs(as.vector(scale(ukb_bioage$kdm_original)))>5) # n=32 with KDM (original) outside 5 SD
table(abs(as.vector(scale(ukb_bioage$kdm_noglu)))>5) # n=56 with KDM (sensitivity analysis) outside 5 SD
table(abs(as.vector(scale(ukb_bioage$phenoage)))>5) # n=45 with PhenoAge outside 5 SD
table(abs(as.vector(scale(ukb_bioage$phenoage_original)))>5) # n=79 with PhenoAge (original) outside 5 SD
table(abs(as.vector(scale(ukb_bioage$phenoage_noglu)))>5) # n=48 with PhenoAge (sensitivity analysis) outside 5 SD


ukb_bioage$kdm <- ifelse(abs(as.vector(scale(ukb_bioage$kdm)))>5, NA, ukb_bioage$kdm)
ukb_bioage$kdm_original <- ifelse(abs(as.vector(scale(ukb_bioage$kdm_original)))>5, NA, ukb_bioage$kdm_original)
ukb_bioage$kdm_noglu <- ifelse(abs(as.vector(scale(ukb_bioage$kdm_noglu)))>5, NA, ukb_bioage$kdm_noglu)
ukb_bioage$phenoage <- ifelse(abs(as.vector(scale(ukb_bioage$phenoage)))>5, NA, ukb_bioage$phenoage)
ukb_bioage$phenoage_original <- ifelse(abs(as.vector(scale(ukb_bioage$phenoage_original)))>5, NA, ukb_bioage$phenoage_original)
ukb_bioage$phenoage_noglu <- ifelse(abs(as.vector(scale(ukb_bioage$phenoage_noglu)))>5, NA, ukb_bioage$phenoage_noglu)

summary(ukb_bioage %>% select(kdm, kdm_original, kdm_noglu, phenoage, phenoage_original, phenoage_noglu)) 

### Regression model for BioAge, using a natural spline of CA with 3 degrees of freedom (to create age residuals)
# Function to generate BA residuals and plot model used
#get_BA_resids <- function(BA){
#  data = ukb_bioage %>% drop_na(BA)
  # Basic model = regress on age alone
#  model <- parse(text = sprintf("lm(%s ~ ns(age, df = 3), data = data)", BA)) %>% eval()
#  model_predict <- ggpredict(model, terms = c("age"))
#  data[,"BA_res"] <- NA
#  data[!is.na(data[BA]),"BA_res"] <- resid(model)
#  return(residuals(model))
#}
#for(BA in c("kdm","kdm_original","kdm_noglu","phenoage","phenoage_original","phenoage_noglu")){
#  BA_res <- paste0(BA, "_res")
#  ukb_bioage[,BA_res] = NA
#  ukb_bioage[!is.na(ukb_bioage[BA]),BA_res] <- get_BA_resids(BA)
#}
write.csv(ukb_bioage,"result/bioage/UKB_BioAge_all_features.csv")
#write.csv(ukb_bioage,"submit/supplymentary/supplymentary_table/predicted_age/supplymentary_table_predicted_age_of_ProteAge_KDM_PhenoAge.csv")





################################ compare BioAge with proteomics aging clock #################################
rm(list=ls())
ukb_bioage <- read.csv("result/bioage/UKB_BioAge_all_features.csv")
ukb_bioage$X <- NULL

########### function
get_BA_resids <- function(ukb_bioage,BA){
  data = ukb_bioage %>% drop_na(BA)
  # Basic model = regress on age alone
  model <- parse(text = sprintf("lm(%s ~ ns(age, df = 3), data = data)", BA)) %>% eval()
  model_predict <- ggpredict(model, terms = c("age"))
  data[,"BA_res"] <- NA
  data[!is.na(data[BA]),"BA_res"] <- resid(model)
  return(residuals(model))
}

################################## female cohort ##################################

sm_statCorr <- function (..., fit.params = list(), corr_method = "pearson", 
    alternative = "two.sided", separate_by = ",", label_x = NULL, 
    label_y = NULL, text_size = 4, show_text = TRUE, borders = TRUE, 
    legends = FALSE, r2 = FALSE, R2) 
{
    if (!missing(R2)) {
        r2 <- R2
    }
    params <- list(...)
    fit.params <- modifyList(params, fit.params)
    fitPlot <- do.call("geom_smooth", modifyList(list(method = "lm", 
        se = FALSE, alpha = 0.2, weight = 0.8), fit.params))
    if (r2 == FALSE) {
        textPlot <- ggpubr::stat_cor(p.accuracy = 0.001, method = corr_method, 
            alternative = alternative, label.sep = separate_by, 
            label.x = label_x, label.y = label_y, size = text_size,digits = 3)
    }
    else {
        if (separate_by == "\n") {
            textPlot <- ggpubr::stat_cor(aes(label = paste(paste0("atop(", 
                after_stat(rr.label), ",", after_stat(p.label), 
                ")"), sep = "~~")), p.accuracy = 0.001, method = corr_method, 
                alternative = alternative, label.x = label_x, 
                label.y = label_y, size = text_size,digits = 3)
        }
        else {
            textPlot <- ggpubr::stat_cor(aes(label = paste(after_stat(rr.label), 
                after_stat(p.label), sep = paste0("~`", separate_by, 
                  "`~"))), p.accuracy = 0.001, method = corr_method, 
                alternative = alternative, label.x = label_x, 
                label.y = label_y, size = text_size,digits = 3)
        }
    }
    if (show_text == FALSE) {
        textPlot <- NULL
    }
    list(fitPlot, textPlot, sm_hvgrid(borders = borders, legends = legends))
}



p_f <- read.csv("result/ProteAge/female_proteage_validation_cohort.csv")
p_f$X <- NULL
ukb_bioage_f <- merge(p_f,ukb_bioage,by = "eid") %>% dplyr::rename(age = age.x) ## 5740 obs, 35 variables
ukb_bioage_f$age.y <- NULL

#}
ukb_bioage_f <- ukb_bioage_f %>% dplyr::select(eid,age,gender,proteage,kdm_original,phenoage_original) %>%
                    dplyr::rename(ProteAge = proteage,
                           KDM = kdm_original,
                           PhenoAge = phenoage_original)
set.seed(777)
for(BA in c("ProteAge","KDM","PhenoAge")){
  BA_res <- paste0(BA, " res")
  ukb_bioage_f[,BA_res] = NA
  ukb_bioage_f[!is.na(ukb_bioage_f[BA]),BA_res] <- get_BA_resids(ukb_bioage_f,BA)
}
write.csv(ukb_bioage_f,"result/bioage/female_cohort.csv")
#show_col(pal_npg("nrc")(10))
ukb_bioage_f <- na.omit(ukb_bioage_f)
# female proteage #992224
# female KDM #E3625D
# female phenoage #F0C284
# male proteage #3E4F94
# male KDM #3E90BF
# male phenoage #58B6E9
mae <- mae(ukb_bioage_f$age,ukb_bioage_f$ProteAge)  ## 2.238
ggplot(ukb_bioage_f,aes(x = age,y = ProteAge)) + 
geom_point(color = "#992224",size = 1,alpha = 0.8) + 
sm_statCorr(color = "#030303", corr_method = "pearson",text_size = 5.5) +
theme_classic(base_size = 14) + 
theme(legend.title = element_blank(),
      legend.position="bottom",
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 18,face = "bold"),
      plot.title = element_text(size = 22,hjust = 0.5,face = "bold")) +
xlim(c(40,75)) +
ylim(c(40,75)) +
labs(x="Chronological Age", y="ProteAge",title = "Female cohort") 
ggview(w = 5.5, h = 5.5) 
ggsave("1st submission NC/Figures JZH 20241226/figure 2/Figure 2A.pdf",w = 5.5, h = 5.5,dpi = 300) 

## proteage and kdm and phenoage
mae <- mae(ukb_bioage_f$ProteAge,ukb_bioage_f$KDM)  ## 3.569
ggplot(ukb_bioage_f,aes(x = ProteAge,y = KDM)) + 
geom_point(aes(colour = age),size = 1,alpha = 0.5) +
sm_statCorr(color = "#030303", corr_method = "pearson",text_size = 5.5) +
theme_classic(base_size = 14) + 
theme(legend.position="bottom",
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 18,face = "bold"),
      plot.title = element_text(size = 22,hjust = 0.5,face = "bold")) +
scale_color_distiller(palette = "RdBu") + 
xlim(c(40,75)) +
ylim(c(40,75)) +
labs(x="ProteAge", y="KDM",color = "Chronological Age",title = "Female cohort") 
ggview(w = 5.5, h = 6.7)
ggsave("1st submission NC/Figures JZH 20241226/figure 2/Figure 2E.pdf",w = 5.5, h = 6.7,dpi = 300)  

mae <- mae(ukb_bioage_f$ProteAge,ukb_bioage_f$PhenoAge)  ## 7.453
ggplot(ukb_bioage_f,aes(x = ProteAge,y = PhenoAge)) + 
geom_point(aes(colour = age),size = 1,alpha = 0.5) +
sm_statCorr(color = "#030303", corr_method = "pearson",text_size = 5.5) +
theme_classic(base_size = 14) + 
theme(legend.position="bottom",
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 18,face = "bold"),
      plot.title = element_text(size = 22,hjust = 0.5,face = "bold")) +
scale_color_distiller(palette = "RdBu") + 
xlim(c(40,75)) +
ylim(c(40,75)) +
labs(x="ProteAge", y="PhenoAge",color = "Chronological Age",title = "Female cohort") 
ggview(w = 5.5, h = 6.7)
ggsave("1st submission NC/Figures JZH 20241226/figure 2/Figure 2F.pdf",w = 5.5, h = 6.7,dpi = 300)  







################################## male cohort ##################################
p_m <- read.csv("result/ProteAge/male_proteage_validation_cohort.csv")
p_m$X <- NULL
ukb_bioage_m <- merge(p_m,ukb_bioage,by = "eid") %>% dplyr::rename(age = age.x) ## 4921 obs, 35 variables
ukb_bioage_m$age.y <- NULL ## 4921 obs, 35 variables

ukb_bioage_m <- ukb_bioage_m %>% dplyr::select(eid,age,gender,proteage,kdm_original,phenoage_original) %>%
                    dplyr::rename(ProteAge = proteage,
                           KDM = kdm_original,
                           PhenoAge = phenoage_original)
set.seed(888)
for(BA in c("ProteAge","KDM","PhenoAge")){
  BA_res <- paste0(BA, " res")
  ukb_bioage_m[,BA_res] = NA
  ukb_bioage_m[!is.na(ukb_bioage_m[BA]),BA_res] <- get_BA_resids(ukb_bioage_m,BA)
}
write.csv(ukb_bioage_m,"result/bioage/male_cohort.csv")
#show_col(pal_npg("nrc")(10))
ukb_bioage_m <- na.omit(ukb_bioage_m)
# female proteage #992224
# female KDM #E3625D
# female phenoage #F0C284
# male proteage #3E4F94
# male KDM #3E90BF
# male phenoage #58B6E9

mae <- mae(ukb_bioage_m$age,ukb_bioage_m$ProteAge) ## 2.444
ggplot(ukb_bioage_m,aes(x = age,y = ProteAge)) + 
geom_point(color = "#3E4F94",size = 1,alpha = 0.8) + 
sm_statCorr(color = "#030303", corr_method = "pearson",text_size = 5.5) +
theme_classic(base_size = 14) + 
theme(legend.title = element_blank(),
      legend.position="bottom",
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 18,face = "bold"),
      plot.title = element_text(size = 22,hjust = 0.5,face = "bold")) +
xlim(c(40,75)) +
ylim(c(40,75)) +
labs(x="Chronological Age", y="ProteAge",title = "Male cohort") 
ggview(w = 5.5, h = 5.5)
ggsave("1st submission NC/Figures JZH 20241226/figure 2/Figure 2B.pdf",w = 5.5, h = 5.5,dpi = 300)      
 


## proteage and kdm and phenoage
mae <- mae(ukb_bioage_m$ProteAge,ukb_bioage_m$KDM) ## 4.554
ggplot(ukb_bioage_m,aes(x = ProteAge,y = KDM)) + 
geom_point(aes(colour = age),size = 1,alpha = 0.5) +
sm_statCorr(color = "#030303", corr_method = "pearson",text_size = 5.5) +
theme_classic(base_size = 14) + 
theme(legend.position="bottom",
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 18,face = "bold"),
      plot.title = element_text(size = 22,hjust = 0.5,face = "bold")) +
scale_color_distiller(palette = "RdBu") + 
xlim(c(40,75)) +
ylim(c(40,75)) +
labs(x="ProteAge", y="KDM",color = "Chronological Age",title = "Male Cohort") 
ggview(w = 5.5, h = 6.7)
ggsave("1st submission NC/Figures JZH 20241226/figure 2/Figure 2G.pdf",w = 5.5, h = 6.7,dpi = 300)  

mae <- mae(ukb_bioage_m$ProteAge,ukb_bioage_m$PhenoAge) ## 6.372
ggplot(ukb_bioage_m,aes(x = ProteAge,y = PhenoAge)) + 
geom_point(aes(colour = age),size = 1,alpha = 0.5) +
sm_statCorr(color = "#030303", corr_method = "pearson",text_size = 5.5) +
theme_classic(base_size = 14) + 
theme(legend.position="bottom",
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 18,face = "bold"),
      plot.title = element_text(size = 22,hjust = 0.5,face = "bold")) +
scale_color_distiller(palette = "RdBu") + 
xlim(c(40,75)) +
ylim(c(40,75)) +
labs(x="ProteAge", y="PhenoAge",color = "Chronological Age",title = "Male Cohort") 
ggview(w = 5.5, h = 6.7)
ggsave("1st submission NC/Figures JZH 20241226/figure 2/Figure 2H.pdf",w = 5.5, h = 6.7,dpi = 300) 


############################# evaluation metrics
#填补缺失值
ukb_bioage_f$KDM <- impute(ukb_bioage_f$KDM, mean)
ukb_bioage_f$PhenoAge <- impute(ukb_bioage_f$PhenoAge, mean)
ukb_bioage_m$KDM <- impute(ukb_bioage_m$KDM, mean)
ukb_bioage_m$PhenoAge <- impute(ukb_bioage_m$PhenoAge, mean)

### MAE
## male
mae_female_proage <- mae(ukb_bioage_f$age,ukb_bioage_f$ProteAge) ## 2.238
mae_female_kdm <- mae(ukb_bioage_f$age,ukb_bioage_f$KDM) ## 2.906
mae_female_phenoage <- mae(ukb_bioage_f$age,ukb_bioage_f$PhenoAge)  ## 7.256
## female
mae_male_proage <- mae(ukb_bioage_m$age,ukb_bioage_m$ProteAge) ## 2.444
mae_male_kdm <- mae(ukb_bioage_m$age,ukb_bioage_m$KDM) ## 3.979
mae_male_phenoage <- mae(ukb_bioage_m$age,ukb_bioage_m$PhenoAge) ## 6.081


MAE <- data.frame(mae_female_proage,mae_female_kdm,mae_female_phenoage,
                  mae_male_proage,mae_male_kdm,mae_male_phenoage)
colnames(MAE) <- c("Female ProteAge","Female KDM","Female PhenoAge",
                   "Male ProteAge","Male KDM","Male PhenoAge")
MAE <- t(MAE) %>% as.data.frame
colnames(MAE) <- c("Value")
MAE$BA <- rownames(MAE)
MAE$group <- c("ProteAge","KDM","PhenoAge","ProteAge","KDM","PhenoAge")
MAE$sex <- c("Female","Female","Female","Male","Male","Male")
MAE$BA <- factor(MAE$BA,
                 levels = c("Female ProteAge","Male ProteAge","Female KDM",
                            "Male KDM","Female PhenoAge","Male PhenoAge"))
MAE$group <- factor(MAE$group,levels = c("ProteAge","KDM","PhenoAge"))

#show_col(pal_npg("nrc")(10))

# female proteage #992224
# female KDM #E3625D
# female phenoage #F0C284
# male proteage #3E4F94
# male KDM #3E90BF
# male phenoage #58B6E9

ggplot(MAE,aes(x = BA,y = Value)) + 
geom_bar(aes(fill = BA),stat = "identity",show.legend = F,lwd = 0.7) + 
theme_classic(base_size = 14) + 
scale_fill_npg() + 
theme(axis.text.y = element_text(size = 15),
      plot.title = element_text(size = 30,hjust = 0.5,face = "bold"),
      axis.text.x = element_text(angle = 15,vjust = 0.6,size = 14.2)) +
scale_fill_manual(values = c(`Female ProteAge` = "#992224",
                              `Male ProteAge` = "#3E4F94",
                              `Female KDM` = "#E3625D",
                              `Male KDM` = "#3E90BF",
                              `Female PhenoAge` = "#F0C284",
                              `Male PhenoAge` = "#58B6E9")) +
labs(x="", y="",title = "Mean Absolute Error") 
ggview(w = 7,h = 6)
ggsave("1st submission NC/Figures JZH 20241226/figure 2/Figure 2K.pdf",w = 7,h = 6, dpi = 300)


#ggplot(MAE[c(1,4),],aes(x = sex,y = Value,fill = BA)) + 
#geom_col(show.legend = F,fill = c("#E71F19","#F4A016"),lwd = 0.7) + 
#theme_classic(base_size = 14) + 
#theme(legend.position="bottom",
#      axis.text = element_text(size = 12),
      #axis.title = element_text(size = 12,face = "bold"),
#      plot.title = element_text(size = 14,hjust = 0.5,face = "bold")) +
#labs(x="", y="",title = "MAE of ProteAge")
#facet_wrap(~ group)
#ggview(w = 3,h = 4.5)
#ggsave("result/bioage/correlation/metric/proteage.pdf",w = 3,h = 4.5, dpi = 300)

#ggplot(MAE[c(2,5),],aes(x = sex,y = Value,fill = BA)) + 
#geom_col(show.legend = F,fill = c("#046586","#28A9A1"),lwd = 0.7) + 
#theme_classic(base_size = 14) + 
#theme(legend.position="bottom",
#      axis.text = element_text(size = 12),
      #axis.title = element_text(size = 12,face = "bold"),
#      plot.title = element_text(size = 14,hjust = 0.5,face = "bold")) +
#labs(x="", y="",title = "MAE of KDM")
#facet_wrap(~ group)
#ggview(w = 3,h = 4.5)
#ggsave("result/bioage/correlation/metric/kdm.pdf",w = 3,h = 4.5, dpi = 300)

#ggplot(MAE[c(3,6),],aes(x = sex,y = Value,fill = BA)) + 
#geom_col(show.legend = F,fill = c("#C9A77C","#F6BBC6"),lwd = 0.7) + 
#theme_classic(base_size = 14) + 
#theme(legend.position="bottom",
#      axis.text = element_text(size = 12),
      #axis.title = element_text(size = 12,face = "bold"),
#      plot.title = element_text(size = 14,hjust = 0.5,face = "bold")) +
#labs(x="", y="",title = "MAE of PhenoAge")
#facet_wrap(~ group)
#ggview(w = 3,h = 4.5)
#ggsave("result/bioage/correlation/metric/phenoage.pdf",w = 3,h = 4.5, dpi = 300)


# female proteage #992224
# female KDM #E3625D
# female phenoage #EF8B67
# male proteage #3E4F94
# male KDM #3E90BF
# male phenoage #A6C0E3

##################################

#### predicted age res
d_f <- melt(ukb_bioage_f[,c(1,2,7:9)],id.vars = c("eid", "age"))

ggplot(d_f,aes(x = age,y = value)) +
stat_smooth(data = d_f,aes(color = variable),method = "loess",size = 1)  +
scale_color_manual(values = c(`ProteAge res` = "#992224",
                              `KDM res` = "#E3625D",
                              `PhenoAge res` = "#EF8B67")) +
#sm_statCorr(color = "#030303", corr_method = "pearson") +
theme_classic(base_size = 14) + 
theme(legend.title = element_blank(),
      legend.position="bottom",
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 18,face = "bold"),
      plot.title = element_text(size = 22,hjust = 0.5,face = "bold"),
      legend.text = element_text(size = 16)) +
labs(x="Chronological age", y="Predicted age res",title = "Female Cohort") 
ggview(w = 6.5, h = 5)
ggsave("1st submission NC/Figures JZH 20241226/figure 2/Figure 2I.pdf",w = 6.5, h = 5,dpi = 300)

d_m <- melt(ukb_bioage_m[,c(1,2,7:9)],id.vars = c("eid", "age"))

ggplot(d_m,aes(x = age,y = value)) +
stat_smooth(data = d_m,aes(color = variable),method = "loess",size = 1)  +
scale_color_manual(values = c(`ProteAge res` = "#3E4F94",
                              `KDM res` = "#3E90BF",
                              `PhenoAge res` = "#A6C0E3")) +
#sm_statCorr(color = "#030303", corr_method = "pearson") +
theme_classic(base_size = 14) + 
theme(legend.title = element_blank(),
      legend.position="bottom",
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 18,face = "bold"),
      plot.title = element_text(size = 22,hjust = 0.5,face = "bold"),
      legend.text = element_text(size = 16)) +
labs(x="Chronological age", y="Predicted age res",title = "Male Cohort") 
ggview(w = 6.5, h = 5)
ggsave("1st submission NC/Figures JZH 20241226/figure 2/Figure 2J.pdf",w = 6.5, h = 5,dpi = 300)


#### predicted age
#d_f <- melt(ukb_bioage_f[,c(1,2,5,6)],id.vars = c("eid", "age"))  

#ggplot(d_f,aes(x = age, y = value)) +
#geom_point(aes(color = variable),size = 1,alpha = 0.5) +
#geom_smooth(aes(color = variable),method = "lm") +
#stat_cor(aes(color = variable), show.legend = F) +
#theme_classic(base_size = 14) + 
#theme(legend.title = element_blank(),
#      legend.position="bottom",
#      axis.text = element_text(size = 14),
#      axis.title = element_text(size = 18,face = "bold"),
#      plot.title = element_text(size = 22,hjust = 0.5,face = "bold")) +
#scale_color_manual(values = c(KDM = "#E3625D",
#                              PhenoAge = "#3E90BF")) +
#xlim(c(40,75)) +
#ylim(c(40,75)) +
#labs(x="Chronological age", y="Predicted age",title = "Female Cohort") 
#ggview(w = 4.5, h = 5.5)
#ggsave("result/bioage/correlation/kdm_phenoage_age_female.pdf",w = 5.5, h = 5.5, dpi = 300)
   
#d_m <- melt(ukb_bioage_m[,c(1,2,5,6)],id.vars = c("eid", "age"))

#ggplot(d_m,aes(x = age, y = value)) +
#geom_point(aes(color = variable),size = 1,alpha = 0.5) +
#geom_smooth(aes(color = variable),method = "lm") +
#stat_cor(aes(color = variable), show.legend = F) +
#theme_classic(base_size = 14) + 
#theme(legend.title = element_blank(),
#      legend.position="bottom",
#      axis.text = element_text(size = 8),
#      axis.title = element_text(size = 12,face = "bold"),
#      plot.title = element_text(size = 16,hjust = 0.5,face = "bold")) +
#scale_color_manual(values = c(KDM = "#EF8B67",
#                              PhenoAge = "#A6C0E3")) +
#labs(x="Chronological age", y="Predicted age",title = "Male Cohort") 
#ggview(w = 4.5, h = 5.5)
#ggsave("result/bioage/correlation/kdm_phenoage_age_male.pdf",w = 4.5, h = 5.5, dpi = 300)

# female proteage #992224
# female KDM #E3625D
# female phenoage #EF8B67
# male proteage #3E4F94
# male KDM #3E90BF
# male phenoage #A6C0E3
#### predicted age
ukb_bioage0 <- rbind(ukb_bioage_f,ukb_bioage_m)

ggplot(ukb_bioage0,aes(x = age, y = KDM)) +
geom_point(aes(color = gender),size = 1,alpha = 0.5) +
geom_smooth(aes(color = gender),method = "lm") +
stat_cor(aes(color = gender,size = 20), show.legend = F,size = 5,digits = 3) +
theme_classic(base_size = 14) + 
theme(legend.title = element_blank(),
      legend.position="bottom",
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 18,face = "bold"),
      plot.title = element_text(size = 22,hjust = 0.5,face = "bold")) +
scale_color_manual(values = c(Female = "#E3625D",
                              Male = "#3E90BF")) +
xlim(c(40,75)) +
ylim(c(40,75)) +
labs(x="Chronological age", y="KDM",title = "KDM Biological Age") 
ggview(w = 5.5, h = 6.5)
ggsave("1st submission NC/Figures JZH 20241226/figure 2/Figure 2C.pdf",w = 5.5, h = 6.5, dpi = 300)

ggplot(ukb_bioage0,aes(x = age, y = PhenoAge)) +
geom_point(aes(color = gender),size = 1,alpha = 0.5) +
geom_smooth(aes(color = gender),method = "lm") +
stat_cor(aes(color = gender,size = 20), show.legend = F,size = 5,digits = 3) +
theme_classic(base_size = 14) + 
theme(legend.title = element_blank(),
      legend.position="bottom",
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 18,face = "bold"),
      plot.title = element_text(size = 22,hjust = 0.5,face = "bold")) +
scale_color_manual(values = c(Female = "#EF8B67",
                              Male = "#A6C0E3")) +
xlim(c(40,75)) +
ylim(c(40,75)) +
labs(x="Chronological age", y="PhenoAge",title = "PhenoAge") 
ggview(w = 5.5, h = 6.5)
ggsave("1st submission NC/Figures JZH 20241226/figure 2/Figure 2D.pdf",w = 5.5, h = 6.5, dpi = 300)


########################  whole cohort
################################ compare BioAge with proteomics aging clock #################################
rm(list=ls())
ukb_bioage <- read.csv("result/bioage/UKB_BioAge_all_features.csv")
ukb_bioage$X <- NULL

########### function
get_BA_resids <- function(ukb_bioage,BA){
  data = ukb_bioage %>% drop_na(BA)
  # Basic model = regress on age alone
  model <- parse(text = sprintf("lm(%s ~ ns(age, df = 3), data = data)", BA)) %>% eval()
  model_predict <- ggpredict(model, terms = c("age"))
  data[,"BA_res"] <- NA
  data[!is.na(data[BA]),"BA_res"] <- resid(model)
  return(residuals(model))
}

################################## female cohort ##################################
p_f <- read.csv("result/ProteAge/female_proteage_whole_female_cohort.csv")
p_f$X <- NULL
ukb_bioage_f <- merge(p_f,ukb_bioage,by = "eid") %>% dplyr::rename(age = age.x) ## 5740 obs, 35 variables
ukb_bioage_f$age.y <- NULL

#}
ukb_bioage_f <- ukb_bioage_f %>% dplyr::select(eid,age,gender,proteage,kdm_original,phenoage_original) %>%
                    dplyr::rename(ProteAge = proteage,
                           KDM = kdm_original,
                           PhenoAge = phenoage_original)
set.seed(777)
for(BA in c("ProteAge","KDM","PhenoAge")){
  BA_res <- paste0(BA, " res")
  ukb_bioage_f[,BA_res] = NA
  ukb_bioage_f[!is.na(ukb_bioage_f[BA]),BA_res] <- get_BA_resids(ukb_bioage_f,BA)
}
write.csv(ukb_bioage_f,"result/bioage/female_cohort_whole.csv")

################################## male cohort ##################################
p_m <- read.csv("result/ProteAge/male_proteage_whole_male_cohort.csv")
p_m$X <- NULL
ukb_bioage_m <- merge(p_m,ukb_bioage,by = "eid") %>% dplyr::rename(age = age.x) ## 4921 obs, 35 variables
ukb_bioage_m$age.y <- NULL ## 4921 obs, 35 variables

ukb_bioage_m <- ukb_bioage_m %>% dplyr::select(eid,age,gender,proteage,kdm_original,phenoage_original) %>%
                    dplyr::rename(ProteAge = proteage,
                           KDM = kdm_original,
                           PhenoAge = phenoage_original)
set.seed(888)
for(BA in c("ProteAge","KDM","PhenoAge")){
  BA_res <- paste0(BA, " res")
  ukb_bioage_m[,BA_res] = NA
  ukb_bioage_m[!is.na(ukb_bioage_m[BA]),BA_res] <- get_BA_resids(ukb_bioage_m,BA)
}
write.csv(ukb_bioage_m,"result/bioage/male_cohort_whole.csv")
#show_col(pal_npg("nrc")(10))

# female proteage #992224
# female KDM #E3625D
# female phenoage #F0C284
# male proteage #3E4F94
# male KDM #3E90BF
# male phenoage #58B6E9

ggplot(ukb_bioage_m,aes(x = age,y = ProteAge)) + 
geom_point(color = "#3E4F94",size = 1,alpha = 0.8) + 
sm_statCorr(color = "#030303", corr_method = "pearson",text_size = 5.5) +
theme_classic(base_size = 14) + 
theme(legend.title = element_blank(),
      legend.position="bottom",
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 18,face = "bold"),
      plot.title = element_text(size = 22,hjust = 0.5,face = "bold")) +
xlim(c(40,75)) +
ylim(c(40,75)) +
labs(x="Chronological Age", y="ProteAge",title = "Male cohort") 
ggview(w = 5.5, h = 5.5)
ggsave("result/bioage/proteAge_age_male.pdf",w = 5.5, h = 5.5,dpi = 300)      
 


## proteage and kdm and phenoage
ggplot(ukb_bioage_m,aes(x = ProteAge,y = KDM)) + 
geom_point(aes(colour = age),size = 1,alpha = 0.5) +
sm_statCorr(color = "#030303", corr_method = "pearson",text_size = 5.5) +
theme_classic(base_size = 14) + 
theme(legend.position="bottom",
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 18,face = "bold"),
      plot.title = element_text(size = 22,hjust = 0.5,face = "bold")) +
scale_color_distiller(palette = "RdBu") + 
xlim(c(40,75)) +
ylim(c(40,75)) +
labs(x="ProteAge", y="KDM",color = "Chronological Age",title = "Male Cohort") 
ggview(w = 5.5, h = 6.7)
ggsave("result/bioage/proteAge_kdm_male.pdf",w = 5.5, h = 6.7,dpi = 300)  

ggplot(ukb_bioage_m,aes(x = ProteAge,y = PhenoAge)) + 
geom_point(aes(colour = age),size = 1,alpha = 0.5) +
sm_statCorr(color = "#030303", corr_method = "pearson",text_size = 5.5) +
theme_classic(base_size = 14) + 
theme(legend.position="bottom",
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 18,face = "bold"),
      plot.title = element_text(size = 22,hjust = 0.5,face = "bold")) +
scale_color_distiller(palette = "RdBu") + 
xlim(c(40,75)) +
ylim(c(40,75)) +
labs(x="ProteAge", y="PhenoAge",color = "Chronological Age",title = "Male Cohort") 
ggview(w = 5.5, h = 6.7)
ggsave("result/bioage/proteAge_phenoage_male.pdf",w = 5.5, h = 6.7,dpi = 300) 


############################# evaluation metrics
#填补缺失值
ukb_bioage_f$KDM <- impute(ukb_bioage_f$KDM, mean)
ukb_bioage_f$PhenoAge <- impute(ukb_bioage_f$PhenoAge, mean)
ukb_bioage_m$KDM <- impute(ukb_bioage_m$KDM, mean)
ukb_bioage_m$PhenoAge <- impute(ukb_bioage_m$PhenoAge, mean)

### MAE
## male
mae_female_proage <- mae(ukb_bioage_f$age,ukb_bioage_f$ProteAge) ## 2.237385
mae_female_kdm <- mae(ukb_bioage_f$age,ukb_bioage_f$KDM) ## 2.906301
mae_female_phenoage <- mae(ukb_bioage_f$age,ukb_bioage_f$PhenoAge)  ## 7.256684
## female
mae_male_proage <- mae(ukb_bioage_m$age,ukb_bioage_m$ProteAge) ## 2.445928
mae_male_kdm <- mae(ukb_bioage_m$age,ukb_bioage_m$KDM) ## 3.987592
mae_male_phenoage <- mae(ukb_bioage_m$age,ukb_bioage_m$PhenoAge) ## 6.082424


MAE <- data.frame(mae_female_proage,mae_female_kdm,mae_female_phenoage,
                  mae_male_proage,mae_male_kdm,mae_male_phenoage)
colnames(MAE) <- c("Female ProteAge","Female KDM","Female PhenoAge",
                   "Male ProteAge","Male KDM","Male PhenoAge")
MAE <- t(MAE) %>% as.data.frame
colnames(MAE) <- c("Value")
MAE$BA <- rownames(MAE)
MAE$group <- c("ProteAge","KDM","PhenoAge","ProteAge","KDM","PhenoAge")
MAE$sex <- c("Female","Female","Female","Male","Male","Male")
MAE$BA <- factor(MAE$BA,
                 levels = c("Female ProteAge","Male ProteAge","Female KDM",
                            "Male KDM","Female PhenoAge","Male PhenoAge"))
MAE$group <- factor(MAE$group,levels = c("ProteAge","KDM","PhenoAge"))

#show_col(pal_npg("nrc")(10))

# female proteage #992224
# female KDM #E3625D
# female phenoage #F0C284
# male proteage #3E4F94
# male KDM #3E90BF
# male phenoage #58B6E9

ggplot(MAE,aes(x = BA,y = Value)) + 
geom_bar(aes(fill = BA),stat = "identity",show.legend = F,lwd = 0.7) + 
theme_classic(base_size = 14) + 
scale_fill_npg() + 
theme(axis.text.y = element_text(size = 15),
      plot.title = element_text(size = 30,hjust = 0.5,face = "bold"),
      axis.text.x = element_text(angle = 15,vjust = 0.6,size = 14.2)) +
scale_fill_manual(values = c(`Female ProteAge` = "#992224",
                              `Male ProteAge` = "#3E4F94",
                              `Female KDM` = "#E3625D",
                              `Male KDM` = "#3E90BF",
                              `Female PhenoAge` = "#F0C284",
                              `Male PhenoAge` = "#58B6E9")) +
labs(x="", y="",title = "Mean Absolute Error") 
ggview(w = 7,h = 6)
ggsave("result/bioage/proteage_kdm_phenoage_MAE.png",w = 7,h = 6, dpi = 300)


#ggplot(MAE[c(1,4),],aes(x = sex,y = Value,fill = BA)) + 
#geom_col(show.legend = F,fill = c("#E71F19","#F4A016"),lwd = 0.7) + 
#theme_classic(base_size = 14) + 
#theme(legend.position="bottom",
#      axis.text = element_text(size = 12),
      #axis.title = element_text(size = 12,face = "bold"),
#      plot.title = element_text(size = 14,hjust = 0.5,face = "bold")) +
#labs(x="", y="",title = "MAE of ProteAge")
#facet_wrap(~ group)
#ggview(w = 3,h = 4.5)
#ggsave("result/bioage/correlation/metric/proteage.pdf",w = 3,h = 4.5, dpi = 300)

#ggplot(MAE[c(2,5),],aes(x = sex,y = Value,fill = BA)) + 
#geom_col(show.legend = F,fill = c("#046586","#28A9A1"),lwd = 0.7) + 
#theme_classic(base_size = 14) + 
#theme(legend.position="bottom",
#      axis.text = element_text(size = 12),
      #axis.title = element_text(size = 12,face = "bold"),
#      plot.title = element_text(size = 14,hjust = 0.5,face = "bold")) +
#labs(x="", y="",title = "MAE of KDM")
#facet_wrap(~ group)
#ggview(w = 3,h = 4.5)
#ggsave("result/bioage/correlation/metric/kdm.pdf",w = 3,h = 4.5, dpi = 300)

#ggplot(MAE[c(3,6),],aes(x = sex,y = Value,fill = BA)) + 
#geom_col(show.legend = F,fill = c("#C9A77C","#F6BBC6"),lwd = 0.7) + 
#theme_classic(base_size = 14) + 
#theme(legend.position="bottom",
#      axis.text = element_text(size = 12),
      #axis.title = element_text(size = 12,face = "bold"),
#      plot.title = element_text(size = 14,hjust = 0.5,face = "bold")) +
#labs(x="", y="",title = "MAE of PhenoAge")
#facet_wrap(~ group)
#ggview(w = 3,h = 4.5)
#ggsave("result/bioage/correlation/metric/phenoage.pdf",w = 3,h = 4.5, dpi = 300)


# female proteage #992224
# female KDM #E3625D
# female phenoage #EF8B67
# male proteage #3E4F94
# male KDM #3E90BF
# male phenoage #A6C0E3

##################################

#### predicted age res
d_f <- melt(ukb_bioage_f[,c(1,2,7:9)],id.vars = c("eid", "age"))

ggplot(d_f,aes(x = age,y = value)) +
stat_smooth(data = d_f,aes(color = variable),method = "loess",size = 1)  +
scale_color_manual(values = c(`ProteAge res` = "#992224",
                              `KDM res` = "#E3625D",
                              `PhenoAge res` = "#EF8B67")) +
#sm_statCorr(color = "#030303", corr_method = "pearson") +
theme_classic(base_size = 14) + 
theme(legend.title = element_blank(),
      legend.position="bottom",
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 18,face = "bold"),
      plot.title = element_text(size = 22,hjust = 0.5,face = "bold"),
      legend.text = element_text(size = 16)) +
labs(x="Chronological age", y="Predicted age res",title = "Female Cohort") 
ggview(w = 6.5, h = 5)
ggsave("result/bioage/proteage_kdm_phenoage_female.pdf",w = 6.5, h = 5,dpi = 300)

d_m <- melt(ukb_bioage_m[,c(1,2,7:9)],id.vars = c("eid", "age"))

ggplot(d_m,aes(x = age,y = value)) +
stat_smooth(data = d_m,aes(color = variable),method = "loess",size = 1)  +
scale_color_manual(values = c(`ProteAge res` = "#3E4F94",
                              `KDM res` = "#3E90BF",
                              `PhenoAge res` = "#A6C0E3")) +
#sm_statCorr(color = "#030303", corr_method = "pearson") +
theme_classic(base_size = 14) + 
theme(legend.title = element_blank(),
      legend.position="bottom",
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 18,face = "bold"),
      plot.title = element_text(size = 22,hjust = 0.5,face = "bold"),
      legend.text = element_text(size = 16)) +
labs(x="Chronological age", y="Predicted age res",title = "Male Cohort") 
ggview(w = 6.5, h = 5)
ggsave("result/bioage/proteage_kdm_phenoage_male.pdf",w = 6.5, h = 5,dpi = 300)


#### predicted age
#d_f <- melt(ukb_bioage_f[,c(1,2,5,6)],id.vars = c("eid", "age"))  

#ggplot(d_f,aes(x = age, y = value)) +
#geom_point(aes(color = variable),size = 1,alpha = 0.5) +
#geom_smooth(aes(color = variable),method = "lm") +
#stat_cor(aes(color = variable), show.legend = F) +
#theme_classic(base_size = 14) + 
#theme(legend.title = element_blank(),
#      legend.position="bottom",
#      axis.text = element_text(size = 14),
#      axis.title = element_text(size = 18,face = "bold"),
#      plot.title = element_text(size = 22,hjust = 0.5,face = "bold")) +
#scale_color_manual(values = c(KDM = "#E3625D",
#                              PhenoAge = "#3E90BF")) +
#xlim(c(40,75)) +
#ylim(c(40,75)) +
#labs(x="Chronological age", y="Predicted age",title = "Female Cohort") 
#ggview(w = 4.5, h = 5.5)
#ggsave("result/bioage/correlation/kdm_phenoage_age_female.pdf",w = 5.5, h = 5.5, dpi = 300)
   
#d_m <- melt(ukb_bioage_m[,c(1,2,5,6)],id.vars = c("eid", "age"))

#ggplot(d_m,aes(x = age, y = value)) +
#geom_point(aes(color = variable),size = 1,alpha = 0.5) +
#geom_smooth(aes(color = variable),method = "lm") +
#stat_cor(aes(color = variable), show.legend = F) +
#theme_classic(base_size = 14) + 
#theme(legend.title = element_blank(),
#      legend.position="bottom",
#      axis.text = element_text(size = 8),
#      axis.title = element_text(size = 12,face = "bold"),
#      plot.title = element_text(size = 16,hjust = 0.5,face = "bold")) +
#scale_color_manual(values = c(KDM = "#EF8B67",
#                              PhenoAge = "#A6C0E3")) +
#labs(x="Chronological age", y="Predicted age",title = "Male Cohort") 
#ggview(w = 4.5, h = 5.5)
#ggsave("result/bioage/correlation/kdm_phenoage_age_male.pdf",w = 4.5, h = 5.5, dpi = 300)

# female proteage #992224
# female KDM #E3625D
# female phenoage #EF8B67
# male proteage #3E4F94
# male KDM #3E90BF
# male phenoage #A6C0E3
#### predicted age
ukb_bioage0 <- rbind(ukb_bioage_f,ukb_bioage_m)

ggplot(ukb_bioage0,aes(x = age, y = KDM)) +
geom_point(aes(color = gender),size = 1,alpha = 0.5) +
geom_smooth(aes(color = gender),method = "lm") +
stat_cor(aes(color = gender,size = 20), show.legend = F,size = 5) +
theme_classic(base_size = 14) + 
theme(legend.title = element_blank(),
      legend.position="bottom",
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 18,face = "bold"),
      plot.title = element_text(size = 22,hjust = 0.5,face = "bold")) +
scale_color_manual(values = c(Female = "#E3625D",
                              Male = "#3E90BF")) +
xlim(c(40,75)) +
ylim(c(40,75)) +
labs(x="Chronological age", y="KDM",title = "KDM Biological Age") 
ggview(w = 5.5, h = 6.5)
ggsave("result/bioage/kdm_age_female_male.pdf",w = 5.5, h = 6.5, dpi = 300)

ggplot(ukb_bioage0,aes(x = age, y = PhenoAge)) +
geom_point(aes(color = gender),size = 1,alpha = 0.5) +
geom_smooth(aes(color = gender),method = "lm") +
stat_cor(aes(color = gender,size = 20), show.legend = F,size = 5) +
theme_classic(base_size = 14) + 
theme(legend.title = element_blank(),
      legend.position="bottom",
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 18,face = "bold"),
      plot.title = element_text(size = 22,hjust = 0.5,face = "bold")) +
scale_color_manual(values = c(Female = "#EF8B67",
                              Male = "#A6C0E3")) +
xlim(c(40,75)) +
ylim(c(40,75)) +
labs(x="Chronological age", y="PhenoAge",title = "PhenoAge") 
ggview(w = 5.5, h = 6.5)
ggsave("result/bioage/phenoage_age_female_male.pdf",w = 5.5, h = 6.5, dpi = 300)
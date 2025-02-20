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
library(cardx)

################
yrs <- read.csv("data/protein_cohort_follow_up_information.csv")
colnames(yrs) <- c("eid","Date_lost_to_follow_up","Date_of_attending_assessment_centre",
                 "Date_of_death","Age_when_attended_assessment_centre","Age_at_recruitment")
yrs$'Death record' <- ifelse(yrs$Date_of_death == "","No","Yes")
yrs$'Death record' <- ifelse(!yrs$Date_lost_to_follow_up == "","Lost to follow up",yrs$'Death record')
yrs <- yrs[,c(1,7)]

age_sex <- read.csv("data/protein_cohort_sex_age.csv") %>% dplyr::rename(sex = p31,age = p21022)
age_sex$X <- NULL


#p <- read.csv("data/protein_data.csv")[,c(1,2)]

d <- read.csv("data/disease_info.csv")
colnames(d) <- c("eid","age","sex","Diabetes","Vascular/heart problems","Cancer")
table(d$`Vascular/heart problems`)
d[d == "Prefer not to answer"] <- "Do not know"
d[d == ""] <- "Do not know"

d$`High blood pressure` <- ifelse(d$`Vascular/heart problems` == "Angina|High blood pressure"|
                                    d$`Vascular/heart problems` == "Angina|Stroke|High blood pressure"|
                                    d$`Vascular/heart problems` == "Heart attack|Angina|High blood pressure"|
                                    d$`Vascular/heart problems` == "Heart attack|Angina|Stroke|High blood pressure"|
                                    d$`Vascular/heart problems` == "Heart attack|High blood pressure"|
                                    d$`Vascular/heart problems` == "Heart attack|Stroke|High blood pressure"|
                                    d$`Vascular/heart problems` == "High blood pressure"|
                                    d$`Vascular/heart problems` == "Stroke|High blood pressure",
                                  "Yes","No")
d$`High blood pressure` <- ifelse(d$`Vascular/heart problems` == "Do not know","Do not know",d$`High blood pressure`)
table(d$`High blood pressure`)
d$`Heart attack` <- ifelse(d$`Vascular/heart problems` == "Heart attack"|
                             d$`Vascular/heart problems` == "Heart attack|Angina"|
                             d$`Vascular/heart problems` == "Heart attack|Angina|High blood pressure"|
                             d$`Vascular/heart problems` == "Heart attack|Angina|Stroke"|
                             d$`Vascular/heart problems` == "Heart attack|Angina|Stroke|High blood pressure"|
                             d$`Vascular/heart problems` == "Heart attack|High blood pressure"|
                             d$`Vascular/heart problems` == "Heart attack|Stroke"|
                             d$`Vascular/heart problems` == "Heart attack|Stroke|High blood pressure",
                             "Yes","No")
d$`Heart attack` <- ifelse(d$`Vascular/heart problems` == "Do not know","Do not know",d$`Heart attack`)
d$Stroke <- ifelse(d$`Vascular/heart problems` == "Angina|Stroke"|
                     d$`Vascular/heart problems` == "Angina|Stroke|High blood pressure"|
                     d$`Vascular/heart problems` == "Heart attack|Angina|Stroke"|
                     d$`Vascular/heart problems` == "Heart attack|Angina|Stroke|High blood pressure"|
                     d$`Vascular/heart problems` == "Heart attack|Stroke"|
                     d$`Vascular/heart problems` == "Heart attack|Stroke|High blood pressure"|
                     d$`Vascular/heart problems` == "Stroke"|
                     d$`Vascular/heart problems` == "Stroke|High blood pressure","Yes","No")
d$Stroke <- ifelse(d$`Vascular/heart problems` == "Do not know","Do not know",d$Stroke)
d$Cancer <- ifelse(d$Cancer == "Yes - you will be asked about this later by an interviewer","Yes",d$Cancer)


cli <- read.csv("data/protein_cohort_clinical_table.csv")
cli <- cli %>% rename(eid = `Participant.ID`,
                      `Townsend deprivation index` = `Townsend.deprivation.index.at.recruitment`,
                      `Smoking Status` = `Smoking.status...Instance.0`,
                      `Alcohol Drinker Status` = `Alcohol.drinker.status...Instance.0`,
                      BMI = `Body.mass.index..BMI....Instance.0`,
                      Ethnic = `Ethnic.background...Instance.0`)
cli$Ethnic <- ifelse(cli$Ethnic == "","Do not know",cli$Ethnic)
cli$`Smoking Status` <- ifelse(cli$`Smoking Status` == "","Do not know",cli$`Smoking Status`)
cli$`Alcohol Drinker Status` <- ifelse(cli$`Alcohol Drinker Status` == "","Do not know",cli$`Alcohol Drinker Status`)

cli[cli == "British"] <- "White"
cli[cli == "Irish"] <- "White"
cli[cli == "Any other white backgroud"] <- "White"

cli[cli == "White and Black Caribbean"] <- "Mixed"
cli[cli == "White and Black African"] <- "Mixed"
cli[cli == "White and Asian"] <- "Mixed"
cli[cli == "Any other mixed background"] <- "Mixed"

cli[cli == "Indian"] <- "Asian or Asian British"
cli[cli == "Pakistani"] <- "Asian or Asian British"
cli[cli == "Bangladeshi"] <- "Asian or Asian British"
cli[cli == "Any other Asian background"] <- "Asian or Asian British"


cli[cli == "Caribbean"] <- "Black or Black British"
cli[cli == "African"] <- "Black or Black British"
cli[cli == "Any other Black background"] <- "Black or Black British"

#f0 <- merge(age_sex,cli,by = "eid")
#f01 <- merge(f0,d,by = "eid")
#f001 <- merge(f01,yrs,by = "eid")
#f1 <- merge(f001,p,by = "eid")
f0 <- merge(d,cli,by = "eid")
f01 <- merge(f0,yrs,by = "eid")


f01$`Age group` <- ifelse(f01$age >= 39 & f01$age <= 50,"39-50 years old","yeas")
f01$`Age group` <- ifelse(f01$age >= 51 & f01$age <= 60,"51-60 years old",f01$`Age group`)
f01$`Age group` <- ifelse(f01$age >= 61 & f01$age <= 70,"61-70 years old",f01$`Age group`)
f01 <- f01 %>% dplyr::rename(Sex = sex)

table3 <- f01 %>%
  select(`Age group`,`Townsend deprivation index`,Sex,BMI,`Smoking Status`,`Alcohol Drinker Status`,Diabetes,`High blood pressure`,`Heart attack`,Stroke,Cancer,'Death record') %>%
  tbl_summary(by = Sex,
              missing = "no",
              statistic = all_continuous() ~ c("{mean} ({sd})",
                                     "{median} ({p25}, {p75})", 
                                     "{min}, {max}"),
              type = list(`Townsend deprivation index` = "continuous2",BMI = "continuous2")) %>%
  #tbl_summary(type = list(age ~ "categorical", gender ~ "categorical")) %>%
  bold_labels() %>%
  add_p() %>%
  modify_header(label = "")
  

table3 %>%
  as_hux_xlsx(file = '1st submission NC/Table JZH 20241227/UKB_PPP_clinical_table JZH 20241227.xlsx')

table3 %>%
  as_flex_table() %>%
  flextable::save_as_docx(table3, path = '1st submission NC/Table JZH 20241227/UKB_PPP_clinical_table JZH 20241227.docx')


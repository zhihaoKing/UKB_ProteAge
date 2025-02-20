setwd("/Users/zhihaojin/Documents/neuroscience_PhD_in_FDU/PhD_project/bioinformatics/UKBiobank_version_2")

rm(list=ls())
library(glmnet)
library(psych)
library(tidyr)
library(Metrics)
library(ggplot2)
library(doMC)
library(ggpubr)
registerDoMC(cores = 4)

#################  whole cohort (female) ##############
rm(list=ls());gc()
#f <- read.csv("data/healthy_cohort_sex_age.csv")
#colnames(f) <- c("EID","Sex","Age")
f <- read.csv("data/protein_cohort_sex_age.csv")
colnames(f) <- c("eid","sex","age")
p <- read.csv("data/protein_data.csv")
colnames(p)[-1] <- toupper(colnames(p)[-1])


data_p <- merge(f,p,by = "eid")
rownames(data_p) <- data_p$eid
data_p$eid <- NULL
data_p_f <- subset(data_p,sex == "Female")

# select protein padj<0.05
### age_protein_cor_fe <- readRDS("result/correlation/healthy_cohort/protein/female/age_protein_cor_fe.rds") ## padj < 0.05
### data_p <- data_p[,c("Sex","Age",age_protein_cor_fe$Protein)]

# select protein |cor|>0.3
#Coef <- read.csv("result/cv_glmnet/healthy_cohort/protein(padj<0.05)/female_coef_alpha0.65.csv")
#colnames(Coef) <- c("Protein","coef")
#Coef <- Coef[order(Coef[,2],decreasing = T),]
#Coef <- Coef[-1,]
#Coef <- subset(Coef,coef < 0|coef > 0)
#protein <- merge(Coef,age_protein_cor_fe,by.x = "Protein",by.y = "Protein")
### age_protein_cor_fe <- readRDS("result/correlation/healthy_cohort/protein/female/age_protein_cor_fe.rds") 
### age_protein_cor_fe <- subset(age_protein_cor_fe,correlation < -0.3|correlation > 0.3)
### data_p <- data_p[,c("Sex","Age",age_protein_cor_fe$Protein)]

#############
for(i in 3:ncol(data_p_f)) {
  data_p_f[,i][is.na(data_p_f[,i])] <- mean(data_p_f[,i], na.rm=TRUE)
}

set.seed(5)
train_f <- sample(nrow(data_p_f),0.7*nrow(data_p_f))
data_train_f <- data_p_f[train_f,]    # n = 20005
data_validation_f	<- data_p_f[-train_f,]    # n = 8574

set.seed(123)
models_f <- list()
for (i in 0:20) {
  name <- paste0("alpha", i/20)
  models_f[[name]] <- cv.glmnet(x = as.matrix(data_train_f[,-c(1:2)]), 
                              y = as.matrix(data_train_f$age), 
                              type.measure = "mse", 
                              alpha = i/20, 
                              family = "gaussian",
                              nfold = 20,
                              parallel = TRUE,
                              trace.it = TRUE)
}
#save(models, file = "result/ProteAge/female_model.RData")
#models_f <- models
save(models_f, file = "result/ProteAge/female_model.RData")

load("result/ProteAge/female_model.RData")
results_f <- data.frame()
for (i in 0:20) {
  name <- paste0("alpha", i/20)
  ## Use each model to predict 'y' given the Testing dataset
  predicted <- predict(models_f[[name]], 
                       s = models_f[[name]]$lambda.min, 
                       newx = as.matrix(data_validation_f[,-c(1:2)]))
  ## Calculate the Mean Squared Error...
  mse <- mean((data_validation_f$age - predicted)^2)
  ## Store the results
  temp <- data.frame(alpha=i/20, mse=mse, name=name)
  results_f <- rbind(results_f, temp)
}
print(results_f) ## alpha 0.05 have the smallest mse value
write.csv(results_f,"result/ProteAge/female_MSE_of_different_alpha_all_protein.csv")


ggplot(results_f, aes(x = alpha,y = mse)) + 
geom_point() +
ggtitle("MSE of different alpha",subtitle = "Female cohort (all proteins)") +
xlab("Alpha value") + ylab("Mean Square Error (MSE)") +
theme(panel.background = element_rect(fill = "white", colour = "gray1"))
ggsave("result/ProteAge/female_MSE_of_different_alpha_all_protein.png",w = 3.5,h = 3.5,dpi = 300)

pdf("result/ProteAge/female_MSE_of_different_lambda.pdf",w = 6,h = 6)
plot(models_f[["alpha0.05"]])
dev.off()

coef_min_f <- coef(models_f[["alpha0.05"]], s = "lambda.min")[1:2924,] %>% as.data.frame
coef_min_f$protein <- rownames(coef_min_f)
colnames(coef_min_f) <- c("coefficient","protein")
write.csv(coef_min_f[-1,],file = "result/ProteAge/female_coef_alpha0.05.csv")

predicted_age_f <- predict(models_f[["alpha0.05"]], 
                         s = models_f[["alpha0.05"]]$lambda.min, 
                         newx = as.matrix(data_validation_f[,-c(1:2)]))
s <- data.frame(predicted_age_f,data_validation_f$age)
colnames(s) <- c("proteage","age")
s$eid <- rownames(s)
write.csv(s,"result/ProteAge/female_proteage_validation_cohort.csv")

predicted_age_whole_cohort <- predict(models_f[["alpha0.05"]], 
                         s = models_f[["alpha0.05"]]$lambda.min, 
                         newx = as.matrix(data_p_f[,-c(1:2)]))
s <- data.frame(predicted_age_whole_cohort,data_p_f$age)                         
colnames(s) <- c("proteage","age")
s$eid <- rownames(s)
write.csv(s,"result/ProteAge/female_proteage_whole_female_cohort.csv")

#ggscatter(s,x = "Chronological_Age", y = "Predicted_Age", add = "reg.line", color = "deeppink3",size = 1,
#               conf.int = T,add.params = list(color = "#030303",fill = "ivory3")) +
#     stat_cor(method = "pearson") + ggtitle("Proteomic Aging Clock",subtitle = "Female cohort") + 
#     xlab("Chronological Age") + ylab("Predicted Age") +
#     theme(plot.title = element_text(size = 15)) 
#ggsave("result/scatter_plot/whole_cohort/predicted_age_chronological_age/whole_female_cohort.jpeg",w = 4,h = 4,dpi = 300)     


#################  whole cohort (male) ##############

data_p_m <- subset(data_p,sex == "Male")

# # select protein padj<0.05
### age_protein_cor_ma <- readRDS("result/correlation/healthy_cohort/protein/male/age_protein_cor_ma.rds")
### data_p <- data_p[,c("Sex","Age",age_protein_cor_ma$Protein)]

#Coef <- read.csv("result/cv_glmnet/healthy_cohort/protein(padj<0.05)/male_coef_alpha0.5.csv")
#colnames(Coef) <- c("Protein","coef")
#Coef <- Coef[order(Coef[,2],decreasing = T),]
#Coef <- Coef[-1,]
#Coef <- subset(Coef,coef < 0|coef > 0)

# select protein |cor|>0.3
### age_protein_cor_ma <- readRDS("result/correlation/healthy_cohort/protein/male/age_protein_cor_ma.rds")
### age_protein_cor_ma <- subset(age_protein_cor_ma,correlation < -0.3|correlation > 0.3)
### data_p <- data_p[,c("Sex","Age",age_protein_cor_ma$Protein)]

for(i in 3:ncol(data_p_m)) {
  data_p_m[,i][is.na(data_p_m[,i])] <- mean(data_p_m[,i], na.rm=TRUE)
}

set.seed(5)
train_m <- sample(nrow(data_p_m),0.7*nrow(data_p_m))
data_train_m <- data_p_m[train_m,]    # n = 17103
data_validation_m	<- data_p_m[-train_m,]  # n = 7331

set.seed(123)
models_m <- list()
for (i in 0:20) {
  name <- paste0("alpha", i/20)
  models_m[[name]] <- cv.glmnet(x = as.matrix(data_train_m[,-c(1:2)]), 
                              y = as.matrix(data_train_m$age), 
                              type.measure = "mse", 
                              alpha = i/20, 
                              family = "gaussian",
                              nfold = 20,
                              parallel = TRUE,
                              trace.it = TRUE)
}
save(models_m, file = "result/ProteAge/male_model.RData")

load("result/ProteAge/male_model.RData")
results_m <- data.frame()
for (i in 0:20) {
  name <- paste0("alpha", i/20)
  ## Use each model to predict 'y' given the Testing dataset
  predicted <- predict(models_m[[name]], 
                       s = models_m[[name]]$lambda.min, 
                       newx = as.matrix(data_validation_m[,-c(1:2)]))
  ## Calculate the Mean Squared Error...
  mse <- mean((data_validation_m$age - predicted)^2)
  ## Store the results
  temp <- data.frame(alpha=i/20, mse=mse, name=name)
  results_m <- rbind(results_m, temp)
}
print(results_m) ## alpha1
#plot(results$alpha, results$mse)
write.csv(results_m,"result/ProteAge/male_MSE_of_different_alpha_all_protein.csv")

ggplot(results_m, aes(x = alpha,y = mse)) + 
geom_point() +
ggtitle("MSE of different alpha",subtitle = "Male cohort (all proteins)") +
xlab("Alpha value") + ylab("Mean Square Error (MSE)") +
theme(panel.background = element_rect(fill = "white", colour = "gray1"))
ggsave("result/ProteAge/male_MSE_of_different_alpha_all_protein.pdf",w = 3.5,h = 3.5,dpi = 300)

pdf("result/ProteAge/male_MSE_of_different_lambda.pdf",w = 6,h = 6)
plot(models_m[["alpha1"]])
dev.off()

coef_min_m <- coef(models_m[["alpha1"]], s = "lambda.min")[1:2924,] %>% as.data.frame
coef_min_m$protein <- rownames(coef_min_m)
colnames(coef_min_m) <- c("coefficient","protein")
write.csv(coef_min_m[-1,],file = "result/ProteAge/male_coef_alpha1.csv")

predicted_age_m <- predict(models_m[["alpha1"]], 
                         s = models_m[["alpha1"]]$lambda.min, 
                         newx = as.matrix(data_validation_m[,-c(1:2)]))

s <- data.frame(predicted_age_m,data_validation_m$age)
colnames(s) <- c("proteage","age")
s$eid <- rownames(s)
write.csv(s,"result/ProteAge/male_proteage_validation_cohort.csv")

predicted_age_whole_cohort <- predict(models_m[["alpha1"]], 
                         s = models_m[["alpha1"]]$lambda.1se, 
                         newx = as.matrix(data_p_m[,-c(1:2)]))
s <- data.frame(predicted_age_whole_cohort,data_p_m$age)                         
colnames(s) <- c("proteage","age")
s$eid <- rownames(s)
write.csv(s,"result/ProteAge/male_proteage_whole_cohort_male.csv")
#ggscatter(s,x = "Chronological_Age", y = "Predicted_Age", add = "reg.line", color = "dodgerblue3",size = 1,
#               conf.int = T,add.params = list(color = "#030303",fill = "ivory3")) +
#     stat_cor(method = "pearson") + ggtitle("Proteomic Aging Clock",subtitle = "Male cohort") + 
#     xlab("Chronological Age") + ylab("Predicted Age") +
#     theme(plot.title = element_text(size = 15)) 
#ggsave("result/scatter_plot/whole_cohort/predicted_age_chronological_age/whole_male_cohort.jpeg",w = 4,h = 4,dpi = 300)  




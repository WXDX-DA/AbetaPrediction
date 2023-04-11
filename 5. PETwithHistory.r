
library(tidyverse)
library(foreign)
library(kableExtra)
library(pROC) ## ROC
library(caret) ## ROC
library(rpart)
library(rpart.plot)
library(gtsummary)
library("scales")  
hex_codes1 <- hue_pal()(6)  

# Data Preparation

ad.raw = read.csv("../data/ad_full_data.csv")%>%
  mutate(PET = as.factor(PET),
         Diagnosis = factor(Diagnosis, levels = c("NC", "SCD", "MCI", "AD", "non-AD")),
         Gender = factor(Gender, levels = c("Male", "Female")),
         APOE_risk = as.factor(APOE_risk))

dim(ad.raw)



# Prediction with/out hypertention

ROC_performance = function(df, Model, ROCFL) {
  
  Prob = predict(Model, newdata = df, type = "prob")
  TH <-
     as.numeric(
       coords(
         ROCFL,
         "best",
         ret = "threshold",
         transpose = FALSE,
         best.method = "youden"
       )[1, 1]
     )
  
  tab = table(Prob[, 2] > TH, df$PET)
  
  
  (PPV = tab[2, 2] / sum(tab[2, ]))
  (NPV = tab[1, 1] / sum(tab[1, ]))
  (TC = sum(diag(tab)) / sum(tab))
  res = matrix(c(TC, PPV, NPV), 3, 1)
  colnames(res) = "Result"
  rownames(res) = c(
    "Accuracy",
    "PPV",
    "NPV"
  )
  return(round(t(res), 3))
}

hypertention = ad.raw %>%
  filter(`高血压` == 1)
non_hypertention = ad.raw %>%
  filter(`高血压` == 0)


get_accuracy = function(df1, df2, df3, Model, Names){
  TC1 = ROC_performance(df1,
                Model$Model,
                Model$ROC) 
  TC2 = ROC_performance(df2,
                Model$Model,
                Model$ROC) 
  TC3 = ROC_performance(df3,
                Model$Model,
                Model$ROC) 
  Model_performance = rbind(TC1, TC2, TC3)
  rownames(Model_performance) = Names
  return(Model_performance)
}
fullData_models = readRDS("../data/fullData_models.RDS")

mutate_df = function(df, Disease, Model){
  result = df %>%
    as.data.frame() %>%
    select(Accuracy) %>%
    mutate(
      Disease = Disease,
      Population = c("All", "Without", "With"),
      Model = Model
    )
  return(result)
}

hyper_full = get_accuracy(
  ad.raw,
  non_hypertention,
  hypertention,
  fullData_models$full,
  c("Complete Data", "Non-Hypertention", "Hypertention")
) 
hyper_full = mutate_df(hyper_full, "Hypertention", "Full Model")
  

hyper_best = get_accuracy(
  ad.raw,
  non_hypertention,
  hypertention,
  fullData_models$Best,
  c("Complete Data", "Non-Hypertention", "Hypertention")
) 

hyper_best = mutate_df(hyper_best, "Hypertention", "Best Model")


hyper_refined = get_accuracy(
  ad.raw,
  non_hypertention,
  hypertention,
  fullData_models$Refined,
  c("Complete Data", "Non-Hypertention", "Hypertention")
) 

hyper_refined = mutate_df(hyper_refined, "Hypertention", "Refined Model")


HLP = ad.raw %>%
  filter(`高脂血症` == 1)
non_HLP = ad.raw %>%
  filter(`高脂血症` == 0)



HLP_full = get_accuracy(
  ad.raw,
  non_HLP,
  HLP,
  fullData_models$full,
  c("Complete Data", "No HLP", "HLP")
) 

HLP_full = mutate_df(HLP_full, "HLP", "Full Model")


HLP_best = get_accuracy(
  ad.raw,
  non_HLP,
  HLP,
  fullData_models$Best,
  c("Complete Data", "No HLP", "HLP")
) 

HLP_best = mutate_df(HLP_best, "HLP", "Best Model")

HLP_refined = get_accuracy(
  ad.raw,
  non_HLP,
  HLP,
  fullData_models$Refined,
  c("Complete Data", "No HLP", "HLP")
) 

HLP_refined = mutate_df(HLP_refined, "HLP", "Refined Model")



# Prediction with/out Diabetes

Diabetes = ad.raw %>%
  filter(`糖尿病` == 1)
non_Diabetes = ad.raw %>%
  filter(`糖尿病` == 0)

Dia_full = get_accuracy(
  ad.raw,
  non_Diabetes,
  Diabetes,
  fullData_models$full,
  c("Complete Data", "No Diabetes", "Diabetes")
) 

Dia_full = mutate_df(Dia_full, "Diabetes", "Full Model")

Dia_best = get_accuracy(
  ad.raw,
  non_Diabetes,
  Diabetes,
  fullData_models$Best,
  c("Complete Data", "No Diabetes", "Diabetes")
) 

Dia_best = mutate_df(Dia_best, "Diabetes", "Best Model")

Dia_refined = get_accuracy(
  ad.raw,
  non_Diabetes,
  Diabetes,
  fullData_models$Refined,
  c("Complete Data", "No Diabetes", "Diabetes")
) 

Dia_refined = mutate_df(Dia_refined, "Diabetes", "Refined Model")




# Prediction with/out AD Family History

AD.fami = ad.raw %>%
  filter(Fami_AD == 1)
non_AD.fami = ad.raw %>%
  filter(Fami_AD == 0)

AD_full = get_accuracy(
  ad.raw,
  non_AD.fami,
  AD.fami,
  fullData_models$full,
  c("Complete Data", "No AD family History", "AD family History")
)

AD_full = mutate_df(AD_full, "AD Family", "Full Model")

AD_best = get_accuracy(
  ad.raw,
  non_AD.fami,
  AD.fami,
  fullData_models$Best,
  c("Complete Data", "No AD family History", "AD family History")
) 

AD_best = mutate_df(AD_best, "AD Family", "Best Model")

AD_refined = get_accuracy(
  ad.raw,
  non_AD.fami,
  AD.fami,
  fullData_models$Refined,
  c("Complete Data", "No AD family History", "AD family History")
) 

AD_refined = mutate_df(AD_refined, "AD Family", "Refined Model")



# Line Plot

Performances = bind_rows(hyper_full, hyper_best, hyper_refined,
                         HLP_full, HLP_best, HLP_refined,
                         Dia_full, Dia_best, Dia_refined,
                         AD_full, AD_best, AD_refined) %>%
  mutate(Model = factor(Model, levels = c("Full Model", "Best Model", "Refined Model")))
Performances

tiff(
  "../output/supplementary figure 2.tiff",
  units = "mm",
  width = 185,
  height = 145,
  res = 500,
  pointsize = 10
)

ggplot(Performances,
       aes(y = Accuracy,
           x = Population,
           color = Disease)) +
  geom_point(size = 2) +
  geom_line(aes(group = Disease)) +
  facet_grid(cols = vars(Model), switch="both") +
  labs(y = "Accuracy") +
  ggtitle("Accuracy of the Established Models Concerning Different Disease History") +
  theme_bw() +
  theme(strip.placement = "outside") +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

dev.off()


ggplot(Performances,
       aes(y = Accuracy,
           x = Population,
           color = Disease)) +
  
  geom_point(size = 2) +
  geom_line(aes(group = Disease)) +
  facet_grid(cols = vars(Model), switch="both") +
  labs(y = "Accuracy") +
  theme_bw() +
  theme(strip.placement = "outside") +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    axis.title.y = element_text(size = 9),
    strip.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

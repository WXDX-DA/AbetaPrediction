
library(dplyr)
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

ADNI = read.csv("../data/ADNI.csv")%>%
  mutate(PET = as.factor(ifelse(PET == "1", "Positive", "Negative")))

table(ADNI$Diagnosis)



# Model Validation in ADNI

ROC_performance = function(df, ROCFL, Prob, cp, index) {
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
  
  ci.bot = ci(ROCFL)[1]
  ci.top = ci(ROCFL)[3]
  AUC = ROCFL$auc
  
  tab = table(Prob[, 2] > TH, df$PET)
  
  (Sen = tab[2, 2] / sum(tab[, 2]))
  (Sep = tab[1, 1] / sum(tab[, 1]))
  (PPV = tab[2, 2] / sum(tab[2, ]))
  (NPV = tab[1, 1] / sum(tab[1, ]))
  (TC = sum(diag(tab)) / sum(tab))
  (CVER = cp[index, 4])
  res = matrix(c(AUC, ci.bot, ci.top, TH, Sen, Sep, PPV, NPV, TC, CVER), 10, 1)
  colnames(res) = "Result"
  rownames(res) = c(
    "AUC",
    "ci.bot",
    "ci.top",
    "Threshold",
    "Sensitivity",
    "Sepecificity",
    "PPV",
    "NPV",
    "Accuracy",
    "CV.ErrorRate"
  )
  return(round(t(res), 3))
}

Plot_AUCFL = function(ROCs) {
  Color = c("black", hex_codes1[1], hex_codes1[4])
  finalplot = ggroc(ROCs, legacy.axes = T, size = 1.2) +
    theme_bw() +
    xlab(label = "1 - Specificity") +
    ylab(label = "Sensitivity") +
    scale_color_manual(
      labels = c("Full Model", "Best Model", "Refined Model"),
      values = Color
    ) +
    geom_segment(
      aes(
        x = 0,
        xend = 1,
        y = 0,
        yend = 1
      ),
      color = "darkgrey",
      linetype = "dashed",
      size = 1.2
    ) +
    
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(
      legend.position = c(0.87, 0.25),
      legend.background = element_rect(fill = "white", color = "black")
    ) +
    theme(legend.title = element_blank())
  
  return(finalplot)
}



## Full Data

ADNI.DTdata = ADNI %>%
  select("PET","Gender","Age","Edu_yrs","APOE","TAU","AB42","AB40",
         "NFL","pTau181","AB42.40", "MMSE")

ROCFL.ADNI.full = Get_ROC_fromTree(
  ADNI.DTdata,
  "PET ~.",
  MAXDEPTH = length(colnames(ADNI.DTdata)),
  MINSPLIT = 8,
  MINBUCKET = 4
)


ADNI.Full_model_performance = ROC_performance(ADNI.DTdata,
                ROCFL.ADNI.full$ROC,
                ROCFL.ADNI.full$prob,
                ROCFL.ADNI.full$cp,
                ROCFL.ADNI.full$Index)  

ROCFL.ADNI.best = Get_ROC_fromTree(
  ADNI.DTdata,
  "PET ~.",
  MAXDEPTH = length(colnames(ADNI.DTdata)),
  MINSPLIT = 8,
  MINBUCKET = 4
)

ADNI.Best_model_performance = ROC_performance(ADNI.DTdata,
                ROCFL.ADNI.best$ROC,
                ROCFL.ADNI.best$prob,
                ROCFL.ADNI.best$cp,
                ROCFL.ADNI.best$Index)  

ADNI.RefinedData = ADNI %>%
  select("PET","MMSE","pTau181","AB42.40","Edu_yrs", "Age")
ROCFL.ADNI_refined = Get_ROC_fromTree(
  ADNI.RefinedData,
  "PET ~.",
  MAXDEPTH = length(colnames(ADNI.RefinedData)),
  MINSPLIT = 8,
  MINBUCKET = 4
)

ADNI.Refined_model_performance = ROC_performance(ADNI.RefinedData,
                ROCFL.ADNI_refined$ROC,
                ROCFL.ADNI_refined$prob,
                ROCFL.ADNI_refined$cp,
                ROCFL.ADNI_refined$Index)  

ROCFL.ADNI_Limited = Get_ROC_fromTree(
  ADNI.RefinedData,
  "PET ~.",
  MAXDEPTH = 3,
  MINSPLIT = 8,
  MINBUCKET = 4
)

ADNI.Limited_model_performance = ROC_performance(ADNI.RefinedData,
                ROCFL.ADNI_Limited$ROC,
                ROCFL.ADNI_Limited$prob,
                ROCFL.ADNI_Limited$cp,
                ROCFL.ADNI_Limited$Index)  

ADNI_performance = rbind(
  ADNI.Full_model_performance,
  ADNI.Best_model_performance,
  ADNI.Refined_model_performance,
  ADNI.Limited_model_performance
)

rownames(ADNI_performance) = c("Full Model", "Best Model", "Refined Model", "Limited Model")

ADNI_performance %>%
  kbl() %>%
  kable_styling()


### AUC plot

ROC = list(ROCFL.ADNI.full$ROC,ROCFL.ADNI_refined$ROC,ROCFL.ADNI_Limited$ROC)
png(
  "../output/Figure4_1.png",
  units = "mm",
  width = 180,
  height = 140,
  res = 300,
  pointsize = 7
)
print(Plot_AUCFL(ROC))
dev.off()

print(Plot_AUCFL(ROC))


### Tree Plot

rpart.plot(ROCFL.ADNI_Limited$Model, roundint = FALSE)


# MCI

## Performance

ADNI.MCI.DTdata = ADNI %>%
  filter(Diagnosis == "MCI") %>%
  select("PET","Gender","Age","Edu_yrs","APOE","TAU","AB42","AB40",
         "NFL","pTau181","AB42.40", "MMSE")

ROCFL.ADNI.MCI.full = Get_ROC_fromTree(
  ADNI.MCI.DTdata,
  "PET ~.",
  MAXDEPTH = length(colnames(ADNI.MCI.DTdata)),
  MINSPLIT = 8,
  MINBUCKET = 4
)

ADNI.MCI.Full_model_performance = ROC_performance(ADNI.MCI.DTdata,
                ROCFL.ADNI.MCI.full$ROC,
                ROCFL.ADNI.MCI.full$prob,
                ROCFL.ADNI.MCI.full$cp,
                ROCFL.ADNI.MCI.full$Index)  

ROCFL.ADNI.MCI.best = Get_ROC_fromTree(
  ADNI.MCI.DTdata,
  "PET ~.",
  MAXDEPTH = length(colnames(ADNI.MCI.DTdata)),
  MINSPLIT = 8,
  MINBUCKET = 4
)

ADNI.MCI.Best_model_performance = ROC_performance(ADNI.MCI.DTdata,
                ROCFL.ADNI.MCI.best$ROC,
                ROCFL.ADNI.MCI.best$prob,
                ROCFL.ADNI.MCI.best$cp,
                ROCFL.ADNI.MCI.best$Index)  

ADNI.MCI.RefinedData = ADNI %>%
  filter(Diagnosis == "MCI") %>%
  select("PET","APOE","pTau181","AB42.40","AB40")

ROCFL.ADNI.MCI_refined = Get_ROC_fromTree(
  ADNI.MCI.RefinedData,
  "PET ~.",
  MAXDEPTH = length(colnames(ADNI.MCI.RefinedData)),
  MINSPLIT = 8,
  MINBUCKET = 4
)

ADNI.MCI.Refined_model_performance = ROC_performance(ADNI.MCI.RefinedData,
                ROCFL.ADNI.MCI_refined$ROC,
                ROCFL.ADNI.MCI_refined$prob,
                ROCFL.ADNI.MCI_refined$cp,
                ROCFL.ADNI.MCI_refined$Index)  


ROCFL.ADNI.MCI_Limited = Get_ROC_fromTree(
  ADNI.MCI.RefinedData,
  "PET ~.",
  MAXDEPTH = 3,
  MINSPLIT = 8,
  MINBUCKET = 4
)

ADNI.MCI.Limited_model_performance = ROC_performance(ADNI.MCI.RefinedData,
                ROCFL.ADNI.MCI_Limited$ROC,
                ROCFL.ADNI.MCI_Limited$prob,
                ROCFL.ADNI.MCI_Limited$cp,
                ROCFL.ADNI.MCI_Limited$Index)  

ADNI.MCI_performance = rbind(
  ADNI.MCI.Full_model_performance,
  ADNI.MCI.Best_model_performance,
  ADNI.MCI.Refined_model_performance,
  ADNI.MCI.Limited_model_performance
)


rownames(ADNI.MCI_performance) = c("Full Model", "Best Model", "Refined Model", "Limited Model")

ADNI.MCI_performance %>%
  kbl() %>%
  kable_styling()


## AUC plot

png(
  "../output/Figure4_2.png",
  units = "mm",
  width = 180,
  height = 140,
  res = 300,
  pointsize = 7
)
ROC = list(ROCFL.ADNI.MCI.full$ROC,
           ROCFL.ADNI.MCI_refined$ROC,
           ROCFL.ADNI.MCI_Limited$ROC)

print(Plot_AUCFL(ROC))
dev.off()

print(Plot_AUCFL(ROC))


## Tree plot

refined.tree = rpart(
  formula = "PET ~.",
  ADNI.MCI.RefinedData,
  control =  rpart.control(
    cp = 0,
    maxdepth = 3,
    minsplit = 8,
    minbucket = 4
  )
)
BT = prune(refined.tree, cp = refined.tree$cptable[4, "CP"])
rpart.plot(BT)





# NC


ADNI.NC.DTdata = ADNI %>%
  filter(Diagnosis == "NC") %>%
  select("PET","Gender","Age","Edu_yrs","APOE","TAU","AB42","AB40",
         "NFL","pTau181","AB42.40", "MMSE")

ROCFL.ADNI.NC.full = Get_ROC_fromTree(
  ADNI.NC.DTdata,
  "PET ~.",
  MAXDEPTH = length(colnames(ADNI.NC.DTdata)),
  MINSPLIT = 8,
  MINBUCKET = 4
)


ADNI.NC.Full_model_performance = ROC_performance(ADNI.NC.DTdata,
                ROCFL.ADNI.NC.full$ROC,
                ROCFL.ADNI.NC.full$prob,
                ROCFL.ADNI.NC.full$cp,
                ROCFL.ADNI.NC.full$Index)  

ROCFL.ADNI.NC.best = Get_ROC_fromTree(
  ADNI.NC.DTdata,
  "PET ~.",
  MAXDEPTH = length(colnames(ADNI.NC.DTdata)),
  MINSPLIT = 8,
  MINBUCKET = 4
)

ADNI.NC.Best_model_performance = ROC_performance(ADNI.NC.DTdata,
                ROCFL.ADNI.NC.best$ROC,
                ROCFL.ADNI.NC.best$prob,
                ROCFL.ADNI.NC.best$cp,
                ROCFL.ADNI.NC.best$Index)  


ADNI.NC.RefinedData = ADNI %>%
  filter(Diagnosis == "NC") %>%
  select("PET","APOE","pTau181","NFL","Age")

ROCFL.ADNI.NC_refined = Get_ROC_fromTree(
  ADNI.NC.RefinedData,
  "PET ~.",
  MAXDEPTH = length(colnames(ADNI.NC.RefinedData)),
  MINSPLIT = 8,
  MINBUCKET = 4
)


ADNI.NC.Refined_model_performance = ROC_performance(ADNI.NC.RefinedData,
                ROCFL.ADNI.NC_refined$ROC,
                ROCFL.ADNI.NC_refined$prob,
                ROCFL.ADNI.NC_refined$cp,
                ROCFL.ADNI.NC_refined$Index)  



ROCFL.ADNI.NC_Limited = Get_ROC_fromTree(
  ADNI.NC.RefinedData,
  "PET ~.",
  MAXDEPTH = 3,
  MINSPLIT = 8,
  MINBUCKET = 4
)

ADNI.NC.Limited_model_performance = ROC_performance(ADNI.NC.RefinedData,
                ROCFL.ADNI.NC_Limited$ROC,
                ROCFL.ADNI.NC_Limited$prob,
                ROCFL.ADNI.NC_Limited$cp,
                ROCFL.ADNI.NC_Limited$Index)  

ADNI.NC_performance = rbind(
  ADNI.NC.Full_model_performance,
  ADNI.NC.Best_model_performance,
  ADNI.NC.Refined_model_performance,
  ADNI.NC.Limited_model_performance
)

rownames(ADNI.NC_performance) = c("Full Model", "Best Model", "Refined Model", "Limited Model")

ADNI.NC_performance %>%
  kbl() %>%
  kable_styling()


## Tree Plot

refined.tree = rpart(
  formula = "PET ~.",
  ADNI.NC.RefinedData,
  control =  rpart.control(
    cp = 0,
    maxdepth = 3,
    minsplit = 8,
    minbucket = 4
  )
)
BT = prune(refined.tree, cp = refined.tree$cptable[2, "CP"])
rpart.plot(BT)



# Plot AUC Difference

fullData_models = readRDS("../data/fullData_models.RDS")
MCIData_models = readRDS("../data/MCIData_models.RDS")
NCData_models = readRDS("../data/NCData_models.RDS")
NCMCIADData_models = readRDS("../data/NCMCIAD.Data_models.RDS")
NCMCIADData_models$full$cls
get_performance = function(Data_models, Population, Data_source) {
  Auc_range_full = c(
    auc(Data_models$full$ROC),
    ci(Data_models$full$ROC)[1],
    ci(Data_models$full$ROC)[3],
    "Full", Population, Data_source
  )
  Auc_range_Best = c(
    auc(Data_models$Best$ROC),
    ci(Data_models$Best$ROC)[1],
    ci(Data_models$Best$ROC)[3],
    "Best", Population, Data_source
  )
  Auc_range_Refined = c(
    auc(Data_models$Refined$ROC),
    ci(Data_models$Refined$ROC)[1],
    ci(Data_models$Refined$ROC)[3],
    "Refined", Population, Data_source
  )
  
  data_performance = data.frame(rbind(Auc_range_full, 
                                      Auc_range_Best,
                                      Auc_range_Refined))
    
    
  colnames(data_performance) = c("AUC", "ci.bot", "ci.top", "Model", "Population", "Data_source")
  data_performance = data_performance %>%
    mutate(AUC = as.numeric(AUC),
           ci.bot = as.numeric(ci.bot),
           ci.top = as.numeric(ci.top))
  return(data_performance)
}
GP.performance.cn = get_performance(fullData_models,
                                    "All Participants",
                                    "Training Cohort (with SCD & non-AD)")
NCMCIAD.performance.cn = get_performance(NCMCIADData_models, "All Participants", "Training Cohort (without SCD & non-AD)")
MCI.performance.cn = get_performance(MCIData_models, "MCI", "Training Cohort (without SCD & non-AD)")
NC.performance.cn = get_performance(NCData_models, "CN", "Training Cohort (without SCD & non-AD)")


GP.performance.ADNI = data.frame(ADNI_performance) %>%
  select(1:3) %>%
  slice(-2) %>%
  mutate(Model = c("Full", "Best", "Refined"),
         Population = "All Participants",
         Data_source = "ADNI")

MCI.performance.ADNI = data.frame(ADNI.MCI_performance) %>%
  select(1:3) %>%
  slice(-2) %>%
  mutate(Model = c("Full", "Best", "Refined"),
         Population = "MCI",
         Data_source = "ADNI")

NC.performance.ADNI = data.frame(ADNI.NC_performance) %>%
  select(1:3) %>%
  slice(-2) %>%
  mutate(
    Model = c("Full", "Best", "Refined"),
    Population = "CN",
    Data_source = "ADNI"
  )

All_performances = bind_rows(
  GP.performance.cn,
  NCMCIAD.performance.cn,
  NC.performance.cn,
  MCI.performance.cn,
  GP.performance.ADNI,
  NC.performance.ADNI,
  MCI.performance.ADNI
  
) %>%
  mutate(
    Model = factor(Model, levels = c("Full", "Best", "Refined")),
    Data_source = factor(
      Data_source,
      levels = c("Training Cohort (with SCD & non-AD)", "Training Cohort (without SCD & non-AD)", "ADNI")
    ),
    Population = factor(Population, levels = c("All Participants", "CN", "MCI"))
  )

All_performances


tiff(
  "../output/Figure 4.tiff",
  units = "mm",
  width = 185,
  height = 145,
  res = 300,
  pointsize = 10
)

ggplot(All_performances,
       aes(y = AUC,
           x = Model)) +
  
  geom_pointrange(aes(
    ymin = ci.bot,
    ymax = ci.top,
    color = Data_source
  ),
  size = 0.5,
  position = position_dodge(0.4)) +
  facet_grid(cols = vars(Population), switch="both") +
  scale_colour_manual(values = c("black",hex_codes1[1], hex_codes1[4]))+
  ylim(0.6, 1)+
  labs(y = "AUROC (95% CI)") +
  ggtitle("Model Validation in ADNI Cohort with AUROC Values")+
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

ggplot(All_performances,
       aes(y = AUC,
           x = Model)) +
  geom_pointrange(aes(
    ymin = ci.bot,
    ymax = ci.top,
    color = Data_source
  ),
  size = 0.5,
  position = position_dodge(0.4)) +
  facet_grid(cols = vars(Population), switch="both") +
  scale_colour_manual(values = c("black",hex_codes1[1], hex_codes1[4]))+
  ylim(0.6, 1)+
  
  labs(y = "AUROC (95% CI)") +
  theme_bw() +
  theme(strip.placement = "outside") +
  theme(
    axis.title.x = element_blank(),
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    axis.title.y = element_text(size = 7),
    strip.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.position = "bottom"
  )

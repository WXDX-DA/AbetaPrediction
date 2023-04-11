
library(dplyr)
library(foreign)
library(kableExtra)
library(pROC) ## ROC
library(caret) ## ROC
library(rpart)
library(rpart.plot)
library(gtsummary)
library(multicon)
library("scales") 
# library(multicon) # scale2
hex_codes1 <- hue_pal()(6)  


# Data Preparation

ad.raw = read.csv("../data/ad_full_data.csv")%>%
  mutate(PET = as.factor(PET),
         Diagnosis = factor(Diagnosis, levels = c("NC", "SCD", "MCI", "AD", "non-AD")),
         Gender = factor(Gender, levels = c("Male", "Female")),
         APOE_risk = as.factor(APOE_risk)) %>%
  mutate_if(is.numeric, scale2, center = TRUE, scale = TRUE)

dim(ad.raw)


# Main Results

## Model in Full Data

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

### Model Performance

Full.DTdata = ad.raw %>%
  select("PET","Gender","Age","Edu_yrs","BMI" ,"APOE_risk","TAU","AB42","AB40",
         "NFL","pTau181","AB42_40", "MMSE","MoCA_B","ACEIII")

ROCFL.full = Get_ROC_fromTree(
  Full.DTdata,
  "PET ~.",
  MAXDEPTH = length(colnames(Full.DTdata)),
  MINSPLIT = 8,
  MINBUCKET = 4
)

full_model_performance = ROC_performance(Full.DTdata,
                ROCFL.full$ROC,
                ROCFL.full$prob,
                ROCFL.full$cp,
                ROCFL.full$Index) 

ROCFL.Best = Get_ROC_fromTree(
  Full.DTdata,
  "PET ~.",
  MAXDEPTH = length(colnames(Full.DTdata)),
  MINSPLIT = 8,
  MINBUCKET = 4
)

Best_model_performance = ROC_performance(Full.DTdata,
                ROCFL.Best$ROC,
                ROCFL.Best$prob,
                ROCFL.Best$cp,
                ROCFL.Best$Index)  


Refined.DTdata = ad.raw %>%
  select("PET","Age","Edu_yrs",
         "pTau181","AB42_40", "MMSE")

ROCFL.Refined = Get_ROC_fromTree(
  Refined.DTdata,
  "PET ~.",
  MAXDEPTH = length(colnames(Refined.DTdata)),
  MINSPLIT = 8,
  MINBUCKET = 4
)

Refined_model_performance = ROC_performance(Refined.DTdata,
                ROCFL.Refined$ROC,
                ROCFL.Refined$prob,
                ROCFL.Refined$cp,
                ROCFL.Refined$Index)  


ROCFL.limited = Get_ROC_fromTree(
  Refined.DTdata,
  "PET ~.",
  MAXDEPTH = 3,
  MINSPLIT = 8,
  MINBUCKET = 4
)

Limited_model_performance = ROC_performance(Refined.DTdata,
                ROCFL.limited$ROC,
                ROCFL.limited$prob,
                ROCFL.limited$cp,
                ROCFL.limited$Index)  

fullData_models = list("full" = ROCFL.full,
                       "Best" = ROCFL.Refined,
                       "Refined" = ROCFL.limited)
saveRDS(fullData_models, "../output/fullData_models.RDS")

full_performance = rbind(
  full_model_performance,
  Best_model_performance,
  Refined_model_performance,
  Limited_model_performance
)

rownames(full_performance) = c("Full Model", "Best Model", "Refined Model", "Limited Model")

full_performance %>%
  kbl() %>%
  kable_styling()


### AUC plot

Plot_AUCFL = function(ROCs) {
  Color = c("black", hex_codes1[1], hex_codes1[4])
  finalplot = ggroc(ROCs, legacy.axes = T, size = 1.2) +
    theme_bw() +
    xlab(label = "1 - Specificity") +
    ylab(label = "Sensitivity") +
    ggtitle("Model Selection Process and Performance in All Participants") + 
    scale_color_manual(
      labels = c("Full Model (AUC=0.935)", "Best Model (AUC=0.830)", "Refined Model (AUC=0.708)"),
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
      legend.position = c(0.75, 0.25),
      legend.background = element_rect(fill = "white", color = "black")
    ) +
    theme(legend.title = element_blank())+
    theme(
    axis.text = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
  
  return(finalplot)
}

ROC = list(ROCFL.full$ROC,ROCFL.Refined$ROC, ROCFL.limited$ROC)



tiff(
  "../output/temp",
  units = "mm",
  width = 250,
  height = 120,
  res = 500,
  pointsize = 10
)
ggplot() +
  ggtitle("Model Selection Process and Performance in All Participants") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
dev.off()


tiff(
  "../output/Full_Data.tiff",
  units = "mm",
  width = 165,
  height = 120,
  res = 500,
  pointsize = 10
)
print(Plot_AUCFL(ROC))
dev.off()
print(Plot_AUCFL(ROC))



rpart.plot(ROCFL.limited$Model, roundint = FALSE,
  extra = 101,
  box.palette="auto",
  branch.lty = 3,
  shadow.col = "gray",
  nn = F,
  yesno = 2,
  faclen = 3,
  digits = 3
)


## Model in Matched Population with ADNI


NCMCIAD.DTdata = ad.raw %>%
  filter(Diagnosis %in% c("NC", "MCI", "AD"))  %>%
  select("PET","Gender","Age","Edu_yrs","BMI" ,"APOE_risk","TAU","AB42","AB40",
         "NFL","pTau181","AB42_40", "MMSE","MoCA_B","ACEIII") 

ROCFL.full.NCMCIAD = Get_ROC_fromTree(
  NCMCIAD.DTdata,
  "PET ~.",
  MAXDEPTH = length(colnames(NCMCIAD.DTdata)),
  MINSPLIT = 8,
  MINBUCKET = 4
)

full_model_performance.NCMCIAD = ROC_performance(NCMCIAD.DTdata,
                ROCFL.full.NCMCIAD$ROC,
                ROCFL.full.NCMCIAD$prob,
                ROCFL.full.NCMCIAD$cp,
                ROCFL.full.NCMCIAD$Index) 

ROCFL.Best.NCMCIAD = Get_ROC_fromTree(
  NCMCIAD.DTdata,
  "PET ~.",
  MAXDEPTH = length(colnames(NCMCIAD.DTdata)),
  MINSPLIT = 8,
  MINBUCKET = 4
)

Best_model_performance.NCMCIAD = ROC_performance(NCMCIAD.DTdata,
                ROCFL.Best.NCMCIAD$ROC,
                ROCFL.Best.NCMCIAD$prob,
                ROCFL.Best.NCMCIAD$cp,
                ROCFL.Best.NCMCIAD$Index)  

NCMCIAD.Refined.DTdata = ad.raw %>%
  filter(Diagnosis %in% c("NC", "MCI", "AD"))  %>%
  select("PET","Age","Edu_yrs",
         "pTau181","AB42_40", "MMSE") 

ROCFL.Refined.NCMCIAD = Get_ROC_fromTree(
  NCMCIAD.Refined.DTdata,
  "PET ~.",
  MAXDEPTH = length(colnames(NCMCIAD.Refined.DTdata)),
  MINSPLIT = 8,
  MINBUCKET = 4
)


Refined_model_performance.NCMCIAD = ROC_performance(NCMCIAD.Refined.DTdata,
                ROCFL.Refined.NCMCIAD$ROC,
                ROCFL.Refined.NCMCIAD$prob,
                ROCFL.Refined.NCMCIAD$cp,
                ROCFL.Refined.NCMCIAD$Index)  

ROCFL.limited.NCMCIAD = Get_ROC_fromTree(
  NCMCIAD.Refined.DTdata,
  "PET ~.",
  MAXDEPTH = 3,
  MINSPLIT = 8,
  MINBUCKET = 4
)

Limited_model_performance.NCMCIAD = ROC_performance(NCMCIAD.Refined.DTdata,
                ROCFL.limited.NCMCIAD$ROC,
                ROCFL.limited.NCMCIAD$prob,
                ROCFL.limited.NCMCIAD$cp,
                ROCFL.limited.NCMCIAD$Index)  

NCMCIAD.Data_models = list("full" = ROCFL.full.NCMCIAD,
                       "Best" = ROCFL.Refined.NCMCIAD,
                       "Refined" = ROCFL.limited.NCMCIAD)
saveRDS(NCMCIAD.Data_models, "../output/NCMCIAD.Data_models.RDS")

full_performance.NCMCIAD = rbind(
  full_model_performance.NCMCIAD,
  Best_model_performance.NCMCIAD,
  Refined_model_performance.NCMCIAD,
  Limited_model_performance.NCMCIAD
)

rownames(full_performance.NCMCIAD) = c("Full Model", "Best Model", "Refined Model", "Limited Model")

full_performance.NCMCIAD %>%
  kbl() %>%
  kable_styling()


## NC

NC.DTdata = ad.raw %>%
  filter(Diagnosis %in% c("NC")) %>%
  select("PET","Gender","Age","Edu_yrs","BMI" ,"APOE_risk","TAU","AB42","AB40",
         "NFL","pTau181","AB42_40", "MMSE","MoCA_B","ACEIII")

ROCFL.NC.full = Get_ROC_fromTree(
  NC.DTdata,
  "PET ~.",
  MAXDEPTH = length(colnames(NC.DTdata)),
  MINSPLIT = 8,
  MINBUCKET = 4
)

NC.Full_model_performance = ROC_performance(NC.DTdata,
                ROCFL.NC.full$ROC,
                ROCFL.NC.full$prob,
                ROCFL.NC.full$cp,
                ROCFL.NC.full$Index)  

ROCFL.NC.best = Get_ROC_fromTree(
  NC.DTdata,
  "PET ~.",
  MAXDEPTH = length(colnames(NC.DTdata)),
  MINSPLIT = 8,
  MINBUCKET = 4
)

NC.Best_model_performance = ROC_performance(NC.DTdata,
                ROCFL.NC.best$ROC,
                ROCFL.NC.best$prob,
                ROCFL.NC.best$cp,
                ROCFL.NC.best$Index)  

NC.RefinedData = ad.raw %>%
  filter(Diagnosis %in% c("NC")) %>%
  select("PET","NFL","Age","pTau181","APOE_risk")

ROCFL.NC_refined = Get_ROC_fromTree(
  NC.RefinedData,
  "PET ~.",
  MAXDEPTH = length(colnames(NC.RefinedData)),
  MINSPLIT = 8,
  MINBUCKET = 4
)

NC.Refined_model_performance = ROC_performance(NC.RefinedData,
                ROCFL.NC_refined$ROC,
                ROCFL.NC_refined$prob,
                ROCFL.NC_refined$cp,
                ROCFL.NC_refined$Index)  

ROCFL.NC_Limited = Get_ROC_fromTree(
  NC.RefinedData,
  "PET ~.",
  MAXDEPTH = 3,
  MINSPLIT = 8,
  MINBUCKET = 4
)

NC.Limited_model_performance = ROC_performance(NC.RefinedData,
                ROCFL.NC_Limited$ROC,
                ROCFL.NC_Limited$prob,
                ROCFL.NC_Limited$cp,
                ROCFL.NC_Limited$Index) 

NCData_models = list("full" = ROCFL.NC.full,
                       "Best" = ROCFL.NC_refined,
                       "Refined" = ROCFL.NC_Limited)
saveRDS(NCData_models, "../output/NCData_models.RDS")

NC_performance = rbind(
  NC.Full_model_performance,
  NC.Best_model_performance,
  NC.Refined_model_performance,
  NC.Limited_model_performance
)

rownames(NC_performance) = c("Full Model", "Best Model", "Refined Model", "Limited Model")

NC_performance %>%
  kbl() %>%
  kable_styling()


### Tree Plot
rpart.plot(ROCFL.NC_Limited$Model, roundint = FALSE,
  extra = 101,
  box.palette="auto",
  branch.lty = 3,
  shadow.col = "gray",
  nn = F,
  yesno = 2,
  faclen = 3,
  digits = 3
)

## SCD

SCD.DTdata = ad.raw %>%
  filter(Diagnosis %in% c("SCD")) %>%
  select("PET","Gender","Age","Edu_yrs","BMI" ,"APOE_risk","TAU","AB42","AB40",
         "NFL","pTau181","AB42_40", "MMSE","MoCA_B","ACEIII")

ROCFL.SCD.full = Get_ROC_fromTree(
  SCD.DTdata,
  "PET ~.",
  MAXDEPTH = length(colnames(SCD.DTdata)),
  MINSPLIT = 8,
  MINBUCKET = 4
)

SCD.Full_model_performance = ROC_performance(SCD.DTdata,
                ROCFL.SCD.full$ROC,
                ROCFL.SCD.full$prob,
                ROCFL.SCD.full$cp,
                ROCFL.SCD.full$Index)  


ROCFL.SCD.best = Get_ROC_fromTree(
  SCD.DTdata,
  "PET ~.",
  MAXDEPTH = length(colnames(SCD.DTdata)),
  MINSPLIT = 8,
  MINBUCKET = 4
)

SCD.Best_model_performance = ROC_performance(SCD.DTdata,
                ROCFL.SCD.best$ROC,
                ROCFL.SCD.best$prob,
                ROCFL.SCD.best$cp,
                ROCFL.SCD.best$Index)  

SCD.RefinedData = ad.raw %>%
  filter(Diagnosis %in% c("SCD")) %>%
  select("PET","AB40","Edu_yrs")

ROCFL.SCD_refined = Get_ROC_fromTree(
  SCD.RefinedData,
  "PET ~.",
  MAXDEPTH = length(colnames(SCD.RefinedData)),
  MINSPLIT = 8,
  MINBUCKET = 4
)

SCD.Refined_model_performance = ROC_performance(SCD.RefinedData,
                ROCFL.SCD_refined$ROC,
                ROCFL.SCD_refined$prob,
                ROCFL.SCD_refined$cp,
                ROCFL.SCD_refined$Index)  


ROCFL.SCD_Limited = Get_ROC_fromTree(
  SCD.RefinedData,
  "PET ~.",
  MAXDEPTH = 3,
  MINSPLIT = 8,
  MINBUCKET = 4
)

SCD.Limited_model_performance = ROC_performance(SCD.RefinedData,
                ROCFL.SCD_Limited$ROC,
                ROCFL.SCD_Limited$prob,
                ROCFL.SCD_Limited$cp,
                ROCFL.SCD_Limited$Index)  

SCD_performance = rbind(
  SCD.Full_model_performance,
  SCD.Best_model_performance,
  SCD.Refined_model_performance,
  SCD.Limited_model_performance
)

rownames(SCD_performance) = c("Full Model", "Best Model", "Refined Model", "Limited Model")

SCD_performance %>%
  kbl() %>%
  kable_styling()


### Tree Plot

rpart.plot(ROCFL.SCD_Limited$Model, roundint = F,
  extra = 101,
  box.palette="auto",
  branch.lty = 3,
  shadow.col = "gray",
  nn = F,
  yesno = 2,
  faclen = 3,
  digits = 3
)


## MCI

MCI.DTdata = ad.raw %>%
  filter(Diagnosis %in% c("MCI")) %>%
  select("PET","Gender","Age","Edu_yrs","BMI" ,"APOE_risk","TAU","AB42","AB40",
         "NFL","pTau181","AB42_40", "MMSE","MoCA_B","ACEIII")

ROCFL.MCI.full = Get_ROC_fromTree(
  MCI.DTdata,
  "PET ~.",
  MAXDEPTH = length(colnames(MCI.DTdata)),
  MINSPLIT = 8,
  MINBUCKET = 4
)


MCI.Full_model_performance = ROC_performance(MCI.DTdata,
                ROCFL.MCI.full$ROC,
                ROCFL.MCI.full$prob,
                ROCFL.MCI.full$cp,
                ROCFL.MCI.full$Index)  


ROCFL.MCI.best = Get_ROC_fromTree(
  MCI.DTdata,
  "PET ~.",
  MAXDEPTH = length(colnames(MCI.DTdata)),
  MINSPLIT = 8,
  MINBUCKET = 4
)

MCI.Best_model_performance = ROC_performance(MCI.DTdata,
                ROCFL.MCI.best$ROC,
                ROCFL.MCI.best$prob,
                ROCFL.MCI.best$cp,
                ROCFL.MCI.best$Index)  


MCI.RefinedData = ad.raw %>%
  filter(Diagnosis %in% c("MCI")) %>%
  select("PET","pTau181","AB40","AB42_40","APOE_risk")

ROCFL.MCI_refined = Get_ROC_fromTree(
  MCI.RefinedData,
  "PET ~.",
  MAXDEPTH = length(colnames(MCI.RefinedData)),
  MINSPLIT = 8,
  MINBUCKET = 4
)

MCI.Refined_model_performance = ROC_performance(MCI.RefinedData,
                ROCFL.MCI_refined$ROC,
                ROCFL.MCI_refined$prob,
                ROCFL.MCI_refined$cp,
                ROCFL.MCI_refined$Index)  


ROCFL.MCI_Limited = Get_ROC_fromTree(
  MCI.RefinedData,
  "PET ~.",
  MAXDEPTH = 3,
  MINSPLIT = 8,
  MINBUCKET = 4
)

MCI.Limited_model_performance = ROC_performance(MCI.RefinedData,
                ROCFL.MCI_Limited$ROC,
                ROCFL.MCI_Limited$prob,
                ROCFL.MCI_Limited$cp,
                ROCFL.MCI_Limited$Index)  
MCIData_models = list("full" = ROCFL.MCI.full,
                       "Best" = ROCFL.MCI_refined,
                       "Refined" = ROCFL.MCI_Limited)
saveRDS(MCIData_models, "../output/MCIData_models.RDS")

MCI_performance = rbind(
  MCI.Full_model_performance,
  MCI.Best_model_performance,
  MCI.Refined_model_performance,
  MCI.Limited_model_performance
)

rownames(MCI_performance) = c("Full Model", "Best Model", "Refined Model", "Limited Model")

MCI_performance %>%
  kbl() %>%
  kable_styling()


### Tree Plot

rpart.plot(ROCFL.MCI_Limited$Model, roundint = FALSE,
  extra = 101,
  box.palette="auto",
  branch.lty = 3,
  shadow.col = "gray",
  nn = F,
  yesno = 2,
  faclen = 3,
  digits = 3
)


## Model in aMCI

### Model Performance

aMCI.DTdata = ad.raw %>%
  filter(Diagnosis_raw %in% c("aMCI")) %>%
  select("PET","Gender","Age","Edu_yrs","BMI" ,"APOE_risk","TAU","AB42","AB40",
         "NFL","pTau181","AB42_40", "MMSE","MoCA_B","ACEIII")

ROCFL.aMCI.full = Get_ROC_fromTree(
  aMCI.DTdata,
  "PET ~.",
  MAXDEPTH = length(colnames(aMCI.DTdata)),
  MINSPLIT = 8,
  MINBUCKET = 4
)

aMCI.Full_model_performance = ROC_performance(aMCI.DTdata,
                ROCFL.aMCI.full$ROC,
                ROCFL.aMCI.full$prob,
                ROCFL.aMCI.full$cp,
                ROCFL.aMCI.full$Index)  

ROCFL.aMCI.best = Get_ROC_fromTree(aMCI.DTdata, "PET ~.", MAXDEPTH = length(colnames(aMCI.DTdata)), MINSPLIT = 8, MINBUCKET = 4)

aMCI.Best_model_performance = ROC_performance(aMCI.DTdata,
                ROCFL.aMCI.best$ROC,
                ROCFL.aMCI.best$prob,
                ROCFL.aMCI.best$cp,
                ROCFL.aMCI.best$Index)  

aMCI.RefinedData = ad.raw %>%
  filter(Diagnosis_raw %in% c("aMCI")) %>%
  select("PET","AB40","AB42","APOE_risk","pTau181")

ROCFL.aMCI_refined = Get_ROC_fromTree(
  aMCI.RefinedData,
  "PET ~.",
  MAXDEPTH = length(colnames(aMCI.RefinedData)),
  MINSPLIT = 8,
  MINBUCKET = 4
)

aMCI.Refined_model_performance = ROC_performance(aMCI.RefinedData,
                ROCFL.aMCI_refined$ROC,
                ROCFL.aMCI_refined$prob,
                ROCFL.aMCI_refined$cp,
                ROCFL.aMCI_refined$Index)  

ROCFL.aMCI_Limited = Get_ROC_fromTree(
  aMCI.RefinedData,
  "PET ~.",
  MAXDEPTH = 3,
  MINSPLIT = 8,
  MINBUCKET = 4
)

aMCI.Limited_model_performance = ROC_performance(aMCI.RefinedData,
                ROCFL.aMCI_Limited$ROC,
                ROCFL.aMCI_Limited$prob,
                ROCFL.aMCI_Limited$cp,
                ROCFL.aMCI_Limited$Index)  

aMCIData_models = list("full" = ROCFL.aMCI.full,
                       "Best" = ROCFL.aMCI_refined,
                       "Refined" = ROCFL.aMCI_Limited)
saveRDS(aMCIData_models, "../output/aMCIData_models.RDS")


aMCI_performance = rbind(
  aMCI.Full_model_performance,
  aMCI.Best_model_performance,
  aMCI.Refined_model_performance,
  aMCI.Limited_model_performance
)

rownames(aMCI_performance) = c("Full Model", "Best Model", "Refined Model", "Limited Model")

aMCI_performance %>%
  kbl() %>%
  kable_styling()

### AUC plot

ROC_aMCI = list(
  ROCFL.aMCI.full$ROC,
  ROCFL.aMCI_refined$ROC,
  ROCFL.aMCI_Limited$ROC
)

png(
  "../output/aMCI_Data.png",
  units = "mm",
  width = 180,
  height = 180,
  res = 300,
  pointsize = 7
)
print(Plot_AUCFL(ROC_aMCI))
dev.off()
print(Plot_AUCFL(ROC_aMCI))


### Tree Plot

png(
  "../output/aMCI_Tree.png",
  units = "mm",
  width = 150,
  height = 150,
  res = 500,
  pointsize = 9
)

ggplot() +
  ggtitle("Decision Tree Diagram of the Refined Model in the aMCI Group")+
  theme(plot.title =  element_text(face = "bold", hjust = 0.5, size = 12))
dev.off()

rpart.plot(
  ROCFL.aMCI_Limited$Model,
  roundint = FALSE,
  extra = 101,
  box.palette="auto",
  branch.lty = 3,
  shadow.col = "gray",
  nn = F,
  yesno = 2,
  faclen = 3,
  digits = 3
)



## Dementia

Dementia.DTdata = ad.raw %>%
  filter(Diagnosis %in% c("non-AD", "AD")) %>%
  select("PET","Gender","Age","Edu_yrs","BMI" ,"APOE_risk","TAU","AB42","AB40",
         "NFL","pTau181","AB42_40", "MMSE","MoCA_B","ACEIII")

ROCFL.Dementia.full = Get_ROC_fromTree(
  Dementia.DTdata,
  "PET ~.",
  MAXDEPTH = length(colnames(Dementia.DTdata)),
  MINSPLIT = 8,
  MINBUCKET = 4
)

Dementia.Full_model_performance = ROC_performance(Dementia.DTdata,
                ROCFL.Dementia.full$ROC,
                ROCFL.Dementia.full$prob,
                ROCFL.Dementia.full$cp,
                ROCFL.Dementia.full$Index)  


ROCFL.Dementia.best = Get_ROC_fromTree(
  Dementia.DTdata,
  "PET ~.",
  MAXDEPTH = length(colnames(Dementia.DTdata)),
  MINSPLIT = 8,
  MINBUCKET = 4
)


Dementia.Best_model_performance = ROC_performance(Dementia.DTdata,
                ROCFL.Dementia.best$ROC,
                ROCFL.Dementia.best$prob,
                ROCFL.Dementia.best$cp,
                ROCFL.Dementia.best$Index)  


Dementia.RefinedData = ad.raw %>%
  filter(Diagnosis %in% c("non-AD", "AD")) %>%
  select("PET","pTau181","NFL","APOE_risk","Age")

ROCFL.Dementia_refined = Get_ROC_fromTree(
  Dementia.RefinedData,
  "PET ~.",
  MAXDEPTH = length(colnames(Dementia.RefinedData)),
  MINSPLIT = 8,
  MINBUCKET = 4
)

Dementia.Refined_model_performance = ROC_performance(Dementia.RefinedData,
                ROCFL.Dementia_refined$ROC,
                ROCFL.Dementia_refined$prob,
                ROCFL.Dementia_refined$cp,
                ROCFL.Dementia_refined$Index)  


ROCFL.Dementia_Limited = Get_ROC_fromTree(
  Dementia.RefinedData,
  "PET ~.",
  MAXDEPTH = 3,
  MINSPLIT = 8,
  MINBUCKET = 4
)

Dementia.Limited_model_performance = ROC_performance(Dementia.RefinedData,
                ROCFL.Dementia_Limited$ROC,
                ROCFL.Dementia_Limited$prob,
                ROCFL.Dementia_Limited$cp,
                ROCFL.Dementia_Limited$Index)  

Dementia_performance = rbind(
  Dementia.Full_model_performance,
  Dementia.Best_model_performance,
  Dementia.Refined_model_performance,
  Dementia.Limited_model_performance
)

rownames(Dementia_performance) = c("Full Model", "Best Model", "Refined Model", "Limited Model")

Dementia_performance %>%
  kbl() %>%
  kable_styling()



rpart.plot(ROCFL.Dementia_Limited$Model, roundint = FALSE,
  extra = 101,
  box.palette="auto",
  branch.lty = 3,
  shadow.col = "gray",
  nn = F,
  yesno = 2,
  faclen = 3,
  digits = 3
)



# Best Model AUC in different groups


ALL_BEST_performance = rbind(
  Refined_model_performance,
  NC.Refined_model_performance,
  SCD.Refined_model_performance,
  MCI.Refined_model_performance,
  Dementia.Refined_model_performance
)

rownames(ALL_BEST_performance) = c("General Population", "NC Group", "SCD Group", 
                 "MCI Group", "Dementia Group")

ALL_BEST_performance %>%
  kbl() %>%
  kable_styling()



model_levels = factor(
  c("All Participants", "CN", "SCD",
    "MCI", "Dementia"),
  levels = c("All Participants", "CN", "SCD",
             "MCI", "Dementia")
)

ALL_Full_performance = data.frame(rbind(
  full_model_performance,
  NC.Full_model_performance,
  SCD.Full_model_performance,
  MCI.Full_model_performance,
  Dementia.Full_model_performance
), Model = "Full Model", Population = model_levels)

rownames(ALL_Full_performance) = NULL

ALL_Best_performance = data.frame(rbind(
  Refined_model_performance,
  NC.Refined_model_performance,
  SCD.Refined_model_performance,
  MCI.Refined_model_performance,
  Dementia.Refined_model_performance
), Model = "Best Model", Population = model_levels)

rownames(ALL_Best_performance) = NULL

ALL_Refined_performance = data.frame(rbind(
  Limited_model_performance,
  NC.Limited_model_performance,
  SCD.Limited_model_performance,
  MCI.Limited_model_performance,
  Dementia.Limited_model_performance
), Model = "Refined Model", Population = model_levels)

rownames(ALL_Refined_performance) = NULL

All_performances1 = bind_rows(
  ALL_Full_performance,
  ALL_Best_performance
)

All_performances = bind_rows(All_performances1,
                             ALL_Refined_performance) %>%
  mutate(Model = factor(Model, 
                        levels = c("Full Model", "Best Model", "Refined Model")))


tiff(
  "../output/Figure 2.tiff",
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
    color = Model
  ),
  size = 0.5,
  position = position_dodge(0.4)) +
  facet_grid(cols = vars(Population), switch="both") +
  scale_colour_manual(values = c("black", hex_codes1[1], hex_codes1[4]))+
  # scale_y_discrete(expand = c(0,2))+
  ylim(0.6, 1)+
  labs(y = "AUROC (95% CI)") +
  ggtitle("AUROC Values and 95% CI of the Established Models")+
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    axis.title.y = element_text(size = 10),
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
    color = Model
  ),
  size = 0.5) +
  facet_grid(cols = vars(Population), switch="both") +
  scale_colour_manual(values = c("black", hex_codes1[1], hex_codes1[4]))+
  ylim(0.5, 1)+
  labs(y = "AUROC (95% CI)") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    axis.title.y = element_text(size = 7),
    strip.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.position = "bottom"
  )



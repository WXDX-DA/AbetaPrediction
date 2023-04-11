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
         Gender = factor(Gender, levels = c("Male", "Female")),
         APOE_risk = as.factor(APOE_risk)
         )
dim(ad.raw)


str(ad.raw)


# Demographic Table

create.summary.table = function(df, separate.var, separate.name) {
  set_gtsummary_theme(list(
    `tbl_summary-fn:percent_fun` = function(x) sprintf("%.1f", x *100)
  ))
  separate.name = paste0("**", separate.name, " Status**")
  table1 = df %>%
    tbl_summary(
      by = separate.var,
      type = list(all_numeric() ~ "continuous2",
                  c(Gender, APOE_risk) ~ "dichotomous"),
      digits = list(all_continuous() ~ 0, BMI ~ 2),
      missing = "no",
      percent = "column"
    ) %>%
    modify_header(label ~ "**Variable**") %>%
    add_p() %>%
    add_overall() %>%
    as_flex_table()
  
  return(table1)
}

demo.data = ad.raw %>%
  select(c("Diagnosis", "Age", "Edu_yrs", Gender, "BMI", APOE_risk,
         MMSE, MoCA_B, ACEIII, AVLT_N5, AVLT_N7, BNT, AFT, Trails1_1, 
         Trails2_1,ADL, FAQ))%>%
  mutate(APOE_risk = as.factor(ifelse((APOE_risk %in% c("1", "2")),0, 1)),
         Gender = as.factor(ifelse(Gender == "Female",0, 1)),
         Diagnosis = factor(Diagnosis, levels = c("NC", "SCD", "MCI", "AD", "non-AD")))

create.summary.table(demo.data, "Diagnosis", "Diagnosis")



## Biomarkers

bio.demodata = ad.raw %>%
  select(c(
    "Diagnosis",
    "AB40",
    "AB42",
    "TAU",
    "pTau181" ,
    "NFL",
    "AB42_40",
    "PET"
  )) %>% 
  mutate(Diagnosis = ifelse(
    ad.raw$Diagnosis %in% c("AD", "non-AD"),
    "Dementia",
    ad.raw$Diagnosis
  )) %>%
  mutate(Diagnosis = factor(Diagnosis, levels = c("NC", "SCD", "MCI", "Dementia"))) 



create.summary.table_bio = function(status, separate.var, separate.name) {
  df = bio.demodata %>%
    filter(Diagnosis %in% status) %>%
    select(-Diagnosis)
  separate.name = paste0("**", separate.name, " Status**")
  table1 = df %>%
    tbl_summary(
      by = separate.var,
      type = list(all_numeric() ~ "continuous2"),
      statistic = list(all_continuous() ~ "{mean} ({sd})"),
      digits = list(all_continuous() ~ c(1, 1),`AB42_40` ~ c(4, 3)),
      missing = "no",
      percent = "column"
    ) %>%
    modify_header(label ~ "**Variable**") %>%
    add_p() %>%
    as_flex_table()
  
  return(table1)
}


### Overall


create.summary.table_bio(unique(bio.demodata$Diagnosis), "PET", "PET")


### NC

create.summary.table_bio("NC", "PET", "PET")


### SCD

create.summary.table_bio("SCD", "PET", "PET")


### MCI

create.summary.table_bio("MCI", "PET", "PET")


### Dementia

create.summary.table_bio("Dementia", "PET", "PET")



## APOE

create.summary.table_APOE = function(df, separate.var, separate.name) {
  set_gtsummary_theme(list(
    `tbl_summary-fn:percent_fun` = function(x) sprintf("%.1f", x *100)
  ))
  separate.name = paste0("**", separate.name, " Status**")
  table1 = df %>%
    tbl_summary(
      by = separate.var,
      type = list(all_numeric() ~ "continuous2"),
      digits = list(all_continuous() ~ 2),
      missing = "no",
      percent = "column"
    ) %>%
    modify_header(label ~ "**Variable**") %>%
    add_overall() %>%
    as_flex_table()
  
  return(table1)
}
APOE.data = ad.raw %>%
  select(c(
    "Diagnosis",
    "APOE_risk"
  )) %>%
  mutate(Diagnosis = factor(Diagnosis, levels = c("NC", "SCD", "MCI", "AD", "non-AD")),
         APOE_risk = as.factor(APOE_risk))

create.summary.table_APOE(APOE.data, "Diagnosis", "Diagnosis")





# Prediction for Dementia Status

subgroup.patients = function(df, out.group, target.group) {
  nonAD.cleaned = df %>%
    filter(!Diagnosis %in% out.group)
  
  result.df = nonAD.cleaned
  
  result.df$Diagnosis = as.factor(ifelse(result.df$Diagnosis %in% target.group, 1, 0))
  return(result.df)
  
}


ROC_performance2 = function(df, ROCFL, Prob, cp, index) {
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
  
  tab = table(Prob[, 2] > TH, df$Diagnosis)
  
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

bio.data = ad.raw %>%
  select(c(
    "PET",
    "Diagnosis",
    "AB40",
    "AB42",
    "TAU",
    "pTau181" ,
    "NFL",
    "AB42_40",
    "APOE_risk"
  ))


ADNC.general = subgroup.patients(bio.data, c("non-AD", "MCI", "SCD"), "AD")
ADSCD.general = subgroup.patients(bio.data, c("non-AD", "MCI", "NC"), "AD")
ADMCI.general = subgroup.patients(bio.data, c("non-AD", "NC", "SCD"), "AD")
SCDMCI.general = subgroup.patients(bio.data, c("non-AD", "NC", "AD"), "MCI")
NCSCD.general = subgroup.patients(bio.data, c("non-AD", "MCI", "AD"), "SCD")
NCMCI.general = subgroup.patients(bio.data, c("non-AD", "SCD", "AD"), "MCI")


general.DT = bio.data %>%
  select(-PET) %>%
  filter(Diagnosis != "non-AD") %>%
  mutate(Diagnosis = factor(Diagnosis, levels = c("NC", "SCD", "MCI", "AD")))

General.ROC1 = Get_ROC_fromTree2(general.DT, "Diagnosis ~.", MINSPLIT = 8, MINBUCKET = 4)


table(General.ROC1$cls, general.DT$Diagnosis)

table( general.DT$Diagnosis)
dim(general.DT)

General.ROC1$ROC




## Performance Measure

ADNC.general.DT = ADNC.general %>%
  select(-PET)
ADNC.ROC = Get_ROC_fromTree2(ADNC.general.DT, "Diagnosis ~.", MAXDEPTH = 5, MINSPLIT = 8, MINBUCKET = 4)

ADNC.Best_model_performance = ROC_performance2(ADNC.general.DT,
                ADNC.ROC$ROC,
                ADNC.ROC$prob,
                ADNC.ROC$cp,
                ADNC.ROC$Index)  

ADSCD.general.DT = ADSCD.general %>%
  select(-PET)
ADSCD.ROC = Get_ROC_fromTree2(ADSCD.general.DT, "Diagnosis ~.", MAXDEPTH = 5, MINSPLIT = 8, MINBUCKET = 4)

ADSCD.Best_model_performance = ROC_performance2(ADSCD.general.DT,
                ADSCD.ROC$ROC,
                ADSCD.ROC$prob,
                ADSCD.ROC$cp,
                ADSCD.ROC$Index)  

ADMCI.general.DT = ADMCI.general %>%
  select(-PET)
ADMCI.ROC = Get_ROC_fromTree2(ADMCI.general.DT, "Diagnosis ~.", MAXDEPTH = 5, MINSPLIT = 8, MINBUCKET = 4)

ADMCI.Best_model_performance = ROC_performance2(ADMCI.general.DT,
                ADMCI.ROC$ROC,
                ADMCI.ROC$prob,
                ADMCI.ROC$cp,
                ADMCI.ROC$Index)  


SCDMCI.general.DT = SCDMCI.general %>%
  select(-PET)
SCDMCI.ROC = Get_ROC_fromTree2(SCDMCI.general.DT, "Diagnosis ~.", MAXDEPTH = 5, MINSPLIT = 8, MINBUCKET = 4)

SCDMCI.Best_model_performance = ROC_performance2(SCDMCI.general.DT,
                SCDMCI.ROC$ROC,
                SCDMCI.ROC$prob,
                SCDMCI.ROC$cp,
                SCDMCI.ROC$Index)  


NCSCD.general.DT = NCSCD.general %>%
  select(-PET)
NCSCD.ROC = Get_ROC_fromTree2(NCSCD.general.DT, "Diagnosis ~.", MAXDEPTH = 5, MINSPLIT = 8, MINBUCKET = 4)

NCSCD.Best_model_performance = ROC_performance2(NCSCD.general.DT,
                NCSCD.ROC$ROC,
                NCSCD.ROC$prob,
                NCSCD.ROC$cp,
                NCSCD.ROC$Index)  


NCMCI.general.DT = NCMCI.general %>%
  select(-PET)
NCMCI.ROC = Get_ROC_fromTree2(NCMCI.general.DT, "Diagnosis ~.", MAXDEPTH = 5, MINSPLIT = 8, MINBUCKET = 4)

NCMCI.Best_model_performance = ROC_performance2(NCMCI.general.DT,
                NCMCI.ROC$ROC,
                NCMCI.ROC$prob,
                NCMCI.ROC$cp,
                NCMCI.ROC$Index)  


phase_performance = rbind(
  ADNC.Best_model_performance,
  ADSCD.Best_model_performance,
  ADMCI.Best_model_performance,
  SCDMCI.Best_model_performance,
  NCSCD.Best_model_performance,
  NCMCI.Best_model_performance
  
) 

rownames(phase_performance) = c("AD vs CN", "AD vs SCD", "AD vs MCI", "SCD vs MCI", "CN vs SCD", "CN vs MCI")

phase_performance %>%
  kbl() %>%
  kable_styling()



## AUC plot

ROCs = list(ADNC.ROC$ROC,
         ADSCD.ROC$ROC,
         ADMCI.ROC$ROC,
         SCDMCI.ROC$ROC,
         NCSCD.ROC$ROC,
         NCMCI.ROC$ROC)

Color = c(
      hex_codes1[1],
      hex_codes1[2],
      hex_codes1[3],
      hex_codes1[4],
      hex_codes1[5],
      hex_codes1[6]
    )
tiff(
  "../output/Figure 5.tiff",
  units = "mm",
  width = 185,
  height = 140,
  res = 300,
  pointsize = 10
)
#------------------------------------#
# plot the rocs
#------------------------------------#
ggroc(ROCs, legacy.axes = T, size = 1.2) +
  theme_bw() +
  xlab(label = "1 - Specificity") +
  ylab(label = "Sensitivity") +
  ggtitle("ROC Plots for Distinguishing the Patient in Different Dementia Status")+
  scale_color_manual(
    labels = c(
      "AD vs CN",
      "AD vs SCD",
      "AD vs MCI",
      "SCD vs MCI",
      "CN vs SCD",
      "CN vs MCI"
    ),
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
    legend.background = element_rect(fill = "white", color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  theme(legend.title = element_blank()) +
  theme(
    axis.text = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

dev.off()

ggroc(ROCs, legacy.axes = T, size = 1.2) +
    theme_bw() +
    xlab(label = "1 - Specificity") +
    ylab(label = "Sensitivity") +
    ggtitle("ROC Plots for Distinguishing the Patient in Different Dementia Status")+
    scale_color_manual(labels = c( "AD vs CN", "AD vs SCD", "AD vs MCI", "SCD vs MCI", "CN vs SCD", "CN vs MCI"), values = Color) +
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

    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    ) +
    theme(legend.position = c(0.87, 0.25),
        legend.background = element_rect(fill = "white", color = "black")) +
    theme(legend.title = element_blank()) 




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

ADNI = read.csv("../data/ADNI.csv")%>%
  mutate(PET = as.factor(ifelse(PET == "1", "Positive", "Negative")))

dim(ADNI)


## Data Wrangling

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

bio.ADNI = ADNI %>%
  select(c(
    "Diagnosis",
    "AB40",
    "AB42",
    "TAU",
    "pTau181" ,
    "NFL",
    "AB42.40",
    "PET"
  )) %>% 
  mutate(Diagnosis = ifelse(Diagnosis == "AD", "Dementia", Diagnosis)) %>%
  mutate(Diagnosis = factor(Diagnosis, levels = c("NC", "MCI", "Dementia")))


get_auc.singleVar = function(response, predictor) {
  roc.biomarker = roc(response,
                      predictor)
  ci.bot = ifelse(auc(roc.biomarker) < 0.5,
                  1 - ci.auc(roc.biomarker)[3],
                  ci.auc(roc.biomarker)[1])%>% round(3)
  
  ci.top = ifelse(auc(roc.biomarker) < 0.5,
                  1 - ci.auc(roc.biomarker)[1],
                  ci.auc(roc.biomarker)[3]) %>% round(3)
    
  var.auc = ifelse(auc(roc.biomarker) < 0.5,
                   1 - auc(roc.biomarker),
                   auc(roc.biomarker)) %>% round(3)
  
  result = c(var.auc, ci.bot, ci.top)
  return(result)
}

subgroup_dementia = function(df, status) {
  subgroup_data = df %>%
    filter(Diagnosis %in% status)
  return(subgroup_data)
}

summarise_auc = function(df, status) {
  
  status.ptau181.auc = get_auc.singleVar(
    subgroup_dementia(df, status)$PET,
    subgroup_dementia(df, status)$pTau181
  )
  
  status.ab4240.auc = get_auc.singleVar(
    subgroup_dementia(df, status)$PET,
    subgroup_dementia(df, status)$AB42_40
  )
  status.ab42.auc = get_auc.singleVar(
    subgroup_dementia(df, status)$PET,
    subgroup_dementia(df, status)$AB42
  )
  
  status.ab40.auc = get_auc.singleVar(
    subgroup_dementia(df, status)$PET,
    subgroup_dementia(df, status)$AB40
  )
  
  status.ttau.auc = get_auc.singleVar(
    subgroup_dementia(df, status)$PET,
    subgroup_dementia(df, status)$TAU
  )
  status.nfl.auc = get_auc.singleVar(
    subgroup_dementia(df, status)$PET,
    subgroup_dementia(df, status)$NFL
  )
  
  result = rbind(status.ab40.auc,
                 status.ab42.auc,
                 status.ttau.auc,
                 status.ptau181.auc,
                 status.nfl.auc,
                 status.ab4240.auc)
  return(result)
}


SingleAUC.table = data.frame(
  summarise_auc(bio.demodata,unique(bio.demodata$Diagnosis)),
  summarise_auc(bio.demodata,"NC"),
  summarise_auc(bio.demodata,"SCD"),
  summarise_auc(bio.demodata,"MCI"),
  summarise_auc(bio.demodata,"Dementia")
)

colnames(SingleAUC.table) = c("All.auc", "All.bot", "All.top",
                              "NC.auc", "NC.bot", "NC.top",
                              "SCD.auc", "SCD.bot", "SCD.top",
                              "MCI.auc", "MCI.bot", "MCI.top",
                              "Dementia.auc", "Dementia.bot", "Dementia.top")
rownames(SingleAUC.table) = c("Aβ40", "Aβ42", "T-tau", "P-tau181", "NfL", "Aβ42/Aβ40")


SingleAUC.table = t(SingleAUC.table) 
SingleAUC.table  = SingleAUC.table %>%
  as.data.frame() %>%
  mutate(Data_source = "Training Cohort")

Pop = factor(c("All Participants", "CN", "SCD", "MCI", "Dementia"),
             levels = c("All Participants", "CN", "SCD", "MCI", "Dementia"))
SingleAUC.table.auc = SingleAUC.table %>%
  slice(c(1,4,7, 10, 13)) %>%
  mutate(Population = Pop)
SingleAUC.table.bot = SingleAUC.table %>%
  slice(c(2,5,8, 11, 14))%>%
  mutate(Population = Pop)
SingleAUC.table.top = SingleAUC.table %>%
  slice(c(3,6,9, 12, 15))%>%
  mutate(Population = Pop)


## ADNI

summarise_auc.ADNI = function(df, status) {
  
  status.ptau181.auc = get_auc.singleVar(
    subgroup_dementia(df, status)$PET,
    subgroup_dementia(df, status)$pTau181
  )
  
  status.ab4240.auc = get_auc.singleVar(
    subgroup_dementia(df, status)$PET,
    subgroup_dementia(df, status)$AB42.40
  )
  
  status.ab42.auc = get_auc.singleVar(
    subgroup_dementia(df, status)$PET,
    subgroup_dementia(df, status)$AB42
  )
  
  status.ab40.auc = get_auc.singleVar(
    subgroup_dementia(df, status)$PET,
    subgroup_dementia(df, status)$AB40
  )
  
  status.ttau.auc = get_auc.singleVar(
    subgroup_dementia(df, status)$PET,
    subgroup_dementia(df, status)$TAU
  )
  status.nfl.auc = get_auc.singleVar(
    subgroup_dementia(df, status)$PET,
    subgroup_dementia(df, status)$NFL
  )
  
  result = rbind(status.ab40.auc,
                 status.ab42.auc,
                 status.ttau.auc,
                 status.ptau181.auc,
                 status.nfl.auc,
                 status.ab4240.auc)
  return(result)
}

SingleAUC.table_ADNI = data.frame(
  summarise_auc.ADNI(bio.ADNI,unique(bio.ADNI$Diagnosis)),
  summarise_auc.ADNI(bio.ADNI,"NC"),
  summarise_auc.ADNI(bio.ADNI,"MCI")
)

colnames(SingleAUC.table_ADNI) = c("All.auc", "All.bot", "All.top",
                              "NC.auc", "NC.bot", "NC.top",
                              "MCI.auc", "MCI.bot", "MCI.top")
rownames(SingleAUC.table_ADNI) = c("Aβ40", "Aβ42", "T-tau", "P-tau181", "NfL", "Aβ42/Aβ40")


SingleAUC.table_ADNI = t(SingleAUC.table_ADNI) 
SingleAUC.table_ADNI  = SingleAUC.table_ADNI %>%
  as.data.frame() %>%
  mutate(Data_source = "ADNI Cohort")



Pop = factor(c("All Participants", "CN", "MCI"), 
             levels = c("All Participants", "CN", "MCI"))

SingleAUC.table_ADNI.auc = SingleAUC.table_ADNI %>%
  slice(c(1,4,7)) %>%
  mutate(Population = Pop)
SingleAUC.table_ADNI.bot = SingleAUC.table_ADNI %>%
  slice(c(2,5,8))%>%
  mutate(Population = Pop)
SingleAUC.table_ADNI.top = SingleAUC.table_ADNI %>%
  slice(c(3,6,9))%>%
  mutate(Population = Pop)



## Combine Results 


aucs = bind_rows(SingleAUC.table.auc, SingleAUC.table_ADNI.auc) 
ci.bots = bind_rows(SingleAUC.table.bot, SingleAUC.table_ADNI.bot)
ci.tops = bind_rows(SingleAUC.table.top, SingleAUC.table_ADNI.top)
aucs = aucs  %>%
  pivot_longer(!c(Data_source, Population), names_to = "Biomarker", values_to = "value")
ci.bots = ci.bots  %>%
  pivot_longer(!c(Data_source, Population), names_to = "Biomarker", values_to = "value")
ci.tops = ci.tops  %>%
  pivot_longer(!c(Data_source, Population), names_to = "Biomarker", values_to = "value")

performances = data.frame(aucs, "ci.bot" = ci.bots$value, "ci.top" = ci.tops$value) %>%
  mutate(
    Data_source = factor(Data_source,
                         levels = c("Training Cohort", "ADNI Cohort")),
    Population = factor(
      Population,
      levels = c("All Participants", "CN", "SCD", "MCI", "Dementia")
    ),
    Biomarker = factor(
      Biomarker,
      levels = c("Aβ40", "Aβ42", "T-tau", "P-tau181", "NfL", "Aβ42/Aβ40")
    )
  ) 
performances %>%
  kbl() %>%
  kable_styling()


hex_codes1 <- hue_pal()(6)

tiff(
  "../output/Supplementary Figure 1.tiff",
  units = "mm",
  width = 185,
  height = 145,
  res = 500,
  pointsize = 10
)
ggplot(performances,
       aes(y = value,
           x = Biomarker)) +
  geom_pointrange(aes(
    ymin = ci.bot,
    ymax = ci.top,
    color = Population
  ),
  size = 0.5,
  position = position_dodge(0.6)) +
  facet_grid(cols = vars(Data_source), switch="both") +
  scale_colour_manual(values = c(
    "black",
    hex_codes1[2],
    hex_codes1[4],
    hex_codes1[1],
    hex_codes1[5]
  ))+
  labs(y = "AUROC (95% CI)") +
  ggtitle("AUROC Values and 95% CI of Each Individual Plasma Biomarkers")+
  theme_bw() +
  theme(strip.placement = "outside") +
  theme(
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

dev.off()

ggplot(performances,
       aes(y = value,
           x = Biomarker)) +
  geom_pointrange(aes(
    ymin = ci.bot,
    ymax = ci.top,
    color = Population
  ),
  size = 0.5,
  position = position_dodge(0.6)) +
  facet_grid(cols = vars(Data_source), switch="both") +
  scale_colour_manual(values = c(
    "black",
    hex_codes1[2],
    hex_codes1[4],
    hex_codes1[1],
    hex_codes1[5]
  ))+

  labs(y = "AUROC (95% CI)") +
  theme_bw() +
  theme(strip.placement = "outside") +
  theme(
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    axis.title.y = element_text(size = 9),
    strip.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.position = "bottom"
  )







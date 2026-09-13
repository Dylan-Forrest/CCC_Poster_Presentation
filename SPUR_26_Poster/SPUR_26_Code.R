##################################
#  Title SPUR_26_Code   
#  Date: 7/22/26
#  Author: Dylan Forrest Seay    
##################################

###### Specific Aims ######

## Aim 1: Does Adverse Childhood Experience score predict later psychopathology 

## Aim 2: Does Emotion Regulation style predict later psychopathology 

## Aim 3: Does Adverse Childhood Experience Score predict later emotion regulation strategy (Suppression and Reappraisal)


library(tidyverse)
library(tidyverse)

library(ggplot2)
library(ggplot2)

library(psych)
library(psych)


datapath <- "your/ABCD/datapath/ABCD_6.0/"

#Demographics
Demo <- read_tsv("abcd_general/ab_g_stc.tsv",
                 na=c("NA","999", "777", "888", "n/a","N/A"))

#ACEs
ACEs <- read.csv("path/to/ACEs/data/ACES_baseline_output.csv",
                 na=c("NA","999", "777", "888", "n/a","N/A"))

#Emotion Regulation
ERQ <- read_tsv("Mental_Health/mh_y_erq.tsv",
                na=c("NA","999", "777", "888", "n/a","N/A"))

#P Factor Measures
Psychopathology <- read_tsv('Mental_Health/mh_y_bpm.tsv',
                            na=c("NA","999", "777", "888", "n/a","N/A"))

#Psychotic-Like Experiences
PLEs_y <- read_tsv('Mental_Health/mh_y_pps.tsv',
                   na=c("NA","999", "777", "888", "n/a","N/A")) 



#########################
##  Demographics Data  ##
#########################

# Pull participant sex #
Participant_Sex <- Demo %>%
  select(participant_id, ab_g_stc__cohort_sex) %>%
  rename(sex = ab_g_stc__cohort_sex)


#################
##  ACEs Data  ##
#################

# Convert baseling age from months to years #
ACEs_clean <- ACEs %>%
  mutate('baseline_age' = interview_age / 12)

# Keep Relevant Columns and Rename 5.0 -> 6.0 #
ACEs_clean <- ACEs_clean %>%
  select(src_subject_id, rel_family_id, baseline_age, demo_sex_hrcode, demo_comb_income_hrcode, race_ethnicity_hrcode, Abuse_ACEs, Neglect_ACEs, ACEs_Total) %>%
  mutate(src_subject_id = str_replace(src_subject_id, "NDAR_INV", "sub-")) %>%
  rename(participant_id = src_subject_id)

# Set relevant variables as factor # 
ACEs_clean$demo_comb_income_hrcode <- as.factor(ACEs_clean$demo_comb_income_hrcode)
ACEs_clean$race_ethnicity_hrcode <- as.factor(ACEs_clean$race_ethnicity_hrcode)


###############################
##  Emotion Regulation Data  ##
###############################

# Select relevant columns and rename for readability #
ERQ_clean <- ERQ %>%
  select(participant_id, session_id, mh_y_erq_age, mh_y_erq__reapp_mean, mh_y_erq__suppr_mean) %>%
  rename('ERQ_age' = mh_y_erq_age, 'Reappraisal_Mean' = mh_y_erq__reapp_mean, 'Suppression_Mean' = mh_y_erq__suppr_mean, "ERQ_session_id" = session_id)

# Pull participants first ERQ session #
ERQ_clean <- ERQ_clean %>%
  mutate(ses_priority = match(ERQ_session_id, c("ses-03A", "ses-04A", "ses-05A", "ses-06A"))) %>%
  arrange(participant_id, ses_priority) %>%
  group_by(participant_id) %>%
  slice(1) %>%
  ungroup() %>%
  select(-ses_priority)

# Remove rows with no ERQ mean data #
ERQ_clean <- ERQ_clean %>%
  filter(Reappraisal_Mean != "n/a") %>%
  drop_na(Reappraisal_Mean) %>%
  filter(Suppression_Mean != "n/a") %>%
  drop_na(Suppression_Mean)

# Make ERQ mean scores #
ERQ_clean <- ERQ_clean %>%
  mutate(Reappraisal_Mean = as.numeric(Reappraisal_Mean)) %>%
  mutate(Suppression_Mean = as.numeric(Suppression_Mean))



#################################################
#  Generalized Psychopathology (P-Factor) Data  #   # Adapted from Shreeya Rao's code #
#################################################
 
# PLEs Dataframe #
PLEs_clean <- PLEs_y %>%
  select(participant_id,
         session_id,
         mh_y_pps__severity_score) #severity score summary

#General Psychopathology of attention, externalizing, internalizing #
Attention_clean <- Psychopathology %>%
  select(participant_id, session_id, 
         mh_y_bpm__attn_sum)  %>% #attention summary
  mutate(mh_y_bpm__attn_sum = if_else(mh_y_bpm__attn_sum =="n/a", NA, as.numeric(mh_y_bpm__attn_sum)))

Externalizing_clean <- Psychopathology %>%
  select(participant_id, session_id, 
         mh_y_bpm__ext_sum) %>% #externalizing summary
  mutate(mh_y_bpm__ext_sum = if_else(mh_y_bpm__ext_sum=="n/a", NA, as.numeric(mh_y_bpm__ext_sum))) 

Internalizing_clean <- Psychopathology %>%
  select(participant_id, session_id, 
         mh_y_bpm__int_sum) %>% #internalizing summary
  mutate(mh_y_bpm__int_sum = if_else(mh_y_bpm__int_sum=="n/a", NA, as.numeric(mh_y_bpm__int_sum)))


# Merge PLE severity scores to one number #
PLEs_wide <- PLEs_clean %>%
  pivot_wider(names_from = "session_id", values_from = "mh_y_pps__severity_score") %>%
  rowwise() %>%
  mutate(PLEs_total = sum(c(`ses-00A`, `ses-01A`, `ses-02A`, `ses-04A`, `ses-05A`, `ses-03A`, `ses-06A`), na.rm = TRUE)) %>%
  select('participant_id', 'PLEs_total')

# Merge Attention to one number #
Attention_wide <- Attention_clean %>%
  pivot_wider(names_from = "session_id", values_from = "mh_y_bpm__attn_sum") %>%
  rowwise() %>%
  mutate(Attention_total = sum(c(`ses-00M`, `ses-01A`, `ses-01M`, `ses-02A`, `ses-02M`, `ses-03A`, `ses-03M`, `ses-04A`, `ses-05A`, `ses-06A`), na.rm = TRUE)) %>% 
  select('participant_id', 'Attention_total')

# Merge Internalizing to one number #
Internalizing_wide <- Internalizing_clean %>%
  pivot_wider(names_from = "session_id", values_from = "mh_y_bpm__int_sum") %>%
  rowwise() %>%
  mutate(Internalizing_total = sum(c(`ses-00M`, `ses-01A`, `ses-01M`, `ses-02A`, `ses-02M`, `ses-03A`, `ses-03M`, `ses-04A`, `ses-05A`, `ses-06A`), na.rm = TRUE)) %>%
  select('participant_id', 'Internalizing_total')

# Merge Externalizing to one number #
Externalizing_wide <- Externalizing_clean %>%
  pivot_wider(names_from = "session_id", values_from = "mh_y_bpm__ext_sum") %>%
  rowwise() %>%
  mutate(Externalizing_total = sum(c(`ses-00M`, `ses-01A`, `ses-01M`, `ses-02A`, `ses-02M`, `ses-03A`, `ses-03M`, `ses-04A`, `ses-05A`, `ses-06A`), na.rm = TRUE)) %>%
  select('participant_id', 'Externalizing_total')

# merge PLE and General Psychopathology #
GeneralPsychopathology_clean <- PLEs_wide
GeneralPsychopathology_clean <- GeneralPsychopathology_clean %>%
  full_join (Attention_wide, by = 'participant_id') %>%
  full_join (Internalizing_wide, by = 'participant_id') %>%
  full_join (Externalizing_wide, by = 'participant_id') %>%
  drop_na()

# Create General Psychopathology score #
GeneralPsychopathology_clean <- GeneralPsychopathology_clean %>%
  mutate(P_Factor_Total = PLEs_total + Externalizing_total + Internalizing_total + Attention_total)

# Create P_Factor Z-score = (data-mean(data))/sd(data) #
GeneralPsychopathology_zscores <- GeneralPsychopathology_clean %>%
  mutate(P_Factor = (P_Factor_Total - mean(P_Factor_Total)) / sd(P_Factor_Total)) %>%
  select(participant_id, P_Factor)

describe(GeneralPsychopathology_zscores)


##################
##  Merge Data  ##
##################

# Create dataframe for all data #
SRP_Data <- ACEs_clean %>%
  left_join(ERQ_clean, by = "participant_id") %>%
  left_join(GeneralPsychopathology_zscores, by = "participant_id") %>%
  left_join(Participant_Sex, by = "participant_id")

SRP_Data <- SRP_Data %>%
  distinct(participant_id, .keep_all = TRUE)

# Remove particpants without relevant # 
SRP_Data <- SRP_Data %>%
  drop_na()

# 11,868 -> 8,636 participants for final analyses


########################
##  Models for Aims   ##
########################

## Aim 1: Does Adverse Childhood Experience score predict later psychopathology 
## Aim 2: Does Emotion Regulation style predict later psychopathology 

Psychopathology_Model <- lm(P_Factor ~ Suppression_Mean + Reappraisal_Mean + ACEs_Total + sex + demo_comb_income_hrcode + race_ethnicity_hrcode, data = SRP_Data)
summary(Psychopathology_Model)
# Supp. t = 21.22 / p = ***
# Reap. t = -6.17 / p = ***
# ACEs  t = 17.02 / p = ***


## Aim 3: Does Adverse Childhood Experience Score predict later emotion regulation strategy (Suppression and Reappraisal)

ACEs_Suppression_Model <- lm(Suppression_Mean ~ ACEs_Total + sex + demo_comb_income_hrcode + race_ethnicity_hrcode, data = SRP_Data)
summary(ACEs_Suppression_Model)
# Total t = 5.97 / p = ***

ACEs_Reappraisal_Model <- lm(Reappraisal_Mean ~ ACEs_Total + sex + demo_comb_income_hrcode + race_ethnicity_hrcode, data = SRP_Data)
summary(ACEs_Reappraisal_Model)
# Total


## Random Descriptives ##
mean(SRP_Data$ERQ_age) # 13.02
mean(SRP_Data$baseline_age) # 9.91
mean(SRP_Data$ERQ_age) - mean(SRP_Data$baseline_age) # 3.11


######################
## Graphing Time :) ##
######################

## Aim 1: ACEs and Psychopathology plot ##
ggplot(SRP_Data, aes (x = ACEs_Total, y = P_Factor)) + 
  labs(x = "Total Number of ACEs",
       y = "General Psychopathology (P Factor)") +
  coord_cartesian(ylim = c(-1,2)) +
  stat_smooth(method = "lm", col  = "black", size = 3)


## Aim 2: Emotion Regulation and Psychopathology plot ## 
ggplot(SRP_Data, aes(x = P_Factor)) +
  geom_smooth(aes(y = Reappraisal_Mean), method = "lm", color = "green", size = 3) +
  geom_smooth(aes(y = Reappraisal_Mean), method = "lm", color = "green", size = 3, linetype = "dotted") +
  geom_smooth(aes(y = Suppression_Mean), method = "lm", color = "orange", size = 3) +
  coord_cartesian(ylim = c(2.5,4.5)) +
  #coord_cartesian(ylim = c(0,5)) +
  labs(x = "General Psychopathology (P-Factor)",
       y = "Emotion Regulation Strategy Usage")


## Aim 3: Emotion Regulation and ACEs plot ##
ggplot(SRP_Data, aes(x = ACEs_Total)) +
  geom_smooth(aes(y = Reappraisal_Mean), method = "lm", color = "green", size = 3) +
  geom_smooth(aes(y = Reappraisal_Mean), method = "lm", color = "green", size = 3, linetype = "dotted") +
  geom_smooth(aes(y = Suppression_Mean), method = "lm", color = "orange", size = 3) +
  coord_cartesian(xlim = c(0,7), ylim = c(2.5,4.5)) +
  #coord_cartesian(xlim = c(0,7), ylim = c(0,5)) +
  labs(x = "Total Number of ACEs",
       y = "Emotion Regulation Strategy Usage")


## One last Sanity Check before you go ##
describe(SRP_Data)




# "     _|_     "
# "    /o o\    "   # goodbye! 
# "    \___/    "
# "      |/     "
# "     /|      "
# "     / \     "
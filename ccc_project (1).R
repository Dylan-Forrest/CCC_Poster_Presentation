##########################################
# Title: C3 Project
# Date: 9/7/2026
# Authors: Emma Lambert + Dylan Seay
##########################################

# Importing packages
library(tidyverse)
library(stats)
library(EnvStats)
library(ggplot2)
library(broom)
library(lme4)
library(lmerTest)
library(lsmeans)
library(psych)

# Setting data path
pathway <- "your/abcd/datapath/"

###############################
# Selecting and Cleaning Data
###############################

gen_path <- file.path(pathway,"abcd_general")

# Reading in Visit/Demographic tsv files
CCC_Data <- read_tsv(paste0(gen_path, '/ab_g_dyn.tsv'))
CCC_Data <-  CCC_Data %>%
  select(participant_id, 
         session_id, 
         ab_g_dyn__visit_age, # Participant age @ each session
         ab_g_dyn__design_mr__serial) # Scanner Serial #
         #ab_g_dyn__design_mr__model) # Scanner Model
age <- CCC_Data %>%
  select(participant_id, 
         session_id, 
         ab_g_dyn__visit_age)
CCC_Data <- filter(CCC_Data, session_id == "ses-00A") # Filtering for only Baseline data
CCC_Data <- CCC_Data %>%
  rename(age = ab_g_dyn__visit_age,
         scanner_id = ab_g_dyn__design_mr__serial)
         #scanner_model = ab_g_dyn__design_mr__model)

# Reading in Sex/Race tsv files
stc <- read_tsv(paste0(gen_path, '/ab_g_stc.tsv'))
stc <- stc %>%
  select(participant_id, 
         ab_g_stc__cohort_sex, # Participant sex
         ab_g_stc__cohort_ethnrace__leg, # Participant ethnicity 
         ab_g_stc__design_famrel) # Family ID 
CCC_Data <-  merge(CCC_Data, stc, by = "participant_id")

CCC_Data <- CCC_Data %>%
  rename(sex = ab_g_stc__cohort_sex,
         ethnrace = ab_g_stc__cohort_ethnrace__leg,
         family_id = ab_g_stc__design_famrel)
CCC_Data <- CCC_Data %>%
  mutate(sex = recode(sex,
                      '1' = 'Male',
                      '2' = 'Female'))
CCC_Data <- CCC_Data %>%
  mutate(ethnrace = recode(ethnrace, 
                           '1'= 'Hispanic',
                           '2' = 'White',
                           '3' = 'Black',
                           '4' = 'Asian',
                           '13' = 'Other'))

# Reading in SES tsv files
demo <- read_tsv(paste0(gen_path, '/ab_p_demo.tsv'))
demo <- demo %>%
  select(participant_id,
         session_id,
         ab_p_demo__income__hhold_001) # Total combined family income for the past 12 months
CCC_Data <- merge(CCC_Data, demo, by = c("participant_id", "session_id"), all.x = TRUE)
CCC_Data <- CCC_Data %>%
  rename(fam_income = ab_p_demo__income__hhold_001)
CCC_Data <- CCC_Data %>%
  mutate(fam_income = na_if(fam_income, 'n/a'),
         fam_income = na_if(fam_income, '999'),
         fam_income = na_if(fam_income, '777'))
CCC_Data <- CCC_Data %>%
  mutate(fam_income = factor(CCC_Data$fam_income, 
                             levels = c("1","2","3","4","5","6","7","8","9","10"), 
                             labels = c("Less than $5,000",
                                        "$5,000 through $11,999",
                                        "$12,000 through $15,999", 
                                        "$16,000 through $24,999", 
                                        "$25,000 through $34,999", 
                                        "$35,000 through $49,999", 
                                        "$50,000 through $74,999", 
                                        "$75,000 through $99,999", 
                                        "$100,000 through $199,999", 
                                        "$200,000 and greater"), ordered=TRUE))

####################################
# Reading in project-specific  files
####################################

############
# ACE DATA
############

#Setting ACE Path
ace_path <- "your/abcd/ACEs/datapath/"

# Reading in ACEs Data
ACEs <- read.csv(paste0(ace_path, "ACES_baseline_output.csv"),
                 na=c("NA","999", "777", "888", "n/a","N/A"))

# Convert baseling age from months to years #
ACEs_clean <- ACEs %>%
  mutate('baseline_age' = interview_age / 12)

# Keep Relevant Columns and Rename 5.0 -> 6.0 #
ACEs_clean <- ACEs_clean %>%
  select(src_subject_id, Abuse_ACEs, Neglect_ACEs, ACEs_Total) %>%
  mutate(src_subject_id = str_replace(src_subject_id, "NDAR_INV", "sub-")) %>%
  rename(participant_id = src_subject_id) 

ACEs_clean <- ACEs_clean %>%
  mutate(total_aces_yn = if_else(ACEs_Total > 0, 1, 0),
         abuse_yn = if_else(Abuse_ACEs > 0, 1, 0),
         neglect_yn = if_else(Neglect_ACEs > 0, 1, 0))


# Merge ACE Data
CCC_Data <- merge(CCC_Data, ACEs_clean, by = c("participant_id"), all.x = TRUE)


######################
# Hippocampal Volumes
######################

# Hippo Path
hippo_path <- "your/abcd/fressurfer/segmentation/datapath/"

# Hippocampal Volumes (Selecting, Renaming, Segmenting)
Hippo_vol <- read_csv(paste0(hippo_path, '/All_Hippo_Volumes.csv'))
Hippo_vol <- Hippo_vol %>%
  mutate(hippo_head = left_parasubiculum + 
           right_parasubiculum + 
           left_presubiculum_head + 
           right_presubiculum_head + 
           left_subiculum_head + 
           right_subiculum_head + 
           left_ca1_head + 
           right_ca1_head + 
           left_ca3_head + 
           right_ca3_head + 
           left_ca4_head + 
           right_ca4_head + 
           left_gc_ml_dg_head +
           right_gc_ml_dg_head +
           left_molecular_layer_hp_head + 
           right_molecular_layer_hp_head + 
           left_hata + 
           right_hata)

Hippo_vol <- Hippo_vol %>%
  mutate(hippo_body = left_presubiculum_body + 
           right_presubiculum_body + 
           left_subiculum_body + 
           right_subiculum_body + 
           left_ca1_body + 
           right_ca1_body + 
           left_ca3_body + 
           right_ca3_body + 
           left_ca4_body + 
           right_ca4_body + 
           left_gc_ml_dg_body +
           right_gc_ml_dg_body +
           left_molecular_layer_hp_body + 
           right_molecular_layer_hp_body + 
           left_fimbria + 
           right_fimbria)

Hippo_vol <- Hippo_vol %>%
  mutate(hippo_tail = left_hippocampal_tail + right_hippocampal_tail)

Hippo_vol <- Hippo_vol %>%
  mutate(hippo_fissure = left_hippocampal_fissure + right_hippocampal_fissure)

Hippo_vol <- Hippo_vol %>%
  mutate(whole_hippo = left_whole_hippocampus + right_whole_hippocampus)

# Drop all NAs 
Hippo_vol <- Hippo_vol %>%
  drop_na() 

# Selecting and merging all data to demo/ACEs
Hippo_vol <- Hippo_vol %>%
  select(participant_id, session_id, hippo_head, hippo_body, hippo_tail, whole_hippo) 

######################
# Amygdala Volumes
######################

# Setting amygdala path
amyg_path <- "your/abcd/fressurfer/segmentation/datapath/"

# Read in amygdala volumes
Amyg_vol <- read_csv(paste0(amyg_path, '/All_Amyg_Volumes.csv'))

# Amygdala Volumes (Selecting, Renaming, Segmenting)
Amyg_vol <- Amyg_vol %>%
  mutate(basolat_complex = `lh_Lateral-nucleus` +
           `lh_Basal-nucleus` +
           `lh_Accessory-Basal-nucleus` +
           `rh_Lateral-nucleus` +
           `rh_Basal-nucleus` +
           `rh_Accessory-Basal-nucleus`)

Amyg_vol <- Amyg_vol %>%
  mutate(centromed_complex = `lh_Central-nucleus` +
           `rh_Central-nucleus` +
           `lh_Medial-nucleus` +
           `rh_Medial-nucleus`)

Amyg_vol <- Amyg_vol %>%
  mutate(whole_amyg = lh_Whole_amygdala +
           rh_Whole_amygdala)

# Merge age data in
Amyg_vol <- merge(Amyg_vol, age, by = c("participant_id", "session_id"), all.x = TRUE)

# Drop NA to aid in slope/traj calculation
Amyg_vol <- Amyg_vol %>%
  drop_na() 

# Select and merge trajectory data to demo/ACEs
Amyg_vol <- Amyg_vol %>%
  select(participant_id, session_id, basolat_complex, centromed_complex, whole_amyg) 

amyg_hippo <- merge(Amyg_vol, Hippo_vol, by = c("participant_id", "session_id"), all.x = TRUE) 
# Before: 30,107 Scans
amyg_hippo <- drop_na(amyg_hippo)
# After: 29,934 Scans

CCC_Data <- merge(CCC_Data, amyg_hippo, by = c("participant_id", "session_id"), all.x = TRUE) %>%
  drop_na()


###########
# Analysis
###########


# Aim 1. Examine the relationship between ACE burden and type in hippocampal subfield volume trajectories. 

hippohead_aces <- lmer(hippo_head ~ ACEs_Total + age + sex +ethnrace + fam_income + (1|family_id:scanner_id), data = CCC_Data)
summary(hippohead_aces)

hippobody_aces <- lmer(hippo_body ~ ACEs_Total + age + sex +ethnrace + fam_income + (1|family_id:scanner_id), data = CCC_Data)
summary(hippobody_aces)

hippotail_aces <- lmer(hippo_tail ~ ACEs_Total + age + sex +ethnrace + fam_income + (1|family_id:scanner_id), data = CCC_Data)
summary(hippotail_aces)

hippohead_type <- lmer(hippo_head ~ Abuse_ACEs + Neglect_ACEs + age + sex +ethnrace + fam_income + (1|family_id:scanner_id), data = CCC_Data)
summary(hippohead_type)

hippobody_type <- lmer(hippo_body ~ Abuse_ACEs + Neglect_ACEs + age + sex +ethnrace + fam_income + (1|family_id:scanner_id), data = CCC_Data)
summary(hippobody_type)

hippotail_type <-lmer(hippo_tail ~ Abuse_ACEs + Neglect_ACEs + age + sex +ethnrace + fam_income + (1|family_id:scanner_id), data = CCC_Data)
summary(hippotail_type)

# Aim 2. Examine the relationship between ACE burden and type in amygdala subfield volume trajectories. 

basolat_aces <- lmer(basolat_complex ~ ACEs_Total + age + sex +ethnrace + fam_income + (1|family_id:scanner_id), data = CCC_Data)
summary(basolat_aces)

centromed_aces <- lmer(centromed_complex ~ ACEs_Total + age + sex +ethnrace + fam_income + (1|family_id:scanner_id), data = CCC_Data)
summary(centromed_aces)

basolat_type <- lmer(basolat_complex ~ Abuse_ACEs + Neglect_ACEs + age + sex +ethnrace + fam_income + (1|family_id:scanner_id), data = CCC_Data)
summary(basolat_type)

centromed_type <- lmer(centromed_complex ~ Abuse_ACEs + Neglect_ACEs + age + sex +ethnrace + fam_income + (1|family_id:scanner_id), data = CCC_Data)
summary(centromed_type)



# Get count of scans, find average
amyg_hippo <- amyg_hippo %>%
group_by(participant_id) %>%
  mutate(scan_count = count(session_id))


############################
# Demographics Information    
############################

describe(CCC_Data)    # Variable / Mean (SD) #
  # Age / 9.95 (0.62)
  # ACEs_Total / 1.19 (1.19)
  # Abuse / 0.45 (0.55)
  # Neglect / 0.14 (0.35)
  # basolateral / 2977.94 (308.11)
  # centromed / 157.50 (27.46)
  # hippo head / 3589.03 (387.43)
  # hippo body / 2525.70 (283.72)
  # hippo tail / 1241.32 (169.44)

table(CCC_Data$sex)
  # Male - 5423
  # Female - 4991

table(CCC_Data$ethnrace)
  # Asian - 190   
  # Black - 1457  
  # Hispanic - 2001    
  # Other - 1064    
  # White - 5702 


# Factoring Errors #
ACEs_clean$ACEs_Total <- as.factor(ACEs_clean$ACEs_Total)
ACEs_clean$Abuse_ACEs <- as.factor(ACEs_clean$Abuse_ACEs)
ACEs_clean$Neglect_ACEs <- as.factor(ACEs_clean$Neglect_ACEs)




###################
# Graphing Time
###################

# Creating sperate loong-form DF to plot ACEs and subfields
plotting_data <- CCC_Data %>%
  select(participant_id, 
         hippo_head, hippo_body, 
         hippo_tail, 
         basolat_complex, 
         centromed_complex, 
         total_aces_yn, 
         abuse_yn, 
         neglect_yn) %>%
  pivot_longer(
    cols = c(total_aces_yn, abuse_yn, neglect_yn),
    names_to = "ACE_Type",
    values_to = "Exposure") %>%
  mutate(Exposure = factor(
    Exposure,
    levels = c(0,1),
    labels = c("No", "Yes"))) # Making Exposure a factor, so it actually plottable



# Hippocampal Head - ACEs
ggplot(plotting_data,aes(
    x = ACE_Type,
    y = hippo_head)) +
  stat_summary( # Creating error bars
    aes(group = Exposure),
    fun.data = mean_se,
    geom = "errorbar",
    color = "black",
    width = 0.25,
    position = position_dodge(width = 0.4)) +
  stat_summary(
    aes(
      color = Exposure,
      group = Exposure),
    fun = mean,
    geom = "point",
    size = 10,
    position = position_dodge(width = 0.4)) +
  scale_color_manual(
    values = c(
      "No" = "darkgreen",
      "Yes" = "darkorange1")) +
  scale_x_discrete(
    labels = c(
      total_aces_yn = "Total ACEs", 
      abuse_yn = "Abuse ACEs",
      neglect_yn = "Neglect ACEs")) +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_blank(),
        legend.position = "none")



# Hippocampal Body - ACEs
ggplot(plotting_data,aes(
  x = ACE_Type,
  y = hippo_body)) +
  stat_summary( # Creating error bars
    aes(group = Exposure),
    fun.data = mean_se,
    geom = "errorbar",
    color = "black",
    width = 0.25,
    position = position_dodge(width = 0.4)) +
  stat_summary(
    aes(
      color = Exposure,
      group = Exposure),
    fun = mean,
    geom = "point",
    size = 10,
    position = position_dodge(width = 0.4)) +
  scale_color_manual(
    values = c(
      "No" = "darkgreen",
      "Yes" = "darkorange")) +
  scale_x_discrete(
    labels = c(
      total_aces_yn = "Total ACEs",
      neglect_yn = "Neglect ACEs")) +
  labs(
    x = "ACE Type",
    y = "Hippocampal Body Volume (mm³)",
    color = "ACE Exposure") +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_blank(),
        legend.position = "none")



# Hippocampal Tail - ACEs
ggplot(plotting_data,aes(
  x = ACE_Type,
  y = hippo_tail)) +
  stat_summary( # Creating error bars
    aes(group = Exposure),
    fun.data = mean_se,
    geom = "errorbar",
    color = "black",
    width = 0.25,
    position = position_dodge(width = 0.4)) +
  stat_summary(
    aes(
      color = Exposure,
      group = Exposure),
    fun = mean,
    geom = "point",
    size = 10,
    position = position_dodge(width = 0.4)) +
  scale_color_manual(
    values = c(
      "No" = "darkgreen",
      "Yes" = "darkorange")) +
  scale_x_discrete(
    labels = c(
      total_aces_yn = "Total ACEs",
      abuse_yn = "Abuse ACEs",
      neglect_yn = "Neglect ACEs")) +
  labs(
    x = "ACE Type",
    y = "Hippocampal Tail Volume (mm³)",
    color = "ACE Exposure") +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_blank(),
        legend.position = "none")




# Basolateral Complex - ACEs
ggplot(plotting_data,aes(
  x = ACE_Type,
  y = basolat_complex)) +
  stat_summary( # Creating error bars
    aes(group = Exposure),
    fun.data = mean_se,
    geom = "errorbar",
    color = "black",
    width = 0.25,
    position = position_dodge(width = 0.4)) +
  stat_summary(
    aes(
      color = Exposure,
      group = Exposure),
    fun = mean,
    geom = "point",
    size = 10,
    position = position_dodge(width = 0.4)) +
  scale_color_manual(
    values = c(
      "No" = "darkgreen",
      "Yes" = "darkorange")) +
  scale_x_discrete(
    labels = c(
      total_aces_yn = "Total ACEs",
      abuse_yn = "Abuse ACEs",
      neglect_yn = "Neglect ACEs")) +
  labs(
    x = "ACE Type",
    y = "Basolateral Amygdala Volume (mm³)",
    color = "ACE Exposure") +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_blank(),
        legend.position = "none")



# Centromedial Complex - ACEs
ggplot(plotting_data,aes(
  x = ACE_Type,
  y = centromed_complex)) +
  stat_summary( # Creating error bars
    aes(group = Exposure),
    fun.data = mean_se,
    geom = "errorbar",
    color = "black",
    width = 0.15,
    position = position_dodge(width = 0.4)) +
  stat_summary(
    aes(
      color = Exposure,
      group = Exposure),
    fun = mean,
    geom = "point",
    size = 3,
    position = position_dodge(width = 0.4)) +
  scale_color_manual(
    values = c(
      "No" = "darkgreen",
      "Yes" = "darkorange")) +
  scale_x_discrete(
    labels = c(
      total_aces_yn = "Total ACEs",
      abuse_yn = "Abuse ACEs",
      neglect_yn = "Neglect ACEs")) +
  ylim(166.1,173) +
  labs(
    x = "ACE Type",
    y = "Centromedial Amygdala Volume (mm³)",
    color = " ACE Exposure") +
    ggtitle("Associations Between ACE Type and Centromedial Amygdala Volume")

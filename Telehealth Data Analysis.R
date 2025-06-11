###

# Telehealth Data
# Kevin Shi


# There are duplicate encounters?? Need to ask Aris/Salman.

###

# Setup ---------------
rm(list=ls())
setwd("~/WIP Research Stuff/Telehealth and Nephrology/Data and Analysis 2/5_31 data")
library(tidyverse)
library(lme4)       # For GLMM with Poisson
library(glmmTMB)    # For GLMM with negative binomial (if overdispersion is present)
library(performance) # For model diagnostics
library(lubridate)
library(DHARMa)
library(broom.mixed)
library(MatchIt)
library(tableone)

# Reading data ---------------
# this data is from 1/2017 to 12/2024
main <- read.csv("telehealth_main_06_04_2025.csv")
labs <- read.csv("telehealth_labs_5_31_2025.csv")
meds <- read.csv("telehealth_current_meds_5_31_2025.csv")
orders <- read.csv("telehealth_orders_5_31_2025.csv")
notes <- read.csv("telehealth_notes_5_31_2025.csv")
med_changes <- read.csv("telehealth_med_changes_5_31_2025.csv")

# Preprocessing ---------------
main$ENCOUNTER_TYPE <- as.factor(main$ENCOUNTER_TYPE)
main$APPT_DATE <- as.Date(main$APPT_DATE, format = "%m/%d/%Y")
main$AGE_AT_ENCOUNTER <- as.numeric(main$AGE_AT_ENCOUNTER)
main$BMI_WITHIN_6_MO <- as.numeric(main$BMI_WITHIN_6_MO)
main$smoking_status <- toupper(main$smoking_status)
colnames(main)[1] <- "PAT_ENC_CSN_ID"
colnames(meds)[1] <- "PAT_ENC_CSN_ID"
colnames(orders)[1] <- "PAT_ENC_CSN_ID"
colnames(med_changes)[1] <- "PAT_ENC_CSN_ID"

dim(main)

# remove duplicates
main <- main %>%
  distinct(PAT_ENC_CSN_ID, .keep_all = TRUE)

# remove not completed/no-show encounters
main <- main %>% filter(ENCOUNTER_STATUS == "Completed")

levels(main$ENCOUNTER_TYPE)

# remove some encounter types we don't want for now

main <- main %>%
  filter(!(ENCOUNTER_TYPE %in% c("AMBULATORY BP 24HR MONITOR", "PD CLINIC", "PHONE CONSULT", "SCHEDULED TELEPHONE FOLLOW UP", "TELEPHONE"))) %>%
  mutate(ENCOUNTER_TYPE = droplevels(ENCOUNTER_TYPE))

levels(main$ENCOUNTER_TYPE)

dim(main) # how many encounters do we have

# Add a telehealth column
main <- main %>% mutate(Telehealth = ifelse(ENCOUNTER_TYPE %in% c("IVV FUP", "IVV NEW", 
                                                                    "VIDEO VISIT", "VIDEO VISIT 15", "VIDEO VISIT 20", 
                                                                    "VIDEO VISIT 30", "VIDEO VISIT 60", 
                                                                    "VIDEO VISIT NEW"), "Telehealth Video", "In-person"))
# more pre-processing
main$Telehealth <- factor(main$Telehealth)
main$SEX_ASSIGNED_AT_BIRTH <- factor(main$SEX_ASSIGNED_AT_BIRTH)
main$RACE_ETHNICITY <- factor(main$RACE_ETHNICITY)

# Calculating eGFR ----------
# Ideally would like eGFR

labs_CSNs <- labs$PAT_ENC_CSN_ID # 53728

main_CSNs <- main$PAT_ENC_CSN_ID # 51882

# every encounter has labs - yay, but some values are NULL
common_values <- intersect(labs_CSNs, main_CSNs) # 51882 - good
unique_to_labs <- setdiff(labs_CSNs, main_CSNs) # 1846
unique_to_main <- setdiff(main_CSNs, labs_CSNs) # empty - good

crs <- labs %>% select(PAT_ENC_CSN_ID,Creatinine_ORD_VALUE,Creatinine_LAB_REFERENCE_UNIT)

# so it looks like the CSNs map directly to each other, which is nice, but there is probably some lookback logic that I don't know
# should ask Salman and Aris

egfr_df <- main %>% select(PAT_ENC_CSN_ID,AGE_AT_ENCOUNTER,SEX_ASSIGNED_AT_BIRTH,Telehealth)

egfr_df <- left_join(egfr_df, crs , by = "PAT_ENC_CSN_ID")

# check if anything needs to be done for units - looks good, all mg/dL or NULL
# units <- as.factor(egfr_df$Creatinine_LAB_REFERENCE_UNIT)
# levels(units)

# eGFR = 142*  min(standardized Scr/K, 1)^?? *  max(standardized Scr/K, 1)^-1.200 *  0.9938^Age *  1.012 [if female]
# eGFR (estimated glomerular filtration rate) = mL/min/ 1.73 m2
# Scr (serum creatinine) = mg/dL
# K = 0.7 (females) or 0.9 (males)
# ?? = -0.241 (females) or -0.302 (males)
# min = indicates the minimum of Scr/K or 1
# max = indicates the maximum of Scr/K or 1

# requires female = 0 if male, 1 if female
egfr_2021 <- function(age, female, cr) {
  return(142 * min(cr / (0.9-0.2*female), 1)^(-0.302+0.061*female) *
           max(cr/(0.9-0.2*female), 1)^(-1.2) * 
           0.9938^age *
           1.012^female)
}

unknowns <- which(main$SEX_ASSIGNED_AT_BIRTH != "Female" & main$SEX_ASSIGNED_AT_BIRTH != "Male")
main[unknowns,]$GENDER_IDENTITY

# some of the unknown genders identify as male, some as female

# Actually all the unknowns identify as males?? Probably best to just go by sex_assigned_at_birth for everything, despite NAs
egfr_df$Creatinine_ORD_VALUE <- as.numeric(egfr_df$Creatinine_ORD_VALUE)
egfr_df <- egfr_df %>% mutate(female = ifelse(SEX_ASSIGNED_AT_BIRTH == "Female", 1, 
                                              ifelse(SEX_ASSIGNED_AT_BIRTH == "Male",0,NA)))
class(egfr_df$female)
class(egfr_df$Creatinine_ORD_VALUE)
class(egfr_df$AGE_AT_ENCOUNTER)
egfr_df <- egfr_df %>%
  mutate(egfr = if_else(!is.na(female) & !is.na(Creatinine_ORD_VALUE),
                        mapply(egfr_2021, AGE_AT_ENCOUNTER, female, Creatinine_ORD_VALUE),
                        NA_real_))

which(is.na(egfr_df$female)) # 19 encs w/ unknown sex at birth

head(egfr_df)

to_join <- egfr_df %>% select(PAT_ENC_CSN_ID,egfr)
dim(to_join)
dim(main)

main <- main %>% left_join(to_join, by = "PAT_ENC_CSN_ID")

length(which(is.na(main$egfr))) # missing 23084 eGFRs

# Med counts ----------

meds_count <- meds %>%
  count(PAT_ENC_CSN_ID, name = "meds_count")

head(meds_count)

dim(meds_count)
dim(main)

meds_CSNs <- meds_count$PAT_ENC_CSN_ID

main_CSNs <- main$PAT_ENC_CSN_ID

common_values <- intersect(meds_CSNs, main_CSNs) # 40851
unique_to_meds <- setdiff(meds_CSNs, main_CSNs) # 970
unique_to_main <- setdiff(main_CSNs, meds_CSNs) # 11031

main <- main %>% left_join(meds_count, by = "PAT_ENC_CSN_ID")

# can we consider all the NAs to be 0? I.e. if there weren't meds listed, then the patients were taking no meds?
length(which(is.na(main$meds_count)))

main <- main %>%
  mutate(meds_count = replace(meds_count, is.na(meds_count), 0))

head(main$meds_count)

# Orders in encounters ---------------
orders$ORDER_TYPE <- as.factor(orders$ORDER_TYPE)
levels(orders$ORDER_TYPE)

# let's remove the order types that I don't think are relevant

orders[which(orders$ORDER_TYPE == "PFT"),]

# Consult
# Cardiac Services
# E-Consult
# ECG
# Echocardiography
# Immunization/Injection
# Outpatient Referral
# Imaging
# Imaging-Vasc
# OT
# PFT
# Pathology and Cytology
# Procedures
# Sleep Center
# Transcribed Referral
# Equipment
# Code Status
# Neurology

orders <- orders %>% filter(ORDER_TYPE %in% c("Consult", "Cardiac Services", "E-Consult", "ECG", 
                           "Echocardiography", "Immunization/Injection", 
                           "Outpatient Referral", "Imaging", "Imaging-Vasc", 
                           "OT", "PFT", "Pathology and Cytology", "Procedures", 
                           "Sleep Center", "Transcribed Referral", "Equipment", 
                           "Code Status", "Neurology"))


dim(orders)
dim(main)

orders_CSNs <- orders$PAT_ENC_CSN_ID

main_CSNs <- main$PAT_ENC_CSN_ID

common_values <- intersect(orders_CSNs, main_CSNs) # 5958 encounters have orders
unique_to_orders <- setdiff(orders_CSNs, main_CSNs) # 120
unique_to_main <- setdiff(main_CSNs, orders_CSNs) # 45924

orders_count <- orders %>%
  count(PAT_ENC_CSN_ID, name = "orders_count")

main <- main %>% left_join(orders_count, by = "PAT_ENC_CSN_ID")

# can we consider all the NAs to be 0? I.e. if there weren't orders listed, then there were no orders placed?
length(which(is.na(main$orders_count))) # 45924 encounters without order of the type above

main <- main %>%
  mutate(orders_count = replace(orders_count, is.na(orders_count), 0))

# Med changes --------------------

# let's remove all the changes that have to deal with historical medications - I think these are patient reported ones
med_changes$ORIGINAL_MED_CLASS <- as.factor(med_changes$ORIGINAL_MED_CLASS)
med_changes$NEW_MED_CLASS <- as.factor(med_changes$NEW_MED_CLASS)
med_changes$STATUS <- as.factor(med_changes$STATUS)
levels(med_changes$NEW_MED_CLASS)
levels(med_changes$ORIGINAL_MED_CLASS)

med_changes <- med_changes %>%
  filter(!(ORIGINAL_MED_CLASS %in% c("Historical Med", "Historical"))) %>%
  mutate(ORIGINAL_MED_CLASS = droplevels(ORIGINAL_MED_CLASS))

med_changes <- med_changes %>%
  filter(!(NEW_MED_CLASS %in% c("Historical Med"))) %>%
  mutate(NEW_MED_CLASS = droplevels(NEW_MED_CLASS))

med_changes <- med_changes %>%
  filter(!(STATUS %in% c("REORDERED"))) %>%
  mutate(STATUS = droplevels(STATUS))

levels(med_changes$NEW_MED_CLASS)
levels(med_changes$ORIGINAL_MED_CLASS)
levels(med_changes$STATUS)

med_changes_count <- med_changes %>%
  count(PAT_ENC_CSN_ID, name = "med_changes_count")

main <- main %>% left_join(med_changes_count, by = "PAT_ENC_CSN_ID")

# if med changes counts doesn't exist, set to 0
main <- main %>%
  mutate(med_changes_count = replace(med_changes_count, is.na(med_changes_count), 0))

# New vs. follow-up appointments --------------

levels(main$ENCOUNTER_TYPE)

# 1 = new visit, 0 = f/u (or at least not explicitly labeled as new)
main$New_Old <- ifelse(grepl("NEW", main$ENCOUNTER_TYPE, ignore.case = TRUE), 1, 0)

# Time since last visit -----------------------
main <- main %>%
  arrange(UCSF_MRN, APPT_DATE) %>%
  group_by(UCSF_MRN) %>%
  mutate(Time_since_last_visit = as.numeric(APPT_DATE - lag(APPT_DATE)),
         Time_since_last_visit = ifelse(is.na(Time_since_last_visit), 0, Time_since_last_visit)) %>%
  ungroup()

# Add notes -----------------------
main <- main %>% left_join(notes, by = "PAT_ENC_CSN_ID")

# Plot of telehealth encounters per year ----------

# First let's find how many telehealth vs. in-person encounters per year
basic <- main %>% select(APPT_DATE, ENCOUNTER_TYPE, Telehealth)

basic_summary <- basic %>%
  mutate(Year = format(as.Date(APPT_DATE), "%Y")) %>%
  group_by(Year, Telehealth) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = Telehealth, values_from = Count, values_fill = 0)

# table of telehealth visits per year (absolute count)
print(basic_summary)

basic_summary$Total <- basic_summary$`In-person` + basic_summary$`Telehealth Video`

basic_summary <- basic_summary %>%  mutate(Percent_Telehealth = (`Telehealth Video` / Total) * 100)

basic_long <- basic_summary %>%
  pivot_longer(cols = c("Telehealth Video", "In-person"), names_to = "Visit_Type", values_to = "Count")

# Plot using ggplot2
ggplot(basic_long, aes(x = Year, y = Count, fill = Visit_Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("In-person" = "lightblue", "Telehealth Video" = "lightpink")) +
  geom_text(data = basic_summary, 
            aes(x = Year, y = Total, label = paste0(round(Percent_Telehealth, 1), "%")), 
            vjust = -0.5,
            inherit.aes = FALSE) +
  labs(x = "Year", y = "Encounter count", fill = "Visit Type", title = "Telehealth Encounters per Year") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# when was the earliest telehealth appointment?
# telehealth_dates <- basic %>% filter(Telehealth == "Telehealth")
# in_person_dates <- basic %>% filter(Telehealth == "In-person")

# min(in_person_dates$APPT_DATE) # 2017-01-03
# min(telehealth_dates$APPT_DATE) # 2018-10-31
# telehealth_dates[which(telehealth_dates$APPT_DATE == min(telehealth_dates$APPT_DATE)),] # was a video visit on Halloween

rm(basic, basic_long, basic_summary)

main2020on <- main %>% filter(year(APPT_DATE) == 2020 | 
                                year(APPT_DATE) == 2021 |
                                year(APPT_DATE) == 2022 |
                                year(APPT_DATE) == 2023 |
                                year(APPT_DATE) == 2024)

table(main2020on$Telehealth)

main3 <- main # dummy for original main (2017-2024, completed encs)

main <- main2020on

# Propensity matching --------------------
main <- main %>%
  mutate(telehealth_binary = ifelse(Telehealth == "Telehealth Video", 1, 0))

# Check distribution
table(main$telehealth_binary)

main2 <- main

# try to impute NAs based on patient medians
main2 <- main2 %>%
  group_by(UCSF_MRN) %>%
  mutate(
    egfr = ifelse(is.na(egfr), median(egfr, na.rm = TRUE), egfr)
  ) %>%
  ungroup()

# still with NAs

length(which(is.na(main2$egfr)==TRUE)) # 10875

# drop the rest

main_clean <- main2 %>%
  filter(!is.na(egfr))

# Fit propensity score model and perform matching
m.out <- matchit(telehealth_binary ~ AGE_AT_ENCOUNTER + SEX_ASSIGNED_AT_BIRTH + 
                   RACE_ETHNICITY + egfr + meds_count + New_Old, 
                 data = main_clean, method = "full", distance = "logit")

# Extract matched data if needed
matched_data <- match.data(m.out)

# Before matching
# Create a dataframe for plotting
ps_df <- data.frame(
  distance = m.out$distance,
  telehealth_binary = main_clean$telehealth_binary
)

# Plot
ggplot(ps_df, aes(x = distance, fill = as.factor(telehealth_binary))) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = c("blue", "red"), labels = c("In-person", "Telehealth Video")) +
  labs(x = "Propensity Score", fill = "Visit Type", title = "Propensity Score Distribution (Before Matching)") +
  theme_minimal()

# After matching
ggplot(data = matched_data, aes(x = distance, fill = as.factor(telehealth_binary))) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = c("blue", "red"), labels = c("In-person", "Telehealth Video")) +
  labs(x = "Propensity Score", fill = "Visit Type", title = "Propensity Score Distribution (After Matching)") +
  theme_minimal()

# After matching - WEIGHTED plot
ggplot(data = matched_data, aes(x = distance, fill = as.factor(telehealth_binary), weight = weights)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = c("blue", "red"), labels = c("In-person", "Telehealth Video")) +
  labs(x = "Propensity Score", fill = "Visit Type", title = "Propensity Score Distribution (After Matching — Weighted)") +
  theme_minimal()

# Table 1 ---------------
# Now let's work on a "Table 1" for telehealth vs. in-person visits
# We are separating things by _encounter type_ (telehealth vs. in-person) - there will be patient duplicates/cross-over

# This package doesn't work for the version of R on the RAE??

df.t1 <- main
# # 1. LIST the variables that will be in table 1.
t1vars <- c("Telehealth", "AGE_AT_ENCOUNTER", "SEX_ASSIGNED_AT_BIRTH", "RACE_ETHNICITY", "egfr",
             "BMI_WITHIN_6_MO", "MYC_ENROLLMENT_STATUS", "INSURANCE_GROUP", "smoking_status", "PMHX_HTN",
             "PMHX_DM","PMHX_CAD","PMHX_CVA","PMHX_CHF", "CPT_CHARGE_CODE", "meds_count", "orders_count", "med_changes_count", "New_Old", "Time_since_last_visit")
# # 2. WHICH are categorical?
t1catvars <- c("Telehealth", "SEX_ASSIGNED_AT_BIRTH", "RACE_ETHNICITY",
                "MYC_ENROLLMENT_STATUS", "INSURANCE_GROUP", "smoking_status", "PMHX_HTN",
                "PMHX_DM","PMHX_CAD","PMHX_CVA","PMHX_CHF", "CPT_CHARGE_CODE", "New_Old")
# # 3. WHICH are non-normal?
# # t1nnormvars <- c("egfr", "baselineUACR_mgg", "kfre1", "kfre2", "kfre5", "eventMonths")
# # none for now
table1 <- CreateTableOne(vars = t1vars, data = df.t1, factorVars = t1catvars, strata = "Telehealth") # gives statistical tests (check documentation)
print(table1, quote = FALSE, contDigits = 5)

df.t1 <- df.t1[df.t1$BMI_WITHIN_6_MO <= 100, ] # remove all the rows with BMI > 100
table1 <- CreateTableOne(vars = t1vars, data = df.t1, factorVars = t1catvars, strata = "Telehealth") # gives statistical tests (check documentation)
print(table1, quote = FALSE, contDigits = 5)

# Well time to start from scratch I guess

# combine_counts_and_proportions <- function(tab, prop) {
#   result <- data.frame(matrix(ncol = ncol(tab), nrow = nrow(tab)))
#   colnames(result) <- colnames(tab)
#   rownames(result) <- rownames(tab)
#   
#   for (col in colnames(tab)) {
#     result[[col]] <- paste(tab[, col], " (", round(prop[, col], 1), ")", sep = "")
#   }
#   
#   return(result)
# }

# make_table_1_rough <- function(df.t1, list_variables, list_categorical_vars) {
#   t1vars <- list_variables
#   t1catvars <- list_categorical_vars
#   
#   # Create the summary table for categorical variables (count and percent)
#   table_cat <- lapply(1:length(t1catvars), function(i) {
#     var <- t1catvars[i]
#     tab <- table(df.t1[[var]], df.t1$Telehealth)  # Cross-tabulate with Telehealth strata
#     prop <- prop.table(tab, 2) * 100  # Calculate percentages within each group
#     res <- combine_counts_and_proportions(tab = tab, prop = prop)  # Combine count and percent
#     colnames(res) <- c("In-person", "Telehealth")  # Label columns
#     # Add variable name as a label to the result
#     rownames(res) <- paste(var, rownames(res), sep = ": ")  
#     return(res)
#   })
#   
#   cat("Categorical Variables Summary:\n")
#   for(i in 1:length(table_cat)) {
#     var_name <- t1catvars[i]  # Get variable name
#     cat("\n", var_name, ":\n")  # Print variable name as a header
#     print(table_cat[[i]])
#   }
#   
#   # Summary for continuous variables (means and SDs)
#   t1contvars <- setdiff(t1vars, t1catvars)  # Identify continuous variables
#   table_cont <- lapply(t1contvars, function(var) {
#     means <- tapply(df.t1[[var]], df.t1$Telehealth, mean, na.rm = TRUE)
#     sds <- tapply(df.t1[[var]], df.t1$Telehealth, sd, na.rm = TRUE)
#     
#     res <- data.frame(
#       Mean_In_Person = round(means[1], 2),
#       Mean_Telehealth = round(means[2], 2),
#       SD_In_Person = round(sds[1], 2),
#       SD_Telehealth = round(sds[2], 2)
#     )
#     return(res)
#   })
#   
#   cat("Continuous Variables Summary:\n")
#   for(i in 1:length(table_cont)) {
#     var_name <- t1contvars[i]  # Get variable name
#     cat("\n", var_name, ":\n")  # Print variable name as a header
#     print(table_cont[[i]])
#   }
  
  # Optional: Statistical tests (Chi-squared for categorical variables)
  # You can uncomment this block if you want statistical tests
  # cat("\nStatistical Tests (Chi-squared for categorical variables):\n")
  # for (var in t1catvars) {
  #   test <- chisq.test(table(df.t1[[var]], df.t1$Telehealth))
  #   cat("\nTest for", var, ":\n")
  #   print(test)
  # }
}

df.t1 <- main # with all data

# 1. List the variables that will be in table 1.
t1vars <- c("Telehealth", "AGE_AT_ENCOUNTER", "SEX_ASSIGNED_AT_BIRTH", "RACE_ETHNICITY",
            "BMI_WITHIN_6_MO", "MYC_ENROLLMENT_STATUS", "smoking_status", "PMHX_HTN",
            "PMHX_DM", "PMHX_CAD", "PMHX_CVA", "PMHX_CHF", "CPT_CHARGE_CODE", "INSURANCE_GROUP", 
            "egfr", "meds_count", "orders_count", "med_changes_count")

# 2. Which are categorical?
t1catvars <- c("Telehealth", "SEX_ASSIGNED_AT_BIRTH", "RACE_ETHNICITY",
               "MYC_ENROLLMENT_STATUS", "smoking_status", "PMHX_HTN",
               "PMHX_DM", "PMHX_CAD", "PMHX_CVA", "PMHX_CHF", "CPT_CHARGE_CODE", "INSURANCE_GROUP")
# 3. WHICH are non-normal?
# t1nnormvars <- c("egfr", "baselineUACR_mgg", "kfre1", "kfre2", "kfre5", "eventMonths")
table1 <- CreateTableOne(vars = t1vars, data = df.t1, factorVars = t1catvars, strata = "Telehealth") # gives statistical tests (check documentation)
print(table1, quote = FALSE, contDigits = 5)


# make_table_1_rough(df.t1,t1vars,t1catvars)
# 
# # since telehealth wasn't really an option over this entire dataframe, can limit data to certain time
# df.t1 <-  df.t1 %>% filter(year(APPT_DATE) == 2021 | year(APPT_DATE) == 2020) # with only 2020-2021 data
# 
# make_table_1_rough(df.t1,t1vars,t1catvars)

# are eGFRs different between groups?

tele <- df.t1 %>% filter(Telehealth == "Telehealth Video")
inp <- df.t1 %>% filter(Telehealth == "In-person")

# a bunch of ad-hoc tests
# t.test(tele$AGE_AT_ENCOUNTER,inp$AGE_AT_ENCOUNTER)
# t.test(tele$BMI_WITHIN_6_MO,inp$BMI_WITHIN_6_MO)
# t.test(tele$egfr, inp$egfr)
# t.test(tele$meds_count, inp$meds_count)
# t.test(tele$orders_count,inp$orders_count)
# t.test(tele$med_changes_count, inp$med_changes_count)
# 
# cont <- table(df.t1$INSURANCE_GROUP, df.t1$Telehealth)
# cont
# chisq.test(cont)

# Other population descriptions ------------------

encounter_count <- main %>%
  group_by(UCSF_MRN) %>%
  summarise(num_encounters = n_distinct(PAT_ENC_CSN_ID)) %>%
  ungroup()

# View the result
print(encounter_count)

hist(encounter_count$num_encounters)
mean(encounter_count$num_encounters)

encounter_count_2020_2021 <- main2020_2021 %>%
  group_by(UCSF_MRN) %>%
  summarise(num_encounters = n_distinct(PAT_ENC_CSN_ID)) %>%
  ungroup()

hist(encounter_count_2020_2021$num_encounters)
mean(encounter_count_2020_2021$num_encounters)

length(unique(main$PROV_NAME)) # number of providers

length(unique(main$UCSF_MRN)) # number of patients

# Build a model -------------------

# poisson_model <- glmer(med_changes_count ~ Telehealth + SEX_ASSIGNED_AT_BIRTH +
#                          AGE_AT_ENCOUNTER + RACE_ETHNICITY + meds_count +
#                          (1 | UCSF_MRN),
#                        data = main,
#                        family = poisson)
# 
# # Check overdispersion
# performance::check_overdispersion(poisson_model)

# pre-process

main$PMHX_DM <- as.factor(main$PMHX_DM)
main$PMHX_HTN <- as.factor(main$PMHX_HTN)
main$AGE_AT_ENCOUNTER_s <- scale(main$AGE_AT_ENCOUNTER) #normalize continuous variables
main$meds_count_s <- scale(main$meds_count)
main$egfr_s <- scale(main$egfr)
main$Time_since_last_visit_s <- scale(main$Time_since_last_visit)
main$years_since_2020 <- as.numeric(format(main$APPT_DATE, "%Y")) - 2020 # does this do anything?

levels(main$RACE_ETHNICITY) # first level is Asian
main$RACE_ETHNICITY <- relevel(main$RACE_ETHNICITY, ref = "White") #relabel to have "White" as ref

# optionally recode race/eth so we only have White, Latinx, Black, and Asian (everything else will be unknown/other)
levels(main$RACE_ETHNICITY)
main <- main %>%
  mutate(RACE_ETHNICITY2 = recode(RACE_ETHNICITY,
                                 "White" = "White",
                                 "Asian" = "Asian",
                                 "Black or African American" = "Black or African American",
                                 "Latinx" = "Latinx",
                                 "Multi-Race/Ethnicity" = "Unknown/Other",
                                 "Native American or Alaska Native" = "Unknown/Other",
                                 "Native Hawaiian or Other Pacific Islander" = "Unknown/Other",
                                 "Other" = "Unknown/Other",
                                 "Southwest Asian and North African" = "Unknown/Other",
                                 "Unknown/Declined" = "Unknown/Other"
  ))

nb_model <- glmmTMB(med_changes_count ~ Telehealth + SEX_ASSIGNED_AT_BIRTH + 
                      AGE_AT_ENCOUNTER_s + RACE_ETHNICITY2 + New_Old + Time_since_last_visit_s
                    + meds_count_s + egfr_s +
                      PMHX_DM + PMHX_HTN +
                      (1 | UCSF_MRN),
                    data = main,
                    family = nbinom2)

nb_model_slope <- glmmTMB(
  med_changes_count ~ Telehealth + SEX_ASSIGNED_AT_BIRTH + 
    AGE_AT_ENCOUNTER_s + RACE_ETHNICITY2 + New_Old + Time_since_last_visit_s + 
    + meds_count_s + egfr_s + PMHX_DM + PMHX_HTN +
    (1 + Telehealth | UCSF_MRN),
  data = main,
  family = nbinom2
)

nb_model_hierarchical <- glmmTMB(
  med_changes_count ~ Telehealth + SEX_ASSIGNED_AT_BIRTH + 
    AGE_AT_ENCOUNTER_s + RACE_ETHNICITY2 + New_Old + 
    Time_since_last_visit_s + meds_count_s + egfr_s + PMHX_DM + PMHX_HTN +
    (1 + Telehealth | UCSF_MRN) +
    (1 | ENCOUNTER_NOTE_AUTHOR),
  data = main,
  family = nbinom2
)

summary(nb_model_hierarchical)
summary(nb_model_hierarchical_full)

nb_model_hierarchical_full <- glmmTMB(
  med_changes_count ~ Telehealth + SEX_ASSIGNED_AT_BIRTH + 
    AGE_AT_ENCOUNTER_s + RACE_ETHNICITY2 + New_Old + 
    Time_since_last_visit_s + meds_count_s + egfr_s + PMHX_DM + PMHX_HTN +
    (1 + Telehealth | UCSF_MRN) +
    (1 + Telehealth | ENCOUNTER_NOTE_AUTHOR),
  data = main,
  family = nbinom2
)

anova(nb_model_hierarchical, nb_model_hierarchical_full)

nb_model_hierarchical <- glmmTMB(
  med_changes_count ~ Telehealth + SEX_ASSIGNED_AT_BIRTH + 
    AGE_AT_ENCOUNTER_s + RACE_ETHNICITY2 + New_Old +
    Time_since_last_visit_s + meds_count_s + egfr_s + PMHX_DM + PMHX_HTN +
    (1 + Telehealth | UCSF_MRN) +
    (1 | ENCOUNTER_NOTE_AUTHOR),
  data = main,
  family = nbinom2
)

# looks like the best model has random patient slopes and intercepts, but only random provider intercepts

library(DHARMa)

# Simulate residuals
sim_res <- simulateResiduals(nb_model_hierarchical, n = 1000)
model_data <- nb_model_hierarchical$frame  # extract cleaned data
plotResiduals(sim_res, form = model_data$AGE_AT_ENCOUNTER_s)
plotResiduals(sim_res, form = model_data$meds_count_s)
plotResiduals(sim_res, form = model_data$Time_since_last_visit_s)
plotResiduals(sim_res, form = model_data$egfr_s)

plotResiduals(sim_res, form = predict(nb_model_hierarchical, type = "response"))
plotResiduals(sim_res, form = predict(nb_model_hierarchical, type = "link"))

testDispersion(sim_res)
testZeroInflation(sim_res)
testUniformity(sim_res)
testOutliers(sim_res)

hist(sim_res$scaledResiduals)
qqnorm(sim_res$scaledResiduals)
qqline(sim_res$scaledResiduals)

# Tidy up the model

tidy_nb_slope <- broom.mixed::tidy(nb_model_hierarchical, effects = "fixed", conf.int = TRUE, exponentiate = TRUE)

irr_table_slope <- broom.mixed::tidy(nb_model_hierarchical) %>%
  filter(effect == "fixed") %>%
  mutate(
    IRR = exp(estimate),
    lower = exp(estimate - 1.96 * std.error),
    upper = exp(estimate + 1.96 * std.error)
  )

irr_table_slope$significance <- case_when(
  irr_table_slope$p.value < 0.001 ~ "***",
  irr_table_slope$p.value < 0.01  ~ "**",
  irr_table_slope$p.value < 0.05  ~ "*",
  TRUE ~ ""
)

# Print
# print(irr_table, digits = 3)
# 
# # Add a column for color based on significance
# irr_table$color <- ifelse(irr_table$p.value < 0.05, "red", "grey")

# Add a column for color based on significance
irr_table_slope$color <- ifelse(irr_table_slope$p.value < 0.05, "red", "grey")

print(irr_table_slope, digits = 3)

# Plot the IRR with significance stars and color highlighting
# ggplot(irr_table, aes(x = reorder(term, IRR), y = IRR, color = color)) +
#   geom_point(size = 3) +
#   geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
#   coord_flip() +
#   theme_minimal() +
#   labs(
#     x = "Predictor",
#     y = "Incidence Rate Ratio (IRR)",
#     title = "IRR (with 95% CI) for Medication Changes"
#   ) +
#   geom_text(aes(label = significance), hjust = -0.2, size = 4) +
#   scale_color_identity()  # Use the color column directly

irr_table_slope$term

# new term names to be prettier

irr_table_slope$term <- c("Intercept", "Telehealth", "Male Sex", "Age", 
                          "Asian Race", "Black or African American Race",
                          "Latinx Race", "Unknown/Other Race", 
                          "First Visit", "Time Since Last Visit", 
                          "Existing Medication Count", "Estimated Glomerular Filtration Rate",
                          "History of Diabetes", "History of Hypertension")

irr_table_slope <- irr_table_slope %>%
  filter(term != "Intercept")

print(irr_table_slope, digits = 3)

# Plot the IRR with significance stars and color highlighting
ggplot(irr_table_slope, aes(x = reorder(term, IRR), y = IRR, color = color)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() +
  theme(
    plot.margin = margin(t = 10, r = 30, b = 10, l = 10),
    plot.title = element_text(hjust = 0.5)  # optional: center the title
  ) + 
  theme_minimal() +
  labs(
    x = "Predictor",
    y = "Incidence Rate Ratio",
    title = "Incidence Rate Ratio\nfor Medication Changes"
  ) +
  #geom_text(aes(label = significance), hjust = -0.2, size = 4) +
  scale_color_identity()  # Use the color column directly

ggplot(irr_table_slope, aes(x = term, y = IRR, color = color)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_minimal() +
  labs(
    x = "Predictor",
    y = "Incidence Rate Ratio (IRR)",
    title = "IRR (with 95% CI) for Medication Changes"
  ) +
  geom_text(aes(label = significance), hjust = -0.2, size = 4) +
  scale_color_identity()

# Extract random effects (conditional modes)
re <- ranef(nb_model_slope)$cond$UCSF_MRN

# Inspect structure
head(re)

# Extract random slopes for Telehealth effect
re_slopes <- re %>%
  tibble::rownames_to_column("UCSF_MRN") %>%
  select(UCSF_MRN, contains("Telehealth"))

# Rename for clarity
colnames(re_slopes)[2] <- "Telehealth_random_slope"

ggplot(re_slopes, aes(x = Telehealth_random_slope)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(
    title = "Distribution of Patient-specific Random Slopes for Telehealth",
    x = "Deviation from Average Telehealth Effect (log-IRR scale)",
    y = "Number of Patients"
  )

# Get fixed effect for Telehealth from the model
fixef_telehealth <- fixef(nb_model_slope)$cond["TelehealthTelehealth Video"]

print(fixef_telehealth)

# Calculate conditional log-IRR and IRR
re_slopes <- re_slopes %>%
  mutate(
    conditional_logIRR = fixef_telehealth + Telehealth_random_slope,
    conditional_IRR = exp(conditional_logIRR)
  )

ggplot(re_slopes, aes(x = conditional_IRR)) +
  geom_histogram(bins = 30, fill = "darkorange", color = "black") +
  geom_vline(xintercept = exp(fixef_telehealth), linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(
    title = "Conditional Telehealth IRRs by Patient",
    x = "Conditional IRR (Telehealth vs In-person)",
    y = "Number of Patients"
  )

# So random slopes seemed to help, but what about time auto-correlation?

main <- main %>%
  group_by(UCSF_MRN) %>%
  arrange(APPT_DATE) %>%
  mutate(APPT_DATE_NUM = as.numeric(APPT_DATE - min(APPT_DATE)))

# Fit model with AR(1)
nb_model_ar1 <- glmmTMB(
  med_changes_count ~ Telehealth + SEX_ASSIGNED_AT_BIRTH + AGE_AT_ENCOUNTER_s +
    RACE_ETHNICITY + meds_count_s + egfr_s + 
    (1 + Telehealth | UCSF_MRN),
  data = main,
  family = nbinom2,
  dispformula = ~ ar1(APPT_DATE_NUM + 0 | UCSF_MRN)
)

# Using a numeric variable for time doesn't seem to work...???

# So we can use what visit number it was

main <- main %>%
  group_by(UCSF_MRN) %>%
  mutate(visit_number = as.factor(row_number())) %>%
  ungroup()

# Now fit the model with AR(1) correlation
nb_model_ar1 <- glmmTMB(
  med_changes_count ~ Telehealth + SEX_ASSIGNED_AT_BIRTH + AGE_AT_ENCOUNTER_s +
    RACE_ETHNICITY + meds_count_s + egfr_s + 
    (1 + Telehealth | UCSF_MRN),
  data = main,
  family = nbinom2,
  dispformula = ~ ar1(visit_number | UCSF_MRN)  # Use the visit_number factor for time
)

# Patient-level analysis --------------------
patients <- unique(main$UCSF_MRN)

main2017 <- main %>% filter(year(APPT_DATE) == 2021)

patient_encounter_types <- main %>%
  group_by(UCSF_MRN) %>%
  summarize(encounter_types = paste(sort(unique(Telehealth)), collapse = "_")) %>%
  mutate(group = case_when(
    encounter_types == "In-person" ~ "In-person only",
    encounter_types == "Telehealth" ~ "Telehealth only",
    encounter_types == "In-person_Telehealth" ~ "Both",
    TRUE ~ "Other"
  ))

main_grouped <- main %>%
  left_join(patient_encounter_types, by = "UCSF_MRN")

# but looking at things, we will probably have only information on the first encounter
# i.e. age at first encounter, DM at first encounter, MyChart at first encounter, which can change
# over time...

summary_stats <- main_grouped %>%
  group_by(group) %>%
  summarize(
    n_patients = n_distinct(UCSF_MRN),
    avg_age = mean(AGE_AT_ENCOUNTER, na.rm = TRUE),
    sd_age = sd(AGE_AT_ENCOUNTER, na.rm = TRUE),
  )

t1vars <- c("Telehealth", "AGE_AT_ENCOUNTER", "SEX_ASSIGNED_AT_BIRTH", "RACE_ETHNICITY",
            "BMI_WITHIN_6_MO", "MYC_ENROLLMENT_STATUS", "smoking_status", "PMHX_HTN",
            "PMHX_DM", "PMHX_CAD", "PMHX_CVA", "PMHX_CHF", "CPT_CHARGE_CODE", "INSURANCE_GROUP", 
            "egfr", "meds_count", "orders_count", "med_changes_count")

patient_summary <- main %>%
  arrange(UCSF_MRN, APPT_DATE) %>%  # Ensure ordered by date
  group_by(UCSF_MRN) %>%
  slice(1) %>%  # First encounter
  ungroup()

# Now summarize by Telehealth at first encounter
patient_stats_by_telehealth <- patient_summary %>%
  group_by(Telehealth) %>%
  summarise(
    n_patients = n(),
    avg_AGE_AT_ENCOUNTER = mean(AGE_AT_ENCOUNTER, na.rm = TRUE),
    sd_AGE_AT_ENCOUNTER = sd(AGE_AT_ENCOUNTER, na.rm = TRUE),
    SEX_counts = list(table(SEX_ASSIGNED_AT_BIRTH)),
    RACE_counts = list(table(RACE_ETHNICITY)),
    avg_BMI = mean(BMI_WITHIN_6_MO, na.rm = TRUE),
    sd_BMI = sd(BMI_WITHIN_6_MO, na.rm = TRUE),
    avg_eGFR = mean(egfr, na.rm = TRUE),
    sd_eGFR = sd(egfr, na.rm = TRUE)
  )

# View
patient_stats_by_telehealth

## MISC Can delete and it won't affect anything ----------------
patient_summary <- main %>%
  group_by(UCSF_MRN) %>%
  summarise(
    avg_age = mean(AGE_AT_ENCOUNTER, na.rm = TRUE),
    avg_egfr = mean(egfr, na.rm = TRUE),
    sex = SEX_ASSIGNED_AT_BIRTH[which.max(tabulate(match(SEX_ASSIGNED_AT_BIRTH, unique(SEX_ASSIGNED_AT_BIRTH))))],
    race = RACE_ETHNICITY2[which.max(tabulate(match(RACE_ETHNICITY2, unique(RACE_ETHNICITY2))))],
    insurance = INSURANCE_GROUP[which.max(tabulate(match(INSURANCE_GROUP, unique(INSURANCE_GROUP))))],
    .groups = "drop"
  )

overall_means <- patient_summary %>%
  summarise(
    mean_age = mean(avg_age, na.rm = TRUE),
    sd_age = sd(avg_age, na.rm = TRUE),
    mean_egfr = mean(avg_egfr, na.rm = TRUE),
    sd_egfr = sd(avg_egfr, na.rm = TRUE)
  )

sex_dist <- patient_summary %>%
  count(sex) %>%
  mutate(percent = round(n / sum(n) * 100, 1))

race_dist <- patient_summary %>%
  count(race) %>%
  mutate(percent = round(n / sum(n) * 100, 1))

insurance_dist <- patient_summary %>%
  count(insurance) %>%
  mutate(percent = round(n / sum(n) * 100, 1))

population_overview <- list(
  "Overall Averages" = overall_means,
  "Sex Distribution (%)" = sex_dist,
  "Race Distribution (%)" = race_dist,
  "Insurance Distribution (%)" = insurance_dist
)

population_overview

main2020_2021 <- main2020_2021 %>% left_join(notes, by = "PAT_ENC_CSN_ID")

patient_ids <- main %>%
  count(UCSF_MRN) %>%
  filter(n >= 3) %>%
  pull(UCSF_MRN)

# Loop through a few patients
for (id in head(patient_ids, 5)) {
  idx <- which(main$UCSF_MRN == id)
  
  # Extract residuals for this patient
  sub_res <- recalculateResiduals(sim_res, select = idx)
  
  # Sort dates to ensure proper order
  sub_dates <- main$APPT_DATE[idx]
  if (anyDuplicated(sub_dates) == 0) {  # Ensure unique dates
    cat("Patient:", id, "\n")
    print(testTemporalAutocorrelation(sub_res, time = sub_dates))
  } else {
    cat("Patient:", id, "has duplicate dates, skipping.\n")
  }
}

nb_model_hierarchical2 <- glmmTMB(
  med_changes_count ~ Telehealth + SEX_ASSIGNED_AT_BIRTH + 
    AGE_AT_ENCOUNTER_s + RACE_ETHNICITY2 + New_Old + PMHX_HTN + PMHX_DM +
    Time_since_last_visit_s + meds_count_s + egfr_s +
    (1 + Telehealth | UCSF_MRN) +
    (1 | ENCOUNTER_NOTE_AUTHOR),
  data = main,
  family = nbinom2
)

tidy_nb_slope <- broom.mixed::tidy(nb_model_hierarchical2, effects = "fixed", conf.int = TRUE, exponentiate = TRUE)

irr_table_slope <- broom.mixed::tidy(nb_model_hierarchical2) %>%
  filter(effect == "fixed") %>%
  mutate(
    IRR = exp(estimate),
    lower = exp(estimate - 1.96 * std.error),
    upper = exp(estimate + 1.96 * std.error)
  )

irr_table_slope$significance <- case_when(
  irr_table_slope$p.value < 0.001 ~ "***",
  irr_table_slope$p.value < 0.01  ~ "**",
  irr_table_slope$p.value < 0.05  ~ "*",
  TRUE ~ ""
)

# Print
# print(irr_table, digits = 3)
# 
# # Add a column for color based on significance
# irr_table$color <- ifelse(irr_table$p.value < 0.05, "red", "grey")

# Add a column for color based on significance
irr_table_slope$color <- ifelse(irr_table_slope$p.value < 0.05, "red", "grey")

print(irr_table_slope, digits = 3)

# Plot the IRR with significance stars and color highlighting
# ggplot(irr_table, aes(x = reorder(term, IRR), y = IRR, color = color)) +
#   geom_point(size = 3) +
#   geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
#   coord_flip() +
#   theme_minimal() +
#   labs(
#     x = "Predictor",
#     y = "Incidence Rate Ratio (IRR)",
#     title = "IRR (with 95% CI) for Medication Changes"
#   ) +
#   geom_text(aes(label = significance), hjust = -0.2, size = 4) +
#   scale_color_identity()  # Use the color column directly

irr_table_slope$term

# new term names to be prettier

irr_table_slope$term <- c("Intercept", "Telehealth", "Male Sex", "Age", 
                          "Asian Race", "Black or African American Race",
                          "Latinx Race", "Unknown/Other Race", "First Visit",
                          "Hypertension", "Diabetes",
                          "Time Since Last Visit", "Existing Medication Count",
                          "Estimated Glomerular Filtration Rate")

irr_table_slope <- irr_table_slope %>%
  filter(term != "Intercept")

print(irr_table_slope, digits = 3)

# Plot the IRR with significance stars and color highlighting
ggplot(irr_table_slope, aes(x = reorder(term, IRR), y = IRR, color = color)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_minimal() +
  labs(
    x = "Predictor",
    y = "Incidence Rate Ratio",
    title = "Incidence Rate Ratio for Medication Changes"
  ) +
  #geom_text(aes(label = significance), hjust = -0.2, size = 4) +
  scale_color_identity()  # Use the color column directly

sim_res <- simulateResiduals(nb_model_hierarchical2, n = 1000)
model_data <- nb_model_hierarchical2$frame  # extract cleaned data
plotResiduals(sim_res, form = model_data$AGE_AT_ENCOUNTER_s)
plotResiduals(sim_res, form = model_data$meds_count_s)
plotResiduals(sim_res, form = model_data$Time_since_last_visit_s)
plotResiduals(sim_res, form = model_data$egfr_s)

plotResiduals(sim_res, form = predict(nb_model_hierarchical2, type = "response"))
plotResiduals(sim_res, form = predict(nb_model_hierarchical2, type = "link"))

testDispersion(sim_res)
testZeroInflation(sim_res)
testUniformity(sim_res)
testOutliers(sim_res)


rm(list=ls())
d = read.csv("/Users/keerthim/Documents/ST4060/ST6090/CLUSTERING DATA/BSH_compound_protein_bile_acid.csv",header=TRUE)

head(d)
dim(d)
str(d)

d$Replicate= as.factor(d$Replicate)
d$Run=as.factor(d$Run)
d$Inhibitor=as.factor(d$Inhibitor)
d$BSH_Enzyme=as.factor(d$BSH_Enzyme)


head(d)
dim(d)
str(d)

numeric_cols <- sapply(d, is.numeric)
d_numeric <- d[, numeric_cols]

pairs(d_numeric[, 1:5], pch = 20)



##################### PERCENTAGE RELATIVE CONCENTRATION TO ANALYSE THE Enzyme  ACTION#############
library(dplyr)
library(tidyr)


bile_acids =  c("TUDCA", "TCDCA", "TDCA", "TCA", "TLCA",
                "GUDCA", "GCDCA", "GDCA", "GCA", "GLCA")

#Convert from wide to long format
df_long = d %>%
  pivot_longer(cols = all_of(bile_acids),
               names_to = "Bile_Acid",
               values_to = "Concentration")

#Extract CT+NI values to use as reference
ct_ni = df_long %>%
  filter(BSH_Enzyme == "Ct", Inhibitor == "NI", Replicate == "Average") %>%
  select(Run, Bile_Acid, CT_Concentration = Concentration)

#Join CT+NI back to data and compute % change
df_relative = df_long %>%
  filter(Replicate == "Average", BSH_Enzyme != "Ct",Inhibitor == "NI") %>%
  left_join(ct_ni, by = c("Run", "Bile_Acid")) %>%
  mutate(Relative_Percent_Change = ((CT_Concentration-Concentration) / CT_Concentration) * 100)

#Summarise mean relative change
df_summary = df_relative %>%
  group_by(Run, BSH_Enzyme, Bile_Acid) %>%
  summarise(Mean_Relative_Change = round(mean(Relative_Percent_Change, na.rm = TRUE), 2), .groups = "drop")

#View the table
print(df_summary)

#Export the table
write.csv(df_summary, "mean_relative_ba_change.csv", row.names = FALSE)


############make it as 6 enzymes , 10 bile acids############
summary(df_summary)

sd(df_summary$Mean_Relative_Change, na.rm = TRUE)

result = t.test(df_summary$Mean_Relative_Change, conf.level = 0.95)
print(result$conf.int)

###confidence interval computation

# Compute mean and standard error
n = length(df_summary$Mean_Relative_Change)
mean_val = mean(df_summary$Mean_Relative_Change)
std_dev = sd(df_summary$Mean_Relative_Change) 
std_error = std_dev / sqrt(n-1)

# Get t-critical value (for 95% CI)
confidence_level = 0.95
alpha = 1 - confidence_level
t_critical = qt(1 - alpha/2, df = n - 1)

# Margin of error
margin_error = t_critical * std_error

# Confidence interval
ci_lower = mean_val - margin_error
ci_upper = mean_val + margin_error

cat(
  "Mean:", mean_val, "\n",
  "95% CI: (", ci_lower, ",", ci_upper, ")"
)


################plotting###########
library(dplyr)

###########pull the negative values:##########

df_summary[df_summary$Mean_Relative_Change <= 0, ]  

#############filter the row#################
subset(d[, c("Run", "BSH_Enzyme", "GCA", "GUDCA")], 
       (Run == "Run 3" & BSH_Enzyme == "B" ) |
         (Run == "Run 3" & BSH_Enzyme == "Ct" ) |
         (Run == "Run 4" & BSH_Enzyme == "B" ) |
         (Run == "Run 4" & BSH_Enzyme == "Ct" ) |
         (Run == "Run 5" & BSH_Enzyme == "1011c" ) |
         (Run == "Run 5" & BSH_Enzyme == "Ct" ) |
         (Run == "Run 6" & BSH_Enzyme == "B" ) |
         (Run == "Run 6" & BSH_Enzyme == "Ct"))



##############summarize accross the run################
df_plot = df_summary %>%
  group_by(BSH_Enzyme, Bile_Acid) %>%
  summarise(
    Mean_Change = mean(Mean_Relative_Change, na.rm = TRUE),
    .groups = "drop"
  )


library(pheatmap)
library(tidyr)
library(tibble)


df_wide = df_plot %>%
  pivot_wider(names_from = Bile_Acid, values_from = Mean_Change) %>%
  column_to_rownames("BSH_Enzyme")


green_white_red = colorRampPalette(c("red", "white", "darkgreen"))

#heatmap
pheatmap(as.matrix(df_wide),
         color = green_white_red(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE,
         number_format = "%.1f",
         main = "Enzyme Action on Bile Acids (% Change)",
         fontsize_number = 10,
         fontsize = 10,
         number_color = "white",
         legend = TRUE,
         fontface_number = "bold")


##dont show

ggplot(df_plot, aes(x = BSH_Enzyme, y = Mean_Change, fill = BSH_Enzyme)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ Bile_Acid, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Mean Relative % Change per Enzyme on Each Bile Acid",
       x = "BSH Enzyme", y = "% Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

####################enzyme +inhibitor ###############


library(dplyr)
library(tidyr)
library(ggplot2)

bile_acids =  c("TUDCA", "TCDCA", "TDCA", "TCA", "TLCA",
                "GUDCA", "GCDCA", "GDCA", "GCA", "GLCA")


df_long = d %>%
  pivot_longer(cols = all_of(bile_acids),
               names_to = "Bile_Acid",
               values_to = "Concentration")

#Ct+NI concentrations (control)
ct_ni = df_long %>%
  filter(BSH_Enzyme == "Ct", Inhibitor == "NI", Replicate == "Average") %>%
  select(Run, Bile_Acid, CT_NI_Concentration = Concentration)

# Enzyme + NI
enzyme_ni = df_long %>%
  filter(BSH_Enzyme != "Ct", Inhibitor == "NI", Replicate == "Average") %>%
  left_join(ct_ni, by = c("Run", "Bile_Acid")) %>%
  mutate(Relative_Percent_Change = ((CT_NI_Concentration-Concentration ) / CT_NI_Concentration) * 100,
         Condition = "Enzyme + NI")

# Enzyme + Inhibitor
enzyme_inhibitor = df_long %>%
  filter(BSH_Enzyme != "Ct", Inhibitor != "NI", Replicate == "Average") %>%
  left_join(ct_ni, by = c("Run", "Bile_Acid")) %>%
  mutate(Relative_Percent_Change = ((CT_NI_Concentration-Concentration ) / CT_NI_Concentration) * 100,
         Condition = "Enzyme + Inhibitor")

# Combine data
df_combined = bind_rows(
  enzyme_ni %>% select(Run, BSH_Enzyme, Bile_Acid, Inhibitor, Condition, Relative_Percent_Change),
  enzyme_inhibitor %>% select(Run, BSH_Enzyme, Bile_Acid, Inhibitor, Condition, Relative_Percent_Change)
)

# Filter out extreme changes (less than zero and greater than 100)
df_extreme = df_combined %>%
  filter(Relative_Percent_Change > 100 | Relative_Percent_Change <=0)

############find in the main data frame##################

# Keep only values between 0 and 100
df_combined_filtered = df_combined %>%
  filter(Relative_Percent_Change >= 0 & Relative_Percent_Change <= 100)

# summary with filtered data
df_summary = df_combined_filtered %>%
  group_by(Run, BSH_Enzyme, Inhibitor, Bile_Acid, Condition) %>%
  summarise(Mean_Relative_Change = round(mean(Relative_Percent_Change, na.rm = TRUE), 2), .groups = "drop")

#extreme values for inspection
write.csv(df_extreme, "excluded_extreme_relative_changes.csv", row.names = FALSE)

#iltered summary
write.csv(df_summary, "mean_relative_ba_change_inhibitor_summary_filtered.csv", row.names = FALSE)


######################6 separate bar plots (one for each enzyme and its inhibitors)#######################
df_JL885 = df_summary %>% filter(BSH_Enzyme == "JL885")
df_A = df_summary %>% filter(BSH_Enzyme == "A")
df_B = df_summary %>% filter(BSH_Enzyme == "B")
df_F = df_summary %>% filter(BSH_Enzyme == "F")
df_1011c = df_summary %>% filter(BSH_Enzyme == "1011c")
df_T2 = df_summary %>% filter(BSH_Enzyme == "T2")



###########enzyme JL885########
#######plotting########
df_JL885$NI_Status <- ifelse(df_JL885$Inhibitor == "NI", "NI", "Inhibitor")

library(dplyr)

# Step 1: Calculate baseline (enzyme+ni) means for each bile acid
baseline_means <- df_JL885 %>%
  filter(Inhibitor == "NI") %>%
  group_by(Bile_Acid) %>%
  summarise(
    Baseline_Mean = mean(Mean_Relative_Change, na.rm = TRUE)
  )

# Step 2: Join baseline data and calculate percentage decrease
df_JL885 <- df_JL885 %>%
  left_join(baseline_means, by = "Bile_Acid") %>%
  mutate(
    Percent_Decrease = ifelse(
      Inhibitor != "NI",
      (Baseline_Mean - Mean_Relative_Change) / Baseline_Mean * 100,
      NA
    )
  ) %>%
  mutate(
    Inhibition_Level = case_when(
      Inhibitor == "NI" ~ "Control (NI)",
      Percent_Decrease < 10 ~ "No Inhibition",
      Percent_Decrease >= 10 & Percent_Decrease < 25 ~ "Weak Inhibition",
      Percent_Decrease >= 25 & Percent_Decrease < 50 ~ "Moderate Inhibition",
      Percent_Decrease >= 50 ~ "Strong Inhibition",
      TRUE ~ "Undefined"
    )
  )



colors <- c(
  "Control (NI)"      = "#1F78B4",  
  "No Inhibition"     = "#1F78B4", 
  "Weak Inhibition"   = "#6BAED6",  
  "Moderate Inhibition" = "#F08080",
  "Strong Inhibition" = "#E31A1C"   
)



ggplot(df_JL885, aes(x = Bile_Acid, y = Mean_Relative_Change, fill = Inhibition_Level)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = colors, name = "Inhibition Level") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~Inhibitor, ncol = 5) +
  labs(title = "JL885 – % Change in Bile Acids vs Ct+NI",
       y = "% Change from Control (Ct+NI)",
       x = "Bile Acid") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")



#########enzyme A##############################

ni_data=df_A%>%  
  filter(Inhibitor == "NI") %>%
  select(Run, BSH_Enzyme, Bile_Acid, Mean_Relative_Change)

CI_summary <- ni_data %>%
  group_by(BSH_Enzyme, Bile_Acid) %>%
  summarise(
    N = n(),
    Mean = mean(Mean_Relative_Change, na.rm = TRUE),
    SD = sd(Mean_Relative_Change, na.rm = TRUE),
    SE = SD / sqrt(N),
    CI_Lower = Mean - qnorm(0.975) * SE,  # 95% CI
    CI_Upper = Mean + qnorm(0.975) * SE,
    .groups = "drop"
  )
CI_summary



#######plotting########
df_A$NI_Status <- ifelse(df_A$Inhibitor == "NI", "NI", "Inhibitor")

#Calculate baseline (enzyme+ni) means for each bile acid
baseline_means <- df_A %>%
  filter(Inhibitor == "NI") %>%
  group_by(Bile_Acid) %>%
  summarise(
    Baseline_Mean = mean(Mean_Relative_Change, na.rm = TRUE)
  )

#Join baseline data and calculate percentage decrease
df_A <- df_A %>%
  left_join(baseline_means, by = "Bile_Acid") %>%
  mutate(
    Percent_Decrease = ifelse(
      Inhibitor != "NI",
      (Baseline_Mean - Mean_Relative_Change) / Baseline_Mean * 100,
      NA
    )
  ) %>%
  mutate(
    Inhibition_Level = case_when(
      Inhibitor == "NI" ~ "Control (NI)",
      Percent_Decrease < 10 ~ "No Inhibition",
      Percent_Decrease >= 10 & Percent_Decrease < 25 ~ "Weak Inhibition",
      Percent_Decrease >= 25 & Percent_Decrease < 50 ~ "Moderate Inhibition",
      Percent_Decrease >= 50 ~ "Strong Inhibition",
      TRUE ~ "Undefined"
    )
  )


colors <- c(
  "Control (NI)"      = "#1F78B4",  
  "No Inhibition"     = "#1F78B4", 
  "Weak Inhibition"   = "#6BAED6",  
  "Moderate Inhibition" = "#F08080",
  "Strong Inhibition" = "#E31A1C"   
)



ggplot(df_A, aes(x = Bile_Acid, y = Mean_Relative_Change, fill = Inhibition_Level)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = colors, name = "Inhibition Level") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~Inhibitor, ncol = 5) +
  labs(title = "Enzyme A – % Change in Bile Acids vs Ct+NI",
       y = "% Change from Control (Ct+NI)",
       x = "Bile Acid") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

##########few discrepant columns for inspection#######
df_A[df_A$Bile_Acid=="GCA",]
d[d$Run=="Run 2"&d$BSH_Enzyme=="A",][,c("Run","BSH_Enzyme","Inhibitor","GCA")]
d[d$Run=="Run 2"&d$BSH_Enzyme=="Ct",][,c("Run","BSH_Enzyme","Inhibitor","GCA")]




#######################Enzyme B#################################
ni_data=df_B%>%  
  filter(Inhibitor == "NI") %>%
  select(Run, BSH_Enzyme, Bile_Acid, Mean_Relative_Change)

CI_summary <- ni_data %>%
  group_by(BSH_Enzyme, Bile_Acid) %>%
  summarise(
    N = n(),
    Mean = mean(Mean_Relative_Change, na.rm = TRUE),
    SD = sd(Mean_Relative_Change, na.rm = TRUE),
    SE = SD / sqrt(N),
    CI_Lower = Mean - qnorm(0.975) * SE,  # 95% CI
    CI_Upper = Mean + qnorm(0.975) * SE,
    .groups = "drop"
  )
CI_summary

#######PLOTTING########
df_B$NI_Status <- ifelse(df_B$Inhibitor == "NI", "NI", "Inhibitor")

#Calculate baseline (enzyme+ni) means for each bile acid
baseline_means <- df_B %>%
  filter(Inhibitor == "NI") %>%
  group_by(Bile_Acid) %>%
  summarise(
    Baseline_Mean = mean(Mean_Relative_Change, na.rm = TRUE)
  )

#Join baseline data and calculate percentage decrease
df_B <- df_B %>%
  left_join(baseline_means, by = "Bile_Acid") %>%
  mutate(
    Percent_Decrease = ifelse(
      Inhibitor != "NI",
      (Baseline_Mean - Mean_Relative_Change) / Baseline_Mean * 100,
      NA
    )
  ) %>%
  mutate(
    Inhibition_Level = case_when(
      Inhibitor == "NI" ~ "Control (NI)",
      Percent_Decrease < 10 ~ "No Inhibition",
      Percent_Decrease >= 10 & Percent_Decrease < 25 ~ "Weak Inhibition",
      Percent_Decrease >= 25 & Percent_Decrease < 50 ~ "Moderate Inhibition",
      Percent_Decrease >= 50 ~ "Strong Inhibition",
      TRUE ~ "Undefined"
    )
  )


colors <- c(
  "Control (NI)"      = "#1F78B4",  
  "No Inhibition"     = "#1F78B4", 
  "Weak Inhibition"   = "#6BAED6",  
  "Moderate Inhibition" = "#F08080",
  "Strong Inhibition" = "#E31A1C"   
)


ggplot(df_B, aes(x = Bile_Acid, y = Mean_Relative_Change, fill = Inhibition_Level)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = colors, name = "Inhibition Level") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~Inhibitor, ncol = 5) +
  labs(title = "Enzyme B – % Change in Bile Acids vs Ct+NI",
       y = "% Change from Control (Ct+NI)",
       x = "Bile Acid") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")



##########inspecting few rows#########
df_B[df_B$Bile_Acid=="GCA" & df_B$Inhibitor=="XI" ,]
df_B[df_B$Bile_Acid=="GCA" & df_B$Inhibitor=="NI" ,]

df_extreme


#######################enzyme F#########################################

ni_data=df_F%>%  
  filter(Inhibitor == "NI") %>%
  select(Run, BSH_Enzyme, Bile_Acid, Mean_Relative_Change)

CI_summary <- ni_data %>%
  group_by(BSH_Enzyme, Bile_Acid) %>%
  summarise(
    N = n(),
    Mean = mean(Mean_Relative_Change, na.rm = TRUE),
    SD = sd(Mean_Relative_Change, na.rm = TRUE),
    SE = SD / sqrt(N),
    CI_Lower = Mean - qnorm(0.975) * SE,  # 95% CI
    CI_Upper = Mean + qnorm(0.975) * SE,
    .groups = "drop"
  )
CI_summary

######inspecting few records
df_F[df_F$Bile_Acid=="TCA" & df_F$Inhibitor=="NI",]

d[d$Run=="Run 2" & d$BSH_Enzyme=="F" & d$Inhibitor=="NI",][,c("Run","BSH_Enzyme","Inhibitor","TCA")]
d[d$Run=="Run 2" & d$BSH_Enzyme=="Ct" & d$Inhibitor=="NI",][,c("Run","BSH_Enzyme","Inhibitor","TCA")]


d[d$Run=="Run 3" & d$BSH_Enzyme=="F" & d$Inhibitor=="NI",][,c("Run","BSH_Enzyme","Inhibitor","TCA")]
d[d$Run=="Run 3" & d$BSH_Enzyme=="Ct" & d$Inhibitor=="NI",][,c("Run","BSH_Enzyme","Inhibitor","TCA")]

#######PLOTTING########
df_F$NI_Status <- ifelse(df_F$Inhibitor == "NI", "NI", "Inhibitor")

#Calculate baseline (enzyme+ni) means for each bile acid
baseline_means <- df_F %>%
  filter(Inhibitor == "NI") %>%
  group_by(Bile_Acid) %>%
  summarise(
    Baseline_Mean = mean(Mean_Relative_Change, na.rm = TRUE)
  )

#Join baseline data and calculate percentage decrease
df_F <- df_F %>%
  left_join(baseline_means, by = "Bile_Acid") %>%
  mutate(
    Percent_Decrease = ifelse(
      Inhibitor != "NI",
      (Baseline_Mean - Mean_Relative_Change) / Baseline_Mean * 100,
      NA
    )
  ) %>%
  mutate(
    Inhibition_Level = case_when(
      Inhibitor == "NI" ~ "Control (NI)",
      Percent_Decrease < 10 ~ "No Inhibition",
      Percent_Decrease >= 10 & Percent_Decrease < 25 ~ "Weak Inhibition",
      Percent_Decrease >= 25 & Percent_Decrease < 50 ~ "Moderate Inhibition",
      Percent_Decrease >= 50 ~ "Strong Inhibition",
      TRUE ~ "Undefined"
    )
  )


colors <- c(
  "Control (NI)"      = "#1F78B4",  
  "No Inhibition"     = "#1F78B4", 
  "Weak Inhibition"   = "#6BAED6",  
  "Moderate Inhibition" = "#F08080",
  "Strong Inhibition" = "#E31A1C"   
)


ggplot(df_F, aes(x = Bile_Acid, y = Mean_Relative_Change, fill = Inhibition_Level)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = colors, name = "Inhibition Level") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~Inhibitor, ncol = 5) +
  labs(title = "Enzyme F – % Change in Bile Acids vs Ct+NI",
       y = "% Change from Control (Ct+NI)",
       x = "Bile Acid") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")



##########inspecting few rows#########
df_F[df_F$Inhibitor=="R",]

df_extreme

##################enzyme 1011c###############################################
ni_data=df_1011c%>%  
  filter(Inhibitor == "NI") %>%
  select(Run, BSH_Enzyme, Bile_Acid, Mean_Relative_Change)

CI_summary <- ni_data %>%
  group_by(BSH_Enzyme, Bile_Acid) %>%
  summarise(
    N = n(),
    Mean = mean(Mean_Relative_Change, na.rm = TRUE),
    SD = sd(Mean_Relative_Change, na.rm = TRUE),
    SE = SD / sqrt(N),
    CI_Lower = Mean - qnorm(0.975) * SE,  # 95% CI
    CI_Upper = Mean + qnorm(0.975) * SE,
    .groups = "drop"
  )
CI_summary

#######PLOTTING########

df_1011c$NI_Status <- ifelse(df_1011c$Inhibitor == "NI", "NI", "Inhibitor")

#Calculate baseline (enzyme+ni) means for each bile acid
baseline_means <- df_1011c %>%
  filter(Inhibitor == "NI") %>%
  group_by(Bile_Acid) %>%
  summarise(
    Baseline_Mean = mean(Mean_Relative_Change, na.rm = TRUE)
  )

#Join baseline data and calculate percentage decrease
df_1011c <- df_1011c %>%
  left_join(baseline_means, by = "Bile_Acid") %>%
  mutate(
    Percent_Decrease = ifelse(
      Inhibitor != "NI",
      (Baseline_Mean - Mean_Relative_Change) / Baseline_Mean * 100,
      NA
    )
  ) %>%
  mutate(
    Inhibition_Level = case_when(
      Inhibitor == "NI" ~ "Control (NI)",
      Percent_Decrease < 10 ~ "No Inhibition",
      Percent_Decrease >= 10 & Percent_Decrease < 25 ~ "Weak Inhibition",
      Percent_Decrease >= 25 & Percent_Decrease < 50 ~ "Moderate Inhibition",
      Percent_Decrease >= 50 ~ "Strong Inhibition",
      TRUE ~ "Undefined"
    )
  )


colors <- c(
  "Control (NI)"      = "#1F78B4",  
  "No Inhibition"     = "#1F78B4", 
  "Weak Inhibition"   = "#6BAED6",  
  "Moderate Inhibition" = "#F08080",
  "Strong Inhibition" = "#E31A1C"   
)


ggplot(df_1011c, aes(x = Bile_Acid, y = Mean_Relative_Change, fill = Inhibition_Level)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = colors, name = "Inhibition Level") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~Inhibitor, ncol = 5) +
  labs(title = "Enzyme 1011c – % Change in Bile Acids vs Ct+NI",
       y = "% Change from Control (Ct+NI)",
       x = "Bile Acid") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")



##########inspecting few rows#########
df_1011c[df_1011c$Inhibitor=="R" ,]

df_extreme

###########################ENZYME T2###########################################
ni_data=df_T2%>%  
  filter(Inhibitor == "NI") %>%
  select(Run, BSH_Enzyme, Bile_Acid, Mean_Relative_Change)

CI_summary <- ni_data %>%
  group_by(BSH_Enzyme, Bile_Acid) %>%
  summarise(
    N = n(),
    Mean = mean(Mean_Relative_Change, na.rm = TRUE),
    SD = sd(Mean_Relative_Change, na.rm = TRUE),
    SE = SD / sqrt(N),
    CI_Lower = Mean - qnorm(0.975) * SE,  # 95% CI
    CI_Upper = Mean + qnorm(0.975) * SE,
    .groups = "drop"
  )
CI_summary

##inspecing
df_T2[df_T2$Bile_Acid=="GUDCA" & df_T2$Run=="Run 6",]



#######PLOTTING########
df_T2$NI_Status <- ifelse(df_T2$Inhibitor == "NI", "NI", "Inhibitor")

#Calculate baseline (enzyme+ni) means for each bile acid
baseline_means <- df_T2 %>%
  filter(Inhibitor == "NI") %>%
  group_by(Bile_Acid) %>%
  summarise(
    Baseline_Mean = mean(Mean_Relative_Change, na.rm = TRUE)
  )

#Join baseline data and calculate percentage decrease
df_T2 <- df_T2 %>%
  left_join(baseline_means, by = "Bile_Acid") %>%
  mutate(
    Percent_Decrease = ifelse(
      Inhibitor != "NI",
      (Baseline_Mean - Mean_Relative_Change) / Baseline_Mean * 100,
      NA
    )
  ) %>%
  mutate(
    Inhibition_Level = case_when(
      Inhibitor == "NI" ~ "Control (NI)",
      Percent_Decrease < 10 ~ "No Inhibition",
      Percent_Decrease >= 10 & Percent_Decrease < 25 ~ "Weak Inhibition",
      Percent_Decrease >= 25 & Percent_Decrease < 50 ~ "Moderate Inhibition",
      Percent_Decrease >= 50 ~ "Strong Inhibition",
      TRUE ~ "Undefined"
    )
  )


colors <- c(
  "Control (NI)"      = "#1F78B4",  
  "No Inhibition"     = "#1F78B4", 
  "Weak Inhibition"   = "#6BAED6",  
  "Moderate Inhibition" = "#F08080",
  "Strong Inhibition" = "#E31A1C"   
)


ggplot(df_T2, aes(x = Bile_Acid, y = Mean_Relative_Change, fill = Inhibition_Level)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = colors, name = "Inhibition Level") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~Inhibitor, ncol = 5) +
  labs(title = "Enzyme T2 – % Change in Bile Acids vs Ct+NI",
       y = "% Change from Control (Ct+NI)",
       x = "Bile Acid") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")


df_extreme
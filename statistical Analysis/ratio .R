rm(list=ls())
d = read.csv("/Users/keerthim/Documents/ST4060/ST6090/CLUSTERING DATA/BSH_compound_protein_bile_acid_detailed.csv",header=TRUE)

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



##################### Log ratio To Analyse THE Enzyme  ACTION#############
library(dplyr)
library(tidyr)

bile_acids <- c("TUDCA", "TCDCA", "TDCA", "TCA", "TLCA",
                "GUDCA", "GCDCA", "GDCA", "GCA", "GLCA")

df_long <- d %>%
  pivot_longer(cols = all_of(bile_acids),
               names_to = "Bile_Acid",
               values_to = "Concentration")

# Extract control (Ct + NI)
ct_ni <- df_long %>%
  filter(BSH_Enzyme == "Ct", Inhibitor == "NI", Replicate == "Average") %>%
  select(Run, Bile_Acid, CT_Concentration = Concentration)

#compute log2 ratio
df_logratio <- df_long %>%
  filter(Replicate == "Average", BSH_Enzyme != "Ct", Inhibitor == "NI") %>%
  left_join(ct_ni, by = c("Run", "Bile_Acid")) %>%
  mutate(
    Log2_Ratio = log2(CT_Concentration/Concentration)
  )

# Summarise log2 ratios
df_summary_log <- df_logratio %>%
  group_by(Run, BSH_Enzyme, Bile_Acid) %>%
  summarise(
    Mean_Log2_Ratio = round(mean(Log2_Ratio, na.rm = TRUE), 4),
    .groups = "drop"
  )


print(df_summary_log)

############make it as 6 enzymes , 10 bile acids############
summary(df_summary_log)

sd(df_summary_log$Mean_Log2_Ratio, na.rm = TRUE)

################plotting###########
###########pull the negative values:##########
df_summary_log[df_summary_log$Mean_Log2_Ratio <= 0, ]

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
df_plot = df_summary_log %>%
  group_by(BSH_Enzyme, Bile_Acid) %>%
  summarise(
    Mean_Change_Ratio = mean(Mean_Log2_Ratio, na.rm = TRUE),
    .groups = "drop"
  )


library(pheatmap)
library(tidyr)
library(tibble)


df_wide = df_plot %>%
  pivot_wider(names_from = Bile_Acid, values_from = Mean_Change_Ratio) %>%
  column_to_rownames("BSH_Enzyme")


green_white_red = colorRampPalette(c("red", "lightgreen", "darkgreen"))

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

####################enzyme +inhibitor ###############
library(tidyr)
library(dplyr)

bile_acids =  c("TUDCA", "TCDCA", "TDCA", "TCA", "TLCA",
                "GUDCA", "GCDCA", "GDCA", "GCA", "GLCA")

# Convert to long format
df_long = d %>%
  pivot_longer(cols = all_of(bile_acids),
               names_to = "Bile_Acid",
               values_to = "Concentration")

# Ct+NI concentrations (control)
ct_ni = df_long %>%
  filter(BSH_Enzyme == "Ct", Inhibitor == "NI", Replicate == "Average") %>%
  select(Run, Bile_Acid, CT_NI_Concentration = Concentration)

# Enzyme + NI condition
enzyme_ni = df_long %>%
  filter(BSH_Enzyme != "Ct", Inhibitor == "NI", Replicate == "Average") %>%
  left_join(ct_ni, by = c("Run", "Bile_Acid")) %>%
  mutate(
    Log2_Ratio = log2(CT_NI_Concentration/Concentration),
    Condition = "Enzyme + NI"
  )

# Enzyme + Inhibitor condition
enzyme_inhibitor = df_long %>%
  filter(BSH_Enzyme != "Ct", Inhibitor != "NI", Replicate == "Average") %>%
  left_join(ct_ni, by = c("Run", "Bile_Acid")) %>%
  mutate(
    Log2_Ratio = log2(CT_NI_Concentration/Concentration),
    Condition = "Enzyme + Inhibitor"
  )

# Combine log2 ratio data
df_combined = bind_rows(
  enzyme_ni %>% select(Run, BSH_Enzyme, Bile_Acid, Inhibitor, Condition, Log2_Ratio),
  enzyme_inhibitor %>% select(Run, BSH_Enzyme, Bile_Acid, Inhibitor, Condition, Log2_Ratio)
)

# Summary mean log2 ratio per group
df_summary = df_combined %>%
  group_by(Run, BSH_Enzyme, Inhibitor, Bile_Acid, Condition) %>%
  summarise(Mean_Log2_Ratio = round(mean(Log2_Ratio, na.rm = TRUE), 3), .groups = "drop")


######################6 separate bar plots (one for each enzyme and its inhibitors)#######################
df_JL885 = df_summary %>% filter(BSH_Enzyme == "JL885")
df_A = df_summary %>% filter(BSH_Enzyme == "A")
df_B = df_summary %>% filter(BSH_Enzyme == "B")
df_F = df_summary %>% filter(BSH_Enzyme == "F")
df_1011c = df_summary %>% filter(BSH_Enzyme == "1011c")
df_T2 = df_summary %>% filter(BSH_Enzyme == "T2")


########### enzyme JL885 ########

df_JL885$NI_Status = ifelse(df_JL885$Inhibitor == "NI", "NI", "Inhibitor")

#Baseline (enzyme + NI) log2 means for each bile acid
baseline_means = df_JL885 %>%
  filter(Inhibitor == "NI") %>%
  group_by(Bile_Acid) %>%
  summarise(
    Baseline_Mean = mean(Mean_Log2_Ratio, na.rm = TRUE)
  )

#inhibition effect
df_JL885 = df_JL885 %>%
  left_join(baseline_means, by = "Bile_Acid") %>%
  mutate(
    Log2_Difference = ifelse(
      Inhibitor != "NI",
      Mean_Log2_Ratio-Baseline_Mean ,
      NA
    )
  ) %>%
   mutate(
      Inhibition_Level = case_when(
        Inhibitor == "NI" ~ "Control (NI)",
        Log2_Difference <= 0 ~ "Strong Inhibition",
        Log2_Difference > 0 & Log2_Difference < 0.5 ~ "Moderate Inhibition",
        Log2_Difference >= 0.5 & Log2_Difference < 1.0 ~ "Weak Inhibition",
        Log2_Difference >= 1.0 ~ "No Inhibition",
        TRUE ~ "Undefined"
      )
  )

colors = c(
  "Control (NI)"        = "#1F78B4",  
  "No Inhibition"       = "#1F78B4", 
  "Weak Inhibition"     = "#6BAED6",  
  "Moderate Inhibition" = "#F08080",
  "Strong Inhibition"   = "#E31A1C"   
)



ggplot(df_JL885, aes(x = Bile_Acid, y = Mean_Log2_Ratio, fill = Inhibition_Level)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = colors, name = "Inhibition Level") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~Inhibitor, ncol = 5) +
  labs(title = "JL885 – Log₂ Fold Change in Bile Acids vs Ct+NI",
       y = "log₂(Ct+NI/Treated))",
       x = "Bile Acid") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")


############RELATION TO H-BOND-DONOR AND H-BOND ACCEPTOR ######
library(dplyr)

inhibitor_features <- d %>%
  select(Inhibitor, Inhibitor_H_Bond_Donor, Inhibitor_H_Bond_Acceptor) %>%
  distinct()

df_JL885_joined <- df_JL885 %>%
  left_join(inhibitor_features, by = "Inhibitor")

lm_simple <- lm(Mean_Log2_Ratio ~ Inhibitor_H_Bond_Donor + Inhibitor_H_Bond_Acceptor,
                data = df_JL885_joined)

summary(lm_simple)



ggplot(df_JL885_joined, aes(x = Inhibitor_H_Bond_Donor, y = Mean_Log2_Ratio, color = Inhibitor)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()



ggplot(df_JL885_joined, aes(x = Inhibitor_H_Bond_Acceptor, y = Mean_Log2_Ratio, color = Inhibitor)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()


########### enzyme A ########

df_A$NI_Status = ifelse(df_A$Inhibitor == "NI", "NI", "Inhibitor")

#Baseline (enzyme + NI) log2 means for each bile acid
baseline_means = df_A %>%
  filter(Inhibitor == "NI") %>%
  group_by(Bile_Acid) %>%
  summarise(
    Baseline_Mean = mean(Mean_Log2_Ratio, na.rm = TRUE)
  )

#inhibition effect
df_A = df_A %>%
  left_join(baseline_means, by = "Bile_Acid") %>%
  mutate(
    Log2_Difference = ifelse(
      Inhibitor != "NI",
      Mean_Log2_Ratio-Baseline_Mean ,
      NA
    )
  ) %>%
  mutate(
    Inhibition_Level = case_when(
      Inhibitor == "NI" ~ "Control (NI)",
      Log2_Difference <= 0 ~ "Strong Inhibition",
      Log2_Difference > 0 & Log2_Difference < 0.5 ~ "Moderate Inhibition",
      Log2_Difference >= 0.5 & Log2_Difference < 1.0 ~ "Weak Inhibition",
      Log2_Difference >= 1.0 ~ "No Inhibition",
      TRUE ~ "Undefined"
    )
  )

colors = c(
  "Control (NI)"        = "#1F78B4",  
  "No Inhibition"       = "#1F78B4", 
  "Weak Inhibition"     = "#6BAED6",  
  "Moderate Inhibition" = "#F08080",
  "Strong Inhibition"   = "#E31A1C"   
)



ggplot(df_A, aes(x = Bile_Acid, y = Mean_Log2_Ratio, fill = Inhibition_Level)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = colors, name = "Inhibition Level") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~Inhibitor, ncol = 5) +
  labs(title = "Enzyme A – Log₂ Fold Change in Bile Acids vs Ct+NI",
       y = "log₂(Ct+NI/Treated)",
       x = "Bile Acid") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

############## enzyme B###############################
df_B$NI_Status = ifelse(df_B$Inhibitor == "NI", "NI", "Inhibitor")

#Baseline (enzyme + NI) log2 means for each bile acid
baseline_means = df_B %>%
  filter(Inhibitor == "NI") %>%
  group_by(Bile_Acid) %>%
  summarise(
    Baseline_Mean = mean(Mean_Log2_Ratio, na.rm = TRUE)
  )

#inhibition effect
df_B = df_B %>%
  left_join(baseline_means, by = "Bile_Acid") %>%
  mutate(
    Log2_Difference = ifelse(
      Inhibitor != "NI",
      Mean_Log2_Ratio-Baseline_Mean ,
      NA
    )
  ) %>%
  mutate(
    Inhibition_Level = case_when(
      Inhibitor == "NI" ~ "Control (NI)",
      Log2_Difference <= 0 ~ "Strong Inhibition",
      Log2_Difference > 0 & Log2_Difference < 0.5 ~ "Moderate Inhibition",
      Log2_Difference >= 0.5 & Log2_Difference < 1.0 ~ "Weak Inhibition",
      Log2_Difference >= 1.0 ~ "No Inhibition",
      TRUE ~ "Undefined"
    )
  )

colors = c(
  "Control (NI)"        = "#1F78B4",  
  "No Inhibition"       = "#1F78B4", 
  "Weak Inhibition"     = "#6BAED6",  
  "Moderate Inhibition" = "#F08080",
  "Strong Inhibition"   = "#E31A1C"   
)



ggplot(df_B, aes(x = Bile_Acid, y = Mean_Log2_Ratio, fill = Inhibition_Level)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = colors, name = "Inhibition Level") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~Inhibitor, ncol = 5) +
  labs(title = "Enzyme B – Log₂ Fold Change in Bile Acids vs Ct+NI",
       y = "log₂(Ct+NI/Treated)",
       x = "Bile Acid") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

#######################enzyme F#########################################
df_F$NI_Status = ifelse(df_F$Inhibitor == "NI", "NI", "Inhibitor")

#Baseline (enzyme + NI) log2 means for each bile acid
baseline_means = df_F %>%
  filter(Inhibitor == "NI") %>%
  group_by(Bile_Acid) %>%
  summarise(
    Baseline_Mean = mean(Mean_Log2_Ratio, na.rm = TRUE)
  )

#inhibition effect
df_F = df_F %>%
  left_join(baseline_means, by = "Bile_Acid") %>%
  mutate(
    Log2_Difference = ifelse(
      Inhibitor != "NI",
      Mean_Log2_Ratio-Baseline_Mean ,
      NA
    )
  ) %>%
  mutate(
    Inhibition_Level = case_when(
      Inhibitor == "NI" ~ "Control (NI)",
      Log2_Difference <= 0 ~ "Strong Inhibition",
      Log2_Difference > 0 & Log2_Difference < 0.5 ~ "Moderate Inhibition",
      Log2_Difference >= 0.5 & Log2_Difference < 1.0 ~ "Weak Inhibition",
      Log2_Difference >= 1.0 ~ "No Inhibition",
      TRUE ~ "Undefined"
    )
  )

colors = c(
  "Control (NI)"        = "#1F78B4",  
  "No Inhibition"       = "#1F78B4", 
  "Weak Inhibition"     = "#6BAED6",  
  "Moderate Inhibition" = "#F08080",
  "Strong Inhibition"   = "#E31A1C"   
)



ggplot(df_F, aes(x = Bile_Acid, y = Mean_Log2_Ratio, fill = Inhibition_Level)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = colors, name = "Inhibition Level") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~Inhibitor, ncol = 5) +
  labs(title = "Enzyme F – Log₂ Fold Change in Bile Acids vs Ct+NI",
       y = "log₂(Ct+NI/Treated)",
       x = "Bile Acid") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

##################enzyme 1011c###############################################
df_1011c$NI_Status = ifelse(df_1011c$Inhibitor == "NI", "NI", "Inhibitor")

#Baseline (enzyme + NI) log2 means for each bile acid
baseline_means = df_1011c %>%
  filter(Inhibitor == "NI") %>%
  group_by(Bile_Acid) %>%
  summarise(
    Baseline_Mean = mean(Mean_Log2_Ratio, na.rm = TRUE)
  )

#inhibition effect
df_1011c = df_1011c %>%
  left_join(baseline_means, by = "Bile_Acid") %>%
  mutate(
    Log2_Difference = ifelse(
      Inhibitor != "NI",
      Mean_Log2_Ratio-Baseline_Mean ,
      NA
    )
  ) %>%
  mutate(
    Inhibition_Level = case_when(
      Inhibitor == "NI" ~ "Control (NI)",
      Log2_Difference <= 0 ~ "Strong Inhibition",
      Log2_Difference > 0 & Log2_Difference < 0.5 ~ "Moderate Inhibition",
      Log2_Difference >= 0.5 & Log2_Difference < 1.0 ~ "Weak Inhibition",
      Log2_Difference >= 1.0 ~ "No Inhibition",
      TRUE ~ "Undefined"
    )
  )

colors = c(
  "Control (NI)"        = "#1F78B4",  
  "No Inhibition"       = "#1F78B4", 
  "Weak Inhibition"     = "#6BAED6",  
  "Moderate Inhibition" = "#F08080",
  "Strong Inhibition"   = "#E31A1C"   
)



ggplot(df_1011c, aes(x = Bile_Acid, y = Mean_Log2_Ratio, fill = Inhibition_Level)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = colors, name = "Inhibition Level") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~Inhibitor, ncol = 5) +
  labs(title = "Enzyme 1011c – Log₂ Fold Change in Bile Acids vs Ct+NI",
       y = "log₂(Ct+NI/Treated)",
       x = "Bile Acid") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

###########################ENZYME T2###########################################
df_T2$NI_Status = ifelse(df_T2$Inhibitor == "NI", "NI", "Inhibitor")

#Baseline (enzyme + NI) log2 means for each bile acid
baseline_means = df_T2 %>%
  filter(Inhibitor == "NI") %>%
  group_by(Bile_Acid) %>%
  summarise(
    Baseline_Mean = mean(Mean_Log2_Ratio, na.rm = TRUE)
  )

#inhibition effect
df_T2 = df_T2 %>%
  left_join(baseline_means, by = "Bile_Acid") %>%
  mutate(
    Log2_Difference = ifelse(
      Inhibitor != "NI",
      Mean_Log2_Ratio-Baseline_Mean ,
      NA
    )
  ) %>%
  mutate(
    Inhibition_Level = case_when(
      Inhibitor == "NI" ~ "Control (NI)",
      Log2_Difference <= 0 ~ "Strong Inhibition",
      Log2_Difference > 0 & Log2_Difference < 0.5 ~ "Moderate Inhibition",
      Log2_Difference >= 0.5 & Log2_Difference < 1.0 ~ "Weak Inhibition",
      Log2_Difference >= 1.0 ~ "No Inhibition",
      TRUE ~ "Undefined"
    )
  )

colors = c(
  "Control (NI)"        = "#1F78B4",  
  "No Inhibition"       = "#1F78B4", 
  "Weak Inhibition"     = "#6BAED6",  
  "Moderate Inhibition" = "#F08080",
  "Strong Inhibition"   = "#E31A1C"   
)



ggplot(df_T2, aes(x = Bile_Acid, y = Mean_Log2_Ratio, fill = Inhibition_Level)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = colors, name = "Inhibition Level") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~Inhibitor, ncol = 5) +
  labs(title = "Enzyme T2 – Log₂ Fold Change in Bile Acids vs Ct+NI",
       y = "log₂(Ct+NI/Treated)",
       x = "Bile Acid") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

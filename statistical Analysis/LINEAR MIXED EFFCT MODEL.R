rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

install.packages("ggplot2")

d = read.csv("/Users/keerthim/Documents/ST4060/ST6090/CLUSTERING DATA/REPLICATE DATA.csv",header=TRUE)

head(d)
dim(d)
str(d)

d$Replicate= as.factor(d$Replicate)
d$Run..=as.factor(d$Run..)
d$Inhibitor=as.factor(d$Inhibitor)
d$BSH.enzyme=as.factor(d$BSH.enzyme)

bile_acids = c("TUDCA", "TCDCA", "TDCA", "TCA", "TLCA",
               "GUDCA", "GCDCA", "GDCA", "GCA", "GLCA")
names(d) <- c(
  "Run",
  "BSH_Enzyme",
  "Inhibitor",
  "Replicate",
  "Sample_ID",
  "TUDCA",
  "TCDCA",
  "TDCA",
  "TCA",
  "TLCA",
  "GUDCA",
  "GCDCA",
  "GDCA",
  "GCA",
  "GLCA"
)

enzyme_colors = c(
  "JL885" = "blue", "1011c" = "brown", "A" = "darkgreen",
  "B" = "red", "F" = "purple", "T2" = "magenta", "Ct" = "#e377c2"
)
#########################INHIBITOR VERSUS RUN FOR EACH BILE ACID########################################
#############TUDCA##########
plot_data_TUDCA = d %>%
  group_by(Run, Replicate, Inhibitor, BSH_Enzyme) %>%
  summarise(TUDCA = first(TUDCA), .groups = "drop") %>%
  
  
  pivot_wider(
    names_from = Inhibitor,
    values_from = TUDCA,
    names_prefix = "TUDCA_"
  ) %>%
  
  
  filter(!is.na(TUDCA_NI)) %>%
  
  pivot_longer(
    cols = starts_with("TUDCA_") & !matches("TUDCA_NI"), # Exclude NI from reshaping
    names_to = "Inhibitor",
    values_to = "TUDCA_with_inhib",
    names_prefix = "TUDCA_"
  ) %>%
  
  filter(!is.na(TUDCA_with_inhib)) %>%
  
  bind_rows(
    d %>%
      filter(Inhibitor == "NI") %>%
      group_by(Run, Replicate, Inhibitor, BSH_Enzyme) %>%
      summarise(TUDCA_with_inhib = first(TUDCA), .groups = "drop") %>%
      mutate(TUDCA_NI = TUDCA_with_inhib, Inhibitor = "NI")
  ) %>%
  
  filter(Inhibitor != "VI") %>%
  
  mutate(
    Reference_NI = ifelse(Inhibitor == "NI", TUDCA_with_inhib, TUDCA_NI)
  )

##PLOT
ggplot(plot_data_TUDCA, aes(x = Reference_NI, y = TUDCA_with_inhib)) +
  geom_point(
    size = 3,
    aes(color = BSH_Enzyme, shape = Replicate),
    alpha = 0.8
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray30") +
  facet_wrap(~Inhibitor) +
  
  # Log scales
  scale_x_log10(
    breaks = c(0.01, 0.1, 1, 10, 100),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10, 100),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  
  scale_color_manual(values = enzyme_colors, name = "BSH Enzyme") +
  
  labs(
    x = "No Inhibitor (NI) Reference Value",
    y = "TUDCA Activity With Inhibitor",
    title = "Inhibitor Effects by BSH Enzyme Variant"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.box.spacing = unit(0.5, "cm"),
    panel.grid.major = element_line(color = "gray90"),
    strip.background = element_rect(fill = "gray95", color = NA)
  ) +
  guides(
    color = guide_legend(nrow = 1, override.aes = list(size = 4)),
    shape = guide_legend(title = "Replicate", override.aes = list(size = 4))
  )



################DATA TRANSFORMATION   : transform all the concentration to log2 ratio##################


df_log2 = d
epsilon = 1e-10

for (i in seq_len(nrow(df_log2))) {
  enzyme = df_log2$BSH_Enzyme[i]
  inhibitor = df_log2$Inhibitor[i]
  run_id = as.character(df_log2$Run[i])
  replicate_id = as.character(df_log2$Replicate[i])
  
  for (acid in bile_acids) {
    original_value = df_log2[[acid]][i]
    original_value = ifelse(original_value == 0, epsilon, original_value)
    
    if (enzyme == "Ct" && inhibitor == "NI") {
      df_log2[[paste0("log2_", acid)]][i] = log2(original_value)
      
    } else if (inhibitor == "NI" && enzyme != "Ct") {
      match_row = df_log2$BSH_Enzyme == "Ct" &
        df_log2$Inhibitor == "NI" &
        as.character(df_log2$Run) == run_id &
        as.character(df_log2$Replicate) == replicate_id
      
      if (sum(match_row) == 1) {
        numerator = df_log2[[acid]][match_row]
        numerator = ifelse(numerator == 0, epsilon, numerator)
        denominator = original_value
        df_log2[[paste0("log2_", acid)]][i] = log2(numerator / denominator)
      } else {
        df_log2[[paste0("log2_", acid)]][i] = NA
      }
      
    } else if (inhibitor != "NI" && enzyme != "Ct") {
      match_row = df_log2$BSH_Enzyme == enzyme &
        df_log2$Inhibitor == "NI" &
        as.character(df_log2$Run) == run_id &
        as.character(df_log2$Replicate) == replicate_id
      
      if (sum(match_row) == 1) {
        numerator = original_value
        denominator = df_log2[[acid]][match_row]
        denominator = ifelse(denominator == 0, epsilon, denominator)
        df_log2[[paste0("log2_", acid)]][i] = log2(numerator / denominator)
      } else {
        df_log2[[paste0("log2_", acid)]][i] = NA
      }
      
    } else {
      df_log2[[paste0("log2_", acid)]][i] = NA
    }
  }
}

###REMOVE THE BILE ACID COLUMNS
df_log2 = df_log2[, !(names(df_log2) %in% bile_acids)]

##############reshape the data#################
library(dplyr)

df_list <- list()

for (ba in bile_acids) {
  log2_col <- paste0("log2_", ba)
  
  df_sub <- df_log2 %>%
    dplyr::select(Run, BSH_Enzyme, Inhibitor, Replicate, Sample_ID, all_of(log2_col)) %>%
    dplyr::rename(log2_concentration = all_of(log2_col)) %>%
    dplyr::mutate(Bile_Acid = ba)
  
  df_list[[ba]] <- df_sub
}


df_combined <- bind_rows(df_list)
df_combined$Bile_Acid=as.factor(df_combined$Bile_Acid)
df_combined$Sample_ID=NULL


##########remove inhibitor VI############
df_combined=df_combined[df_combined$Inhibitor !='VI',]
df_combined$Inhibitor = droplevels(df_combined$Inhibitor)
str(df_combined)
  
df_combined$BSH_Enzyme=relevel(df_combined$BSH_Enzyme,ref="Ct")
df_combined$Inhibitor=relevel(df_combined$Inhibitor,ref="NI")
df_combined$Bile_Acid=relevel(df_combined$Bile_Acid,ref="TUDCA")

####sanity check for any NA /inf values
clean_data <- df_combined %>%
  filter(!is.na(log2_concentration),
         is.finite(log2_concentration)) %>%
  mutate(across(where(is.factor), droplevels)) %>%
  group_by(BSH_Enzyme, Inhibitor) %>%
  filter(n() > 0) %>%
  ungroup()


 anti_join(df_combined, clean_data)
 
 
 #############################PLOT THE TRANSFORMED DATA ##########################
 
 #Get reference (no inhibitor) log2 concentrations for each enzyme/run/rep
 ref_data_TUDCA <- df_combined %>%
   filter(Bile_Acid == "TUDCA", Inhibitor == "NI") %>%
   select(Run, Replicate, BSH_Enzyme, log2_concentration) %>%
   rename(log2_conc_NI = log2_concentration)
 
 #Join reference with all inhibitor data (except VI)
 df_plot_TUDCA <- df_combined %>%
   filter(Bile_Acid == "TUDCA", Inhibitor != "VI") %>%
   left_join(ref_data_TUDCA, by = c("Run", "Replicate", "BSH_Enzyme")) %>%
   filter(!(Inhibitor == "NI" & BSH_Enzyme == "Ct")) %>%
   filter(!is.na(log2_conc_NI))
 
 #Plot (log2 conc vs log2 conc)
 ggplot(df_plot_TUDCA, aes(x = log2_conc_NI, y = log2_concentration)) +
   geom_hline(yintercept = 0, linetype = "dotted", color = "gray30") +
   
   geom_point(aes(color = BSH_Enzyme, shape = Replicate), size = 3, alpha = 0.85) +
   
   facet_wrap(~Inhibitor, scales = "fixed") +
   
   scale_x_continuous(
     name = expression(Log[2]*" Ratio (Ct+NI / Enzyme+NI) TUDCA Concentration (No Inhibitor - NI)"),
     breaks = seq(-4, 16, by = 4),
     limits = c(-4, 16)
   ) +
   scale_y_continuous(
     name = expression(Log[2]*" Ratio (Enzyme+Inhibitor / Enzyme+NI) TUDCA Concentration (With Inhibitor)"),
     breaks = seq(-4, 16, by = 4),
     limits = c(-4, 16)
   ) +
   
   scale_color_manual(values = enzyme_colors) +
   
   labs(
     title = expression(Log[2]~"TUDCA Ratios: Effect of Inhibitors on BSH Enzyme Activity"),
     caption = "Above 0: Inhibition | Below 0: No Inhibition | At 0: No Effect"
   ) +
   
   theme_minimal(base_size = 13) +
   theme(
     axis.text.x = element_text(angle = 45, hjust = 1),
     axis.text.y = element_text(angle = 45),
     legend.position = "bottom",
     strip.background = element_rect(fill = "gray95", color = NA),
     plot.caption = element_text(size = 10, hjust = 0.5, face = "italic", color = "gray30")
   ) +
   
   guides(
     color = guide_legend(nrow = 1, override.aes = list(size = 4)),
     shape = guide_legend(title = "Replicate", override.aes = list(size = 4))
   )
 
 
###############Linear mixed model######################
 

library(lme4)
library(lmerTest)

set.seed(6090)
###RUN AND REPLICATE AS RANDOM EFFECT , ENZYME_INIBITOR(main effect plus interaction) AS FIXED EFFECT##############

############model_1 LINEAR MODEL##################
model_1=lm(log2_concentration ~ BSH_Enzyme + Inhibitor + Bile_Acid
           ,data = df_combined)
summary(model_1)

## model is not account for interactions

##############MODEL 2 ENZYME , INHIBITOR INTERACTION###############
model_2=lm(log2_concentration ~ BSH_Enzyme * Inhibitor + Bile_Acid ,data = df_combined)
summary(model_2)
library(tidyverse)

MISSING=df_combined %>%
  count(BSH_Enzyme, Inhibitor) %>%
  complete(BSH_Enzyme, Inhibitor, fill = list(n = 0)) %>%
  filter(n == 0)

print(MISSING,n=22)

###RANK DEFICINECY ############
X <- model.matrix(~ BSH_Enzyme * Inhibitor, data = df_combined)
cat("Rank of design matrix:", qr(X)$rank, "/", ncol(X), "columns")



##############BETTER MODEL ACCOUNTING FOR THE RAMDOM EFFECTS######################

######LINEAR MIXED EFFECT MODEL (ADDICTIVE MODEL)
model_3 = lmer(log2_concentration ~ BSH_Enzyme + Inhibitor + Bile_Acid +
                 (1|Run) + (1|Run:Replicate), 
               data = df_combined)

summary(model_3)$varcor
lme4::isSingular(model_3)
summary(model_3)
######LINEAR MIXED EFFECT MODEL WITH INTERACTIONS( ENZYME INHIBITOR INTERACTION)###############
model_4 = lmer(log2_concentration ~ BSH_Enzyme * Inhibitor + Bile_Acid +
                (1|Run) + (1|Run:Replicate), 
              data = df_combined)

summary(model_4)$varcor
lme4::isSingular(model_4)
summary(model_4)
######comparing model 3 and model 4########
anova(model_3,model_4)

######LINEAR MIXED EFFECT MODEL WITH INTERACTIONS(INHIBITOR BILE ACID INTERACTION)###############
model_5 = lmer(log2_concentration ~ BSH_Enzyme + Inhibitor * Bile_Acid +
                 (1|Run) + (1|Run:Replicate), 
               data = df_combined)

summary(model_5)$varcor
lme4::isSingular(model_5)
summary(model_5)

######LINEAR MIXED EFFECT MODEL WITH INTERACTIONS( ENZYME Bile Acid INTERACTION)###############
model_6 = lmer(log2_concentration ~ BSH_Enzyme * Bile_Acid + Inhibitor  +
                 (1|Run) + (1|Run:Replicate), 
               data = df_combined)

summary(model_6)$varcor
lme4::isSingular(model_6)
summary(model_6)

#########compare model##############INTERACTION TERMS###########################
anova(model_3,model_5)
anova(model_3,model_6)
anova(model_3,model_4,model_5,model_6)

#############COMBINE BOTH INTERACTION TERMS WHICH WE GOT SIGNIFICANT#########
model_7 = lmer(log2_concentration ~ BSH_Enzyme * Bile_Acid + BSH_Enzyme * Inhibitor  +
                 (1|Run) + (1|Run:Replicate), 
               data = df_combined)

summary(model_7)$varcor
lme4::isSingular(model_7)
summary(model_7)

anova(model_3,model_4,model_6,model_7)
anova(model_4,model_6,model_7)
############### run as random effect , replicate as a dummy variable ############
##Replicates are not interchangeable across runs.
###Run-to-run differences matter (e.g., reagent batches, day of experiment).
df_model_8 = df_combined %>%
  mutate(Replicate = factor(Replicate))

model_8 = lmer(log2_concentration ~ BSH_Enzyme * Inhibitor + Replicate + Bile_Acid +
                  (1|Run),  
                data = df_combined)

summary(model_8)$varcor
lme4::isSingular(model_8)
summary(model_8)

######comparing model 1 and model 3########
anova(model_4,model_7,model_8)

###lets save df


################model 9  enzyme , inhibitor, bile acid interaction SATURATED MODEL####################
model_9 = lmer(log2_concentration ~ BSH_Enzyme * Inhibitor * Bile_Acid +
                 (1|Run) + (1|Run:Replicate), 
               data = df_combined)

summary(model_9)$varcor
lme4::isSingular(model_9)
summary(model_9)

######comparing the models
result=anova(model_9,model_4,model_7)
result
#######MODEL 10 #######################
df_combined$BSH_Enzyme = relevel(as.factor(df_combined$BSH_Enzyme), ref = "Ct")
df_combined$Inhibitor  = relevel(as.factor(df_combined$Inhibitor), ref = "NI")
df_combined$Bile_Acid  = relevel(as.factor(df_combined$Bile_Acid), ref = "TUDCA")

model_10 = lmer(log2_concentration ~ (BSH_Enzyme * Inhibitor)+ (BSH_Enzyme *  Bile_Acid) +
                  (Inhibitor *  Bile_Acid) +
                 (1|Run) + (1|Run:Replicate), 
               data = df_combined)

summary(model_10)$varcor
lme4::isSingular(model_10)
summary(model_10)

library(broom.mixed)

fixed_effects <- tidy(model_10, effects = "fixed")
dropped <- fixed_effects[is.na(fixed_effects$estimate), ]
print(dropped$term)

#######################MODEL 10 COMPARISION###########################
result=anova(model_9,model_7,model_10)
result
#Get the full summary
model_summary <- summary(model_9)
coef_table <- as.data.frame(model_summary$coefficients)
coef_table$Term <- rownames(coef_table)

#Filter for significant terms
sig_coefs <- coef_table %>%
  filter(`Pr(>|t|)` < 0.05) %>%
  select(Term, Estimate, `Std. Error`, `Pr(>|t|)`)
print(sig_coefs)

write.csv(sig_coefs, "/Users/keerthim/Documents/ST4060/ST6090/R CODE/significant_coefficients_model_4.csv", row.names = FALSE)

#############VISUALIZING THE INTERACTIONS OF MODEL 7####################
library(emmeans)
library(ggplot2)
library(dplyr)
emm_inhib <- emmeans(model_10, ~ BSH_Enzyme * Inhibitor,pbkrtest.limit = 4950)
emm_bile <- emmeans(model_10, ~ BSH_Enzyme * Bile_Acid,pbkrtest.limit = 4950)
emm_df1 <- as.data.frame(emm_inhib)

ggplot(emm_df1, aes(x = Inhibitor, y = emmean, color = BSH_Enzyme, group = BSH_Enzyme)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  labs(
    title = "Interaction: BSH Enzyme × Inhibitor",
    x = "Inhibitor",
    y = expression(Log[2]*" Concentration")
  ) +
  theme_minimal()

###############
library(dplyr)
library(ggplot2)

#Calculate average effect per inhibitor
broad_summary <- emm_df1 %>%
  group_by(Inhibitor) %>%
  summarise(avg_effect = mean(emmean, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(avg_effect))

#Select top 3 broad-spectrum inhibitors
broad_spectrum_inhibs <- broad_summary %>%
  slice_max(avg_effect, n = 3) %>%
  pull(Inhibitor)

#Identify strongest inhibitor per enzyme
top_inhibs <- emm_df1 %>%
  group_by(BSH_Enzyme) %>%
  filter(emmean == max(emmean)) %>%
  ungroup()

#Plot
ggplot(emm_df1, aes(x = Inhibitor, y = emmean, group = BSH_Enzyme, color = BSH_Enzyme)) +

  geom_line() +
  geom_point(size = 2) +

  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +

  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.7) +
  geom_point(
    data = subset(emm_df1, Inhibitor %in% broad_spectrum_inhibs),
    aes(x = Inhibitor, y = emmean),
    shape = 21, stroke = 1.2, fill = NA, color = "black", size = 3,
    inherit.aes = FALSE
  ) +
  
  geom_text(
    data = top_inhibs,
    aes(label = round(emmean, 2)),
    vjust = -1.2,
    size = 3,
    show.legend = FALSE
  ) +
  
  labs(
    title = "Interaction: BSH Enzyme × Inhibitor",
    subtitle = "Dashed line at y = 0 → No effect. Above 0 = Inhibition, Below 0 = No Inhibition",
    x = "Inhibitor",
    y = expression(Log[2]*" Concentration Change (vs NI)")
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )


############ENZYME BILE ACID INTERACTION #####################
emm_df2 <- as.data.frame(emm_bile)

ggplot(emm_df2, aes(x = Bile_Acid, y = emmean, color = BSH_Enzyme, group = BSH_Enzyme)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  labs(
    title = "Interaction: BSH Enzyme × Bile Acid",
    x = "Bile Acid",
    y = expression(Log[2]*" Concentration")
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


####################SHOW THIS PLOT FOR INTERACTION MODEL 7###################
library(dplyr)
library(ggplot2)

# Categorize effect direction for legend
emm_df1_filtered <- emm_df1 %>%
  filter(BSH_Enzyme != "Ct", Inhibitor != "NI") %>%
  filter(!is.na(emmean)) %>%
  mutate(
    Effect = case_when(
      emmean > 0 ~ "Above 0 (Inhibition)",
      emmean < 0 ~ "Below 0 (No Inhibition)",
      TRUE ~ "No Effect"
    )
  )


# Plot
ggplot(emm_df1_filtered , aes(x = Inhibitor, y = emmean, color = Effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.6) +
  
  geom_point(size = 2) +
  geom_line(aes(group = BSH_Enzyme), linewidth = 0.6) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.25, color = "gray40") +
  
  facet_wrap(~ BSH_Enzyme) +
  
  scale_color_manual(
    values = c(
      "Above 0 (Inhibition)" = "firebrick",
      "Below 0 (No Inhibition)" = "darkgreen",
      "No Effect" = "black"
    ),
    name = "Effect Type"
  ) +
  
  labs(
    title = "Estimated Interaction: BSH Enzyme × Inhibitor",
    subtitle = "Dashed line at y = 0 shows baseline (NI). Above = Inhibition. Below = No inhibition.",
    x = "Inhibitor",
    y = expression(Log[2]*" Concentration (Estimated)")
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )


emm_df2 <- as.data.frame(emmeans(model_7, ~ BSH_Enzyme * Bile_Acid,pbkrtest.limit = 4950))

# Cleaned data remove NA and "Ct"
emm_df2_filtered <- emm_df2 %>%
  filter(!is.na(emmean)) %>%
  filter(BSH_Enzyme != "Ct")

# Plot
ggplot(emm_df2_filtered, aes(x = Bile_Acid, y = emmean)) +
  geom_point(color = "firebrick", size = 2) +
  geom_line(aes(group = 1), color = "firebrick") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  facet_wrap(~ BSH_Enzyme) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    title = "Estimated Interaction: BSH Enzyme × Bile Acid",
    subtitle = "Dashed line at 0: No change; ↑ Log2 = ↑ inhibition (↓ enzyme activity)",
    x = "Bile Acid",
    y = expression(Log[2]*" Concentration (Estimated)")
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "gray95", color = NA),
    strip.text = element_text(face = "bold")
  )

##################PREDICTED VERSUS OBSERVED#########################
#######predicted versus observed MODEL 7#####################
df_combined$fitted_values_7 <- as.numeric(predict(model_7))

ggplot(df_combined, aes(x = fitted_values_7, y = log2_concentration)) +
  geom_point(alpha = 0.6, color = "darkblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  annotate("text", x = min(df_combined$fitted_values_7, na.rm = TRUE),
           y = max(df_combined$log2_concentration, na.rm = TRUE),
           label = "y = x (Perfect Prediction)",
           hjust = 0, vjust = 1.2, size = 4, color = "gray30") +
  labs(
    title = expression("Observed vs Predicted " * log[2] * " Concentration"),  
    subtitle = "Dashed line: y = x — perfect prediction line",
    x = expression("Predicted " * log[2] * " Concentration"),
    y = expression("Observed " * log[2] * " Concentration")
  )+
  theme_minimal()

###########MODEL 9
df_combined$fitted_values_9 <- as.numeric(predict(model_9))

ggplot(df_combined, aes(x = fitted_values_9, y = log2_concentration)) +
  geom_point(alpha = 0.6, color = "darkblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  annotate("text", x = min(df_combined$fitted_values_9, na.rm = TRUE),
           y = max(df_combined$log2_concentration, na.rm = TRUE),
           label = "y = x (Perfect Prediction)",
           hjust = 0, vjust = 1.2, size = 4, color = "gray30") +
  labs(
    title = expression("Observed vs Predicted " * log[2] * " Concentration"),  
    subtitle = "Dashed line: y = x — perfect prediction line",
    x = expression("Predicted " * log[2] * " Concentration"),
    y = expression("Observed " * log[2] * " Concentration")
  )+
  theme_minimal()

#################MODEL 10 ##################
df_combined$fitted_values_10 <- predict(model_10)

ggplot(df_combined, aes(x = fitted_values_10, y = log2_concentration)) +
  geom_point(alpha = 0.6, color = "darkblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  annotate("text", x = min(df_combined$fitted_values, na.rm = TRUE),
           y = max(df_combined$log2_concentration, na.rm = TRUE),
           label = "y = x (Perfect Prediction)",
           hjust = 0, vjust = 1.2, size = 4, color = "gray30") +
  labs(
    title = expression("Observed vs Predicted " * log[2] * " Concentration"),  
    subtitle = "Dashed line: y = x — perfect prediction line",
    x = expression("Predicted " * log[2] * " Concentration"),
    y = expression("Observed " * log[2] * " Concentration")
  )+
  theme_minimal()


#####################MODEL DIAGNOSTICS##################################
######RESIDUAL VERSUS FITTED
plot(resid(model_7) ~ fitted(model_7),
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, lty = 2, col = "gray")

plot(resid(model_9) ~ fitted(model_9),
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, lty = 2, col = "gray")

plot(resid(model_10) ~ fitted(model_10),
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, lty = 2, col = "gray")

########Q-Q PLOT###############
qqnorm(resid(model_7))
qqline(resid(model_7), col = "red")
hist(resid(model_7)) 
ranef(model_7)  

qqnorm(resid(model_9))
qqline(resid(model_9), col = "red")
hist(resid(model_9)) 
ranef(model_9)  

qqnorm(resid(model_10))
qqline(resid(model_10), col = "red")

install.packages(c("gtable", "performance", "see", "patchwork"))
library(performance)
library(see)


check_model(model_10)
check_model(model_7)

############FITTED VALUES MODEL 7#######################


#########extract the fitted value from model 7##############

##########add the standard error############# STANDARD ERROR AND FITTED VALUE#############
library(merTools)
library(dplyr)

pred_ci_fitted <- predictInterval(
  model_7,
  newdata = df_combined,
  n.sims = 1000,
  level = 0.95,
  include.resid.var = FALSE,
  stat = "mean"
)

pred_ci_fitted <- pred_ci_fitted %>%
  mutate(fitted_se = (upr - lwr) / (2 * 1.96))

df_combined <- df_combined %>%
  mutate(
    fitted_values_7 = pred_ci_fitted$fit,
    fitted_se_7 = pred_ci_fitted$fitted_se,
    ci_lwr_7 = pred_ci_fitted$lwr,
    ci_upr_7 = pred_ci_fitted$upr
  )

df_combined

#############overall plot tudca############
library(MASS)
library(dplyr)  
library(tidyr)
fitted_TUDCA = df_combined %>%
  filter(Bile_Acid == "TUDCA", Inhibitor != "VI", BSH_Enzyme != "Ct")

ref_TUDCA = df_combined %>%
  filter(Bile_Acid == "TUDCA", Inhibitor == "NI", BSH_Enzyme != "Ct") %>%
  dplyr::select(Run, Replicate, BSH_Enzyme, log2_concentration) %>%
  rename(log2_conc_NI = log2_concentration)

fitted_TUDCA_plot = fitted_TUDCA %>%
  left_join(ref_TUDCA, by = c("Run", "Replicate", "BSH_Enzyme")) %>%
  filter(!is.na(log2_conc_NI))

library(ggplot2)
observed_df <- fitted_TUDCA_plot |> 
  mutate(PointType = "Observed")

fitted_df <- fitted_TUDCA_plot |>
  mutate(
    log2_concentration = fitted_values_7,
    PointType = "Fitted"
  )

plot_df <- rbind(observed_df, fitted_df)

ggplot(plot_df, aes(x = log2_conc_NI, y = log2_concentration)) +
  geom_abline(
    intercept = 0, slope = 1,
    linetype = "dotted",
    color = "black",
    size = 1,
    show.legend = TRUE
  ) +
  geom_smooth(
    data = fitted_df,
    method = "lm",
    aes(linetype = "Regression (Fitted)", color = "Regression (Fitted)"),
    se = FALSE,
    size = 1
  ) +
  geom_point(
    aes(color = PointType),
    size = 2,
    alpha = 0.8
  ) +
  scale_color_manual(
    name = "Data",
    values = c(
      "Observed" = "black",
      "Fitted" = "blue",
      "Regression (Fitted)" = "blue"
    )
  )+ 
  
  labs(
    x = expression(Log[2]*" TUDCA Concentration (No Inhibitor – NI)"),
    y = expression(Log[2]*" TUDCA Concentration (With Inhibitor)"),
    title = "Observed vs Fitted TUDCA Values with Identity & Model Regression Lines"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )

############CONFIDENCE INTERVAL#############

#########CI for the fitted values###########
library(merTools)

###for fitted data
fitted_ci <- predictInterval(model_7, 
                             newdata = df_combined,
                             level = 0.95,
                             n.sims = 1000,include.resid.var = FALSE)
pred_ci <- predictInterval(model_7, include.resid.var = TRUE)

mean(fitted_ci$upr - fitted_ci$lwr)  # Narrower
mean(pred_ci$upr - pred_ci$lwr)      # Wider


###############################PLOT THE CI FOR TUDCA########################
library(merTools)
library(ggplot2)
library(dplyr)

tudca_data = df_combined %>%
  filter(Bile_Acid == "TUDCA", Inhibitor != "VI", BSH_Enzyme != "Ct")

tudca_ci = predictInterval(
  model_7,
  newdata = tudca_data,
  n.sims = 1000,
  include.resid.var = FALSE
) %>%
  bind_cols(tudca_data)

##FOR EACH INHIBITOR
tudca_summary <- tudca_ci %>%
  group_by(Inhibitor) %>%
  summarise(
    avg_fit = mean(fit),
    avg_upr = mean(upr),
    avg_lwr = mean(lwr),
    n = n()
  ) %>%
  ungroup()

#############CI FOR EACH INHIBITOR#############
print(tudca_summary,n=20)

###INHIBITOR VERSUS NO INHIBITOR###########
tudca_summary <- tudca_ci %>%
  mutate(Inhibitor_Status = ifelse(Inhibitor == "NI", "No Inhibitor", "With Inhibitor")) %>%
  group_by(Inhibitor_Status) %>%
  summarise(
    avg_fit = mean(fit),
    avg_upr = mean(upr),
    avg_lwr = mean(lwr)
  )

print(tudca_summary)
with_inhib_avg <- tudca_summary %>% 
  filter(Inhibitor_Status == "With Inhibitor")

ggplot(plot_df, aes(x = log2_conc_NI, y = log2_concentration)) +
  # Reference line (y=x)
  geom_abline(
    intercept = 0, slope = 1,
    linetype = "dotted",
    color = "black",
    linewidth = 1
  ) +
  
  # Regression line
  geom_smooth(
    data = fitted_df,
    method = "lm",
    aes(linetype = "Regression (Fitted)", color = "Regression (Fitted)"),
    se = FALSE,
    linewidth = 1
  ) +
  
  # Points
  geom_point(
    aes(color = PointType),
    size = 2,
    alpha = 0.8
  ) +
  
  # CI lines for With Inhibitor
  geom_hline(
    data = with_inhib_avg,
    aes(yintercept = avg_fit, linetype = "Avg (With Inhibitor)"),
    color = "red3",
    linewidth = 1
  ) +
  geom_hline(
    data = with_inhib_avg,
    aes(yintercept = avg_upr, linetype = "95% CI"),
    color = "red3",
    linewidth = 0.6
  ) +
  geom_hline(
    data = with_inhib_avg,
    aes(yintercept = avg_lwr, linetype = "95% CI"),
    color = "red3",
    linewidth = 0.6
  )+
  
  # Text label
  geom_text(
    data = with_inhib_avg,
    aes(
      x = max(plot_df$log2_conc_NI),
      y = avg_fit,
      label = paste0("Avg = ", round(avg_fit, 2))
    ),
    color = "red3",
    hjust = 1.1,
    vjust = -0.5,
    size = 3
  )+
  
  # Color and line type scales
  scale_color_manual(
    name = "Data",
    values = c(
      "Observed" = "black",
      "Fitted" = "blue",
      "Regression (Fitted)" = "blue"
    )
  ) +
  scale_linetype_manual(
    name = "Reference Lines",
    values = c(
      "Regression (Fitted)" = "solid",
      "Avg (With Inhibitor)" = "solid",
      "95% CI" = "dashed"
    ),
    guide = guide_legend(
      override.aes = list(
        color = c("red3", "red3", "blue"),
        linewidth = c(1, 1, 0.6)
      )
    )
  ) +
  scale_y_continuous(
    name = expression(Log[2]*" TUDCA Concentration (With Inhibitor)"),
    breaks = c(
      pretty(range(plot_df$log2_concentration), n = 5),
      round(with_inhib_avg$avg_fit,3),
      round(with_inhib_avg$avg_lwr,3),
      round(with_inhib_avg$avg_upr,3)
    ),
    labels = function(x) {
      ifelse(
        x %in% round(c(with_inhib_avg$avg_fit, with_inhib_avg$avg_lwr, with_inhib_avg$avg_upr), 3),
        paste0(x, "  ◀"),
        x
      )
    }
  )+
  
  labs(
    x = expression(Log[2]*" TUDCA Concentration (No Inhibitor – NI)"),
    y = expression(Log[2]*" TUDCA Concentration (With Inhibitor)"),
    title = "Observed vs Fitted TUDCA Values with Identity & Model Regression Lines"
  )+ 
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(0.2, "cm")
  )


###################BETTER CI PLOT FOR ENZYME AND BILE ACID##########################
library(dplyr)

fitted_ci_df = bind_cols(df_combined, fitted_ci)
summary_df <- fitted_ci_df %>%
  group_by(BSH_Enzyme, Inhibitor, Bile_Acid) %>%
  summarise(
    mean_fit = mean(fit),
    lwr = mean(lwr),
    upr = mean(upr),
    .groups = "drop"
  )

library(ggplot2)

plot_by_enzyme <- function(enzyme_name) {
  enzyme_df <- summary_df %>% filter(BSH_Enzyme == enzyme_name)
  
  ggplot(enzyme_df, aes(x = Bile_Acid, y = mean_fit, fill = Bile_Acid)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = lwr, ymax = upr),
                  position = position_dodge(width = 0.8), width = 0.2) +
    facet_wrap(~ Inhibitor) +
    labs(
      title = paste("Estimated BSH Activity for Enzyme:", enzyme_name),
      y = expression(Log[2]*" Fitted Concentration"),
      x = "Bile Acid"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

for (e in unique(summary_df$BSH_Enzyme)) {
  print(plot_by_enzyme(e))
}

#############overall PLOT for CI for the entire dataset#################
library(dplyr)
library(merTools)
library(ggplot2)

# Reference (NI): exclude Ct
ref_NI <- df_combined %>%
  filter(Inhibitor == "NI", BSH_Enzyme != "Ct") %>%
  dplyr::select(Run, Replicate, BSH_Enzyme, Bile_Acid, log2_concentration) %>%
  rename(log2_conc_NI = log2_concentration)
with_inhib <- df_combined %>%
  filter(Inhibitor != "NI", Inhibitor != "VI", BSH_Enzyme != "Ct") %>%
  dplyr::select(Run, Replicate, BSH_Enzyme, Bile_Acid, Inhibitor, log2_concentration)

# Join reference + with-inhibitor observed
plot_df <- inner_join(with_inhib, ref_NI,
                      by = c("Run", "Replicate", "BSH_Enzyme", "Bile_Acid")) %>%
  mutate(PointType = "Observed")

fitted_ci <- predictInterval(
  model_10, 
  newdata = df_combined,
  n.sims = 1000,
  include.resid.var = FALSE
)

df_combined_pred <- bind_cols(df_combined, fitted_ci)

# Fitted data: with inhibitor only
fitted_df <- df_combined_pred %>%
  filter(Inhibitor != "NI", Inhibitor != "VI", BSH_Enzyme != "Ct") %>%
  dplyr::select(Run, Replicate, BSH_Enzyme, Bile_Acid, Inhibitor, fit) %>%
  inner_join(ref_NI, by = c("Run", "Replicate", "BSH_Enzyme", "Bile_Acid")) %>%
  mutate(log2_concentration = fit, PointType = "Fitted")
plot_all <- bind_rows(plot_df, fitted_df)

with_inhib_avg <- df_combined_pred %>%
  filter(Inhibitor != "NI", Inhibitor != "VI", BSH_Enzyme != "Ct") %>%
  summarise(
    avg_fit = mean(fit),
    avg_lwr = mean(lwr),
    avg_upr = mean(upr)
  )

ggplot(plot_all, aes(x = log2_conc_NI, y = log2_concentration)) +
  # Identity line
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black", linewidth = 1) +
  
  # Regression line from fitted values
  geom_smooth(
    data = filter(plot_all, PointType == "Fitted"),
    method = "lm",
    aes(linetype = "Regression (Fitted)", color = "Regression (Fitted)"),
    se = FALSE,
    linewidth = 1
  ) +
  
  # Points
  geom_point(
    aes(color = PointType),
    size = 2,
    alpha = 0.6
  ) +
  
  # CI lines for With Inhibitor
  geom_hline(data = with_inhib_avg, aes(yintercept = avg_fit, linetype = "Avg (With Inhibitor)"),
             color = "red3", linewidth = 1) +
  geom_hline(data = with_inhib_avg, aes(yintercept = avg_upr, linetype = "95% CI"),
             color = "red3", linewidth = 0.6) +
  geom_hline(data = with_inhib_avg, aes(yintercept = avg_lwr, linetype = "95% CI"),
             color = "red3", linewidth = 0.6) +
  
  # Label the mean
  geom_text(data = with_inhib_avg,
            aes(x = max(plot_df$log2_conc_NI), y = avg_fit),
            label = "Mean (With Inhibitor)", color = "red3",
            hjust = 1.1, vjust = -0.5, size = 3) +
  
  # Custom scales
  scale_color_manual(
    name = "Data",
    values = c("Observed" = "black", "Fitted" = "blue", "Regression (Fitted)" = "blue")
  ) +
  scale_linetype_manual(
    name = "Reference Lines",
    values = c("Regression (Fitted)" = "solid", "Avg (With Inhibitor)" = "solid", "95% CI" = "dashed"),
    guide = guide_legend(
      override.aes = list(color = c("red3", "red3", "blue"), linewidth = c(1, 0.6, 1))
    )
  ) +
  scale_y_continuous(
    name = expression(Log[2]*" Concentration (With Inhibitor)"),
    breaks = c(
      pretty(range(plot_all$log2_concentration), n = 5),
      round(with_inhib_avg$avg_fit, 3),
      round(with_inhib_avg$avg_lwr, 3),
      round(with_inhib_avg$avg_upr, 3)
    ),
    labels = function(x) {
      ifelse(
        x %in% round(c(with_inhib_avg$avg_fit, with_inhib_avg$avg_lwr, with_inhib_avg$avg_upr), 3),
        paste0(x, "  ◀"),
        x
      )
    }
  ) +
  labs(
    x = expression(Log[2]*" Concentration (No Inhibitor – NI)"),
    y = expression(Log[2]*" Concentration (With Inhibitor)"),
    title = "Observed vs Fitted Values for All Bile Acids and Enzymes"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(0.2, "cm")
  )+
   annotate("text",
             x = max(plot_all$log2_conc_NI),
             y = min(plot_all$log2_concentration),
             label = "Interpretation:\nAbove 0 → Inhibition\nAt 0 → No Effect\nBelow 0 → No Inhibition",
             hjust = 1.1, vjust = -0.5,
             size = 4,
             fontface = "italic",
             color = "gray30")



#######plot for TUDCA#######
ggplot(fitted_TUDCA_plot, aes(x = log2_conc_NI, y = log2_concentration)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray30") +
  geom_point(aes(color = BSH_Enzyme, shape = Replicate), size = 3, alpha = 0.85) +
  geom_point(aes(y = fitted_value), shape = 18, color = "black", size = 3) +
  facet_wrap(~Inhibitor, scales = "free") +
  scale_x_continuous(
    name = expression(Log[2]*" Ratio (Ct+NI / Enzyme+NI) TUDCA Concentration (No Inhibitor - NI)")
  ) +
  scale_y_continuous(
    name = expression(Log[2]*" Ratio (Enzyme+Inhibitor / Enzyme+NI) TUDCA Concentration (With Inhibitor)")
  ) +
  scale_color_manual(values = enzyme_colors) +
  labs(
    title = expression(Log[2]~"TUDCA Ratios: Effect of Inhibitors on BSH Enzyme Activity")
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(angle = 45),
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray95", color = NA)
  ) +
  guides(
    color = guide_legend(nrow = 1, override.aes = list(size = 4)),
    shape = guide_legend(title = "Replicate", override.aes = list(size = 4))
  )



########better colours#####################
ggplot(fitted_TUDCA_plot,
       aes(x = log2_conc_NI, y = log2_concentration)) +
  
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray30") +
  
 
  geom_point(
    aes(colour = BSH_Enzyme, shape = Replicate),
    size  = 3,
    alpha = 0.7
  ) +
  
  geom_point(
    aes(y = fitted_value, shape = Replicate),
    colour = "black",
    size   = 3,      # a bit bigger
    stroke = 1.1,
    alpha  = 1,
    show.legend = FALSE   # avoid duplicate legend
  ) +
  
  geom_point(
    aes(y = fitted_value, colour = BSH_Enzyme, shape = Replicate),
    size   = 2,
    alpha  = 1,
    show.legend = FALSE
  ) +
  
  facet_wrap(~Inhibitor, scales = "free") +
  scale_colour_manual(values = enzyme_colors) +
  labs(
    x = expression(Log[2]*" Ratio (Ct+NI / Enzyme+NI) TUDCA Concentration (No Inhibitor – NI)"),
    y = expression(Log[2]*" Ratio (Enzyme+Inhibitor / Enzyme+NI) TUDCA Concentration (With Inhibitor)"),
    title = expression(Log[2]~"TUDCA Ratios: Effect of Inhibitors on BSH Enzyme Activity")
  )


############EXTRACT COEFFICIENTS VALUES###############
coef_df <- data.frame(
  term = names(fixef(model_7)),
  estimate = as.numeric(fixef(model_7))
)
coef_df

str(coef_df)

library(lme4)
library(broom.mixed)

set.seed(6090)
K = 5
n = nrow(df_combined)
folds = sample(rep(1:K, length.out = n))
cv_coefs = list()

for (k in 1:K) {
  i.train = which(folds != k)
  dat.train = df_combined[i.train, ]
  dat.test  = df_combined[-i.train, ]

  model_k = lmer(log2_concentration ~ 
                    BSH_Enzyme * Inhibitor + 
                    BSH_Enzyme * Bile_Acid + 
                    Inhibitor * Bile_Acid + 
                    (1|Run) + (1|Run:Replicate),
                  data = dat.train)
  
  coefs_k = broom.mixed::tidy(model_k, effects = "fixed")
  coefs_k$Fold = k
  cv_coefs[[k]] = coefs_k
}
cv_results = do.call(rbind, cv_coefs)

library(dplyr)

cv_summary <- cv_results %>%
  group_by(term) %>%
  summarise(mean_estimate = mean(estimate), .groups = "drop")
##########extract TUDCA value#########
library(emmeans)

#Get enzyme × TUDCA
enzyme_tudca_direct <- emmeans(model_7, ~ BSH_Enzyme * Bile_Acid,pbkrtest.limit = 4950)

enzyme_tudca_df <- as.data.frame(enzyme_tudca_direct) %>%
  filter(Bile_Acid == "TUDCA") %>%
  mutate(
    term = paste0("BSH_Enzyme", BSH_Enzyme, ":Bile_AcidTUDCA"),
    mean_estimate = emmean
  ) %>%
  dplyr::select(term, mean_estimate)
enzyme_tudca_df=na.omit(enzyme_tudca_df)

################NOT WORKING ##################
# 2. Get inhibitor × TUDCA

emm_inhib_bile <- emmeans(model_7, ~ Inhibitor * Bile_Acid,pbkrtest.limit = 4950)
emm_inhib_bile=as.data.frame(emm_inhib_bile)
tudca_inhib_means <- as.data.frame(emm_inhib_bile) %>%
  filter(Bile_Acid == "TUDCA") %>%
  mutate(
    term = paste0("Inhibitor", Inhibitor, ":Bile_AcidTUDCA"),
    mean_estimate = emmean
  ) %>%
  dplyr::select(term, mean_estimate)



#Combine and bind to cv_summary

cv_summary_extended = bind_rows(cv_summary,enzyme_tudca_df)
intercept_value = cv_summary_extended$mean_estimate[cv_summary_extended$term == "(Intercept)"]

get_coef <- function(term, coef_table) {
  if (is.na(term)) return(0)
  val <- coef_table$mean_estimate[coef_table$term == term]
  if (length(val) == 0) return(0)
  return(val)
}

df_combined <- df_combined %>%
  mutate(
    enzyme_term = if_else(BSH_Enzyme != "Ct", paste0("BSH_Enzyme", BSH_Enzyme), NA_character_),
    inhibitor_term = if_else(Inhibitor != "NI", paste0("Inhibitor", Inhibitor), NA_character_),
    bileacid_term = paste0("Bile_Acid", Bile_Acid),
    
    ei_term = if_else(!is.na(enzyme_term) & !is.na(inhibitor_term),
                      paste0(enzyme_term, ":", inhibitor_term), NA_character_),
    
    eb_term = if_else(!is.na(enzyme_term),
                      paste0(enzyme_term, ":", bileacid_term), NA_character_)
  ) %>%
  rowwise() %>%
  mutate(
    enzyme_baseline_7 = if_else(BSH_Enzyme == "Ct", 0, get_coef(enzyme_term, cv_summary_extended)),
    inhibitor_baseline_7 = if_else(Inhibitor == "NI", 0, 
                                 get_coef(paste0("Inhibitor", Inhibitor), cv_summary_extended)),
    
    # Calculate inhibitor potency as enzyme + enzyme:inhibitor interaction
    enzyme_inhibitor_interaction_7 = if_else(BSH_Enzyme == "Ct" | Inhibitor == "NI", 0,
                                           get_coef(ei_term, cv_summary_extended)),
    inhibitor_potency_7 = if_else(Inhibitor != "NI",
                                  enzyme_baseline_7 - enzyme_inhibitor_interaction_7,
                                  0),
    
    enzyme_bileacid_interaction_7 = if_else(BSH_Enzyme == "Ct", 0, get_coef(eb_term, cv_summary_extended))
  ) %>%
  ungroup() %>%
  dplyr::select(-enzyme_term, -inhibitor_term, -bileacid_term, -ei_term, -eb_term)



nrow(df_combined)  # Should be 4950


write.csv(df_combined, "/Users/keerthim/Documents/ST4060/ST6090/R CODE/df_combined_with_components.csv", row.names = FALSE)
write.csv(cv_summary_extended, "/Users/keerthim/Documents/ST4060/ST6090/R CODE/cv_summary_coefficients.csv", row.names = FALSE)


################INSIGHTS################################

library(ggplot2)
library(dplyr)


df_combined %>%
  filter(Inhibitor == "NI", BSH_Enzyme != "Ct") %>%
  group_by(BSH_Enzyme, Bile_Acid) %>%
  summarise(mean_fitted = mean(fitted_values_7, na.rm = TRUE), .groups = "drop") %>%
  group_by(Bile_Acid) %>%
  mutate(highlight = BSH_Enzyme == BSH_Enzyme[which.max(mean_fitted)]) %>%
  ggplot(aes(x = reorder(BSH_Enzyme, -mean_fitted), y = mean_fitted, fill = highlight)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "steelblue"), guide = "none") +
  labs(
    title = "Top BSH Enzyme per Bile Acid (Inhibitor = NI)",
    x = "BSH Enzyme", y = "Mean log₂ concentration"
  ) +
  facet_wrap(~ Bile_Acid, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#############INHIBITOR POTENCY PLOT#########
library(ggplot2)
library(dplyr)


df_combined %>%
  filter(Inhibitor != "NI") %>%  # Exclude reference inhibitor
  group_by(BSH_Enzyme, Inhibitor, Bile_Acid) %>%
  summarise(mean_potency = mean(inhibitor_potency_7, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = Inhibitor, y = BSH_Enzyme, fill = mean_potency)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    name = "Potency"
  ) +
  facet_wrap(~ Bile_Acid, ncol = 3) +
  labs(
    title = "Inhibitor Potency per Enzyme Across Bile Acids (Excluding NI)",
    x = "Inhibitor",
    y = "BSH Enzyme"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )



library(dplyr)
library(ggplot2)

df_combined %>%
  filter(Inhibitor != "NI") %>%  # Exclude the reference inhibitor
  group_by(BSH_Enzyme, Inhibitor, Bile_Acid) %>%
  summarise(mean_potency = mean(inhibitor_potency_7, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = reorder(BSH_Enzyme, -mean_potency), y = mean_potency)) +
  geom_col(fill = "darkorange") +
  facet_grid(Bile_Acid ~ Inhibitor, scales = "free_y") +
  labs(
    title = "Inhibitor Potency per Enzyme, Faceted by Bile Acid and Inhibitor",
    x = "BSH Enzyme",
    y = "Mean Inhibitor Potency"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 14)
  )

library(dplyr)
library(ggplot2)

df_combined %>%
  group_by(BSH_Enzyme, Inhibitor, Bile_Acid) %>%
  summarise(mean_concentration = mean(fitted_values_7, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = reorder(BSH_Enzyme, -mean_concentration), y = mean_concentration)) +
  geom_col(fill = "steelblue") +
  facet_grid(Bile_Acid ~ Inhibitor, scales = "free_y") +
  labs(
    title = "Mean log2 Bile Acid Concentration per Enzyme\nFaceted by Inhibitor and Bile Acid",
    x = "BSH Enzyme",
    y = "Mean log2 Concentration"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 14)
  )
############FITTED VALUES MODEL 9#######################


#########extract the fitted value from model 9##############

##########add the standard error############# STANDARD ERROR AND FITTED VALUE#############
library(merTools)
library(dplyr)

pred_ci_fitted <- predictInterval(
  model_9,
  newdata = df_combined,
  n.sims = 200,
  level = 0.95,
  include.resid.var = FALSE,
  stat = "mean"
)

pred_ci_fitted <- pred_ci_fitted %>%
  mutate(fitted_se = (upr - lwr) / (2 * 1.96))

df_combined <- df_combined %>%
  mutate(
    fitted_values_9 = pred_ci_fitted$fit,
    fitted_se_9 = pred_ci_fitted$fitted_se,
    ci_lwr_9 = pred_ci_fitted$lwr,
    ci_upr_9 = pred_ci_fitted$upr
  )

df_combined

#############overall plot tudca############
library(MASS)
library(dplyr)  
library(tidyr)
fitted_TUDCA = df_combined %>%
  filter(Bile_Acid == "TUDCA", Inhibitor != "VI", BSH_Enzyme != "Ct")

ref_TUDCA = df_combined %>%
  filter(Bile_Acid == "TUDCA", Inhibitor == "NI", BSH_Enzyme != "Ct") %>%
  dplyr::select(Run, Replicate, BSH_Enzyme, log2_concentration) %>%
  rename(log2_conc_NI = log2_concentration)

fitted_TUDCA_plot = fitted_TUDCA %>%
  left_join(ref_TUDCA, by = c("Run", "Replicate", "BSH_Enzyme")) %>%
  filter(!is.na(log2_conc_NI))

library(ggplot2)
observed_df <- fitted_TUDCA_plot |> 
  mutate(PointType = "Observed")

fitted_df <- fitted_TUDCA_plot |>
  mutate(
    log2_concentration = fitted_values_9,
    PointType = "Fitted"
  )

plot_df <- rbind(observed_df, fitted_df)

ggplot(plot_df, aes(x = log2_conc_NI, y = log2_concentration)) +
  geom_abline(
    intercept = 0, slope = 1,
    linetype = "dotted",
    color = "black",
    size = 1,
    show.legend = TRUE
  ) +
  geom_smooth(
    data = fitted_df,
    method = "lm",
    aes(linetype = "Regression (Fitted)", color = "Regression (Fitted)"),
    se = FALSE,
    size = 1
  ) +
  geom_point(
    aes(color = PointType),
    size = 2,
    alpha = 0.8
  ) +
  scale_color_manual(
    name = "Data",
    values = c(
      "Observed" = "black",
      "Fitted" = "blue",
      "Regression (Fitted)" = "blue"
    )
  )+ 
  
  labs(
    x = expression(Log[2]*" TUDCA Concentration (No Inhibitor – NI)"),
    y = expression(Log[2]*" TUDCA Concentration (With Inhibitor)"),
    title = "Observed vs Fitted TUDCA Values with Identity & Model Regression Lines"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )

############CONFIDENCE INTERVAL#############

#########CI for the fitted values###########
library(merTools)

###for fitted data
fitted_ci <- predictInterval(model_9, 
                             newdata = df_combined,
                             level = 0.95,
                             n.sims = 200,include.resid.var = FALSE)
pred_ci <- predictInterval(model_9, include.resid.var = TRUE)

mean(fitted_ci$upr - fitted_ci$lwr)  # Narrower
mean(pred_ci$upr - pred_ci$lwr)      # Wider


###############################PLOT THE CI FOR TUDCA########################
library(merTools)
library(ggplot2)
library(dplyr)

tudca_data = df_combined %>%
  filter(Bile_Acid == "TUDCA", Inhibitor != "VI", BSH_Enzyme != "Ct")

tudca_ci = predictInterval(
  model_9,
  newdata = tudca_data,
  n.sims = 200,
  include.resid.var = FALSE
) %>%
  bind_cols(tudca_data)

##FOR EACH INHIBITOR
tudca_summary <- tudca_ci %>%
  group_by(Inhibitor) %>%
  summarise(
    avg_fit = mean(fit),
    avg_upr = mean(upr),
    avg_lwr = mean(lwr),
    n = n()
  ) %>%
  ungroup()

#############CI FOR EACH INHIBITOR#############
print(tudca_summary,n=20)

###INHIBITOR VERSUS NO INHIBITOR###########
tudca_summary <- tudca_ci %>%
  mutate(Inhibitor_Status = ifelse(Inhibitor == "NI", "No Inhibitor", "With Inhibitor")) %>%
  group_by(Inhibitor_Status) %>%
  summarise(
    avg_fit = mean(fit),
    avg_upr = mean(upr),
    avg_lwr = mean(lwr)
  )

print(tudca_summary)
with_inhib_avg <- tudca_summary %>% 
  filter(Inhibitor_Status == "With Inhibitor")

ggplot(plot_df, aes(x = log2_conc_NI, y = log2_concentration)) +
  # Reference line (y=x)
  geom_abline(
    intercept = 0, slope = 1,
    linetype = "dotted",
    color = "black",
    linewidth = 1
  ) +
  
  # Regression line
  geom_smooth(
    data = fitted_df,
    method = "lm",
    aes(linetype = "Regression (Fitted)", color = "Regression (Fitted)"),
    se = FALSE,
    linewidth = 1
  ) +
  
  # Points
  geom_point(
    aes(color = PointType),
    size = 2,
    alpha = 0.8
  ) +
  
  # CI lines for With Inhibitor
  geom_hline(
    data = with_inhib_avg,
    aes(yintercept = avg_fit, linetype = "Avg (With Inhibitor)"),
    color = "red3",
    linewidth = 1
  ) +
  geom_hline(
    data = with_inhib_avg,
    aes(yintercept = avg_upr, linetype = "95% CI"),
    color = "red3",
    linewidth = 0.6
  ) +
  geom_hline(
    data = with_inhib_avg,
    aes(yintercept = avg_lwr, linetype = "95% CI"),
    color = "red3",
    linewidth = 0.6
  )+
  
  # Text label
  geom_text(
    data = with_inhib_avg,
    aes(
      x = max(plot_df$log2_conc_NI),
      y = avg_fit,
      label = paste0("Avg = ", round(avg_fit, 2))
    ),
    color = "red3",
    hjust = 1.1,
    vjust = -0.5,
    size = 3
  )+
  
  # Color and line type scales
  scale_color_manual(
    name = "Data",
    values = c(
      "Observed" = "black",
      "Fitted" = "blue",
      "Regression (Fitted)" = "blue"
    )
  ) +
  scale_linetype_manual(
    name = "Reference Lines",
    values = c(
      "Regression (Fitted)" = "solid",
      "Avg (With Inhibitor)" = "solid",
      "95% CI" = "dashed"
    ),
    guide = guide_legend(
      override.aes = list(
        color = c("red3", "red3", "blue"),
        linewidth = c(1, 1, 0.6)
      )
    )
  ) +
  scale_y_continuous(
    name = expression(Log[2]*" TUDCA Concentration (With Inhibitor)"),
    breaks = c(
      pretty(range(plot_df$log2_concentration), n = 5),
      round(with_inhib_avg$avg_fit,3),
      round(with_inhib_avg$avg_lwr,3),
      round(with_inhib_avg$avg_upr,3)
    ),
    labels = function(x) {
      ifelse(
        x %in% round(c(with_inhib_avg$avg_fit, with_inhib_avg$avg_lwr, with_inhib_avg$avg_upr), 3),
        paste0(x, "  ◀"),
        x
      )
    }
  )+
  
  labs(
    x = expression(Log[2]*" TUDCA Concentration (No Inhibitor – NI)"),
    y = expression(Log[2]*" TUDCA Concentration (With Inhibitor)"),
    title = "Observed vs Fitted TUDCA Values with Identity & Model Regression Lines"
  )+ 
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(0.2, "cm")
  )


###################BETTER CI PLOT FOR ENZYME AND BILE ACID##########################
library(dplyr)

fitted_ci_df = bind_cols(df_combined, fitted_ci)
summary_df <- fitted_ci_df %>%
  group_by(BSH_Enzyme, Inhibitor, Bile_Acid) %>%
  summarise(
    mean_fit = mean(fit),
    lwr = mean(lwr),
    upr = mean(upr),
    .groups = "drop"
  )

library(ggplot2)

plot_by_enzyme <- function(enzyme_name) {
  enzyme_df <- summary_df %>% filter(BSH_Enzyme == enzyme_name)
  
  ggplot(enzyme_df, aes(x = Bile_Acid, y = mean_fit, fill = Bile_Acid)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = lwr, ymax = upr),
                  position = position_dodge(width = 0.8), width = 0.2) +
    facet_wrap(~ Inhibitor) +
    labs(
      title = paste("Estimated BSH Activity for Enzyme:", enzyme_name),
      y = expression(Log[2]*" Fitted Concentration"),
      x = "Bile Acid"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

for (e in unique(summary_df$BSH_Enzyme)) {
  print(plot_by_enzyme(e))
}

#############overall PLOT for CI for the entire dataset#################
library(dplyr)
library(merTools)
library(ggplot2)

# Reference (NI): exclude Ct
ref_NI <- df_combined %>%
  filter(Inhibitor == "NI", BSH_Enzyme != "Ct") %>%
  dplyr::select(Run, Replicate, BSH_Enzyme, Bile_Acid, log2_concentration) %>%
  rename(log2_conc_NI = log2_concentration)
with_inhib <- df_combined %>%
  filter(Inhibitor != "NI", Inhibitor != "VI", BSH_Enzyme != "Ct") %>%
  dplyr::select(Run, Replicate, BSH_Enzyme, Bile_Acid, Inhibitor, log2_concentration)

# Join reference + with-inhibitor observed
plot_df <- inner_join(with_inhib, ref_NI,
                      by = c("Run", "Replicate", "BSH_Enzyme", "Bile_Acid")) %>%
  mutate(PointType = "Observed")

fitted_ci <- predictInterval(
  model_9, 
  newdata = df_combined,
  n.sims = 200,
  include.resid.var = FALSE
)

df_combined_pred <- bind_cols(df_combined, fitted_ci)

# Fitted data: with inhibitor only
fitted_df <- df_combined_pred %>%
  filter(Inhibitor != "NI", Inhibitor != "VI", BSH_Enzyme != "Ct") %>%
  dplyr::select(Run, Replicate, BSH_Enzyme, Bile_Acid, Inhibitor, fit) %>%
  inner_join(ref_NI, by = c("Run", "Replicate", "BSH_Enzyme", "Bile_Acid")) %>%
  mutate(log2_concentration = fit, PointType = "Fitted")
plot_all <- bind_rows(plot_df, fitted_df)

with_inhib_avg <- df_combined_pred %>%
  filter(Inhibitor != "NI", Inhibitor != "VI", BSH_Enzyme != "Ct") %>%
  summarise(
    avg_fit = mean(fit),
    avg_lwr = mean(lwr),
    avg_upr = mean(upr)
  )

ggplot(plot_all, aes(x = log2_conc_NI, y = log2_concentration)) +
  # Identity line
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black", linewidth = 1) +
  
  # Regression line from fitted values
  geom_smooth(
    data = filter(plot_all, PointType == "Fitted"),
    method = "lm",
    aes(linetype = "Regression (Fitted)", color = "Regression (Fitted)"),
    se = FALSE,
    linewidth = 1
  ) +
  
  # Points
  geom_point(
    aes(color = PointType),
    size = 2,
    alpha = 0.6
  ) +
  
  # CI lines for With Inhibitor
  geom_hline(data = with_inhib_avg, aes(yintercept = avg_fit, linetype = "Avg (With Inhibitor)"),
             color = "red3", linewidth = 1) +
  geom_hline(data = with_inhib_avg, aes(yintercept = avg_upr, linetype = "95% CI"),
             color = "red3", linewidth = 0.6) +
  geom_hline(data = with_inhib_avg, aes(yintercept = avg_lwr, linetype = "95% CI"),
             color = "red3", linewidth = 0.6) +
  
  # Label the mean
  geom_text(data = with_inhib_avg,
            aes(x = max(plot_df$log2_conc_NI), y = avg_fit),
            label = "Mean (With Inhibitor)", color = "red3",
            hjust = 1.1, vjust = -0.5, size = 3) +
  
  # Custom scales
  scale_color_manual(
    name = "Data",
    values = c("Observed" = "black", "Fitted" = "blue", "Regression (Fitted)" = "blue")
  ) +
  scale_linetype_manual(
    name = "Reference Lines",
    values = c("Regression (Fitted)" = "solid", "Avg (With Inhibitor)" = "solid", "95% CI" = "dashed"),
    guide = guide_legend(
      override.aes = list(color = c("red3", "red3", "blue"), linewidth = c(1, 0.6, 1))
    )
  ) +
  scale_y_continuous(
    name = expression(Log[2]*" Concentration (With Inhibitor)"),
    breaks = c(
      pretty(range(plot_all$log2_concentration), n = 5),
      round(with_inhib_avg$avg_fit, 3),
      round(with_inhib_avg$avg_lwr, 3),
      round(with_inhib_avg$avg_upr, 3)
    ),
    labels = function(x) {
      ifelse(
        x %in% round(c(with_inhib_avg$avg_fit, with_inhib_avg$avg_lwr, with_inhib_avg$avg_upr), 3),
        paste0(x, "  ◀"),
        x
      )
    }
  ) +
  labs(
    x = expression(Log[2]*" Concentration (No Inhibitor – NI)"),
    y = expression(Log[2]*" Concentration (With Inhibitor)"),
    title = "Observed vs Fitted Values for All Bile Acids and Enzymes"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(0.2, "cm")
  )+
  annotate("text",
           x = max(plot_all$log2_conc_NI),
           y = min(plot_all$log2_concentration),
           label = "Interpretation:\nAbove 0 → Inhibition\nAt 0 → No Effect\nBelow 0 → No Inhibition",
           hjust = 1.1, vjust = -0.5,
           size = 4,
           fontface = "italic",
           color = "gray30")



#######plot for TUDCA#######
ggplot(fitted_TUDCA_plot, aes(x = log2_conc_NI, y = log2_concentration)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray30") +
  geom_point(aes(color = BSH_Enzyme, shape = Replicate), size = 3, alpha = 0.85) +
  geom_point(aes(y = fitted_values_9), shape = 18, color = "black", size = 3) +
  facet_wrap(~Inhibitor, scales = "free") +
  scale_x_continuous(
    name = expression(Log[2]*" Ratio (Ct+NI / Enzyme+NI) TUDCA Concentration (No Inhibitor - NI)")
  ) +
  scale_y_continuous(
    name = expression(Log[2]*" Ratio (Enzyme+Inhibitor / Enzyme+NI) TUDCA Concentration (With Inhibitor)")
  ) +
  scale_color_manual(values = enzyme_colors) +
  labs(
    title = expression(Log[2]~"TUDCA Ratios: Effect of Inhibitors on BSH Enzyme Activity")
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(angle = 45),
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray95", color = NA)
  ) +
  guides(
    color = guide_legend(nrow = 1, override.aes = list(size = 4)),
    shape = guide_legend(title = "Replicate", override.aes = list(size = 4))
  )



########better colours#####################
ggplot(fitted_TUDCA_plot,
       aes(x = log2_conc_NI, y = log2_concentration)) +
  
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray30") +
  
  
  geom_point(
    aes(colour = BSH_Enzyme, shape = Replicate),
    size  = 3,
    alpha = 0.7
  ) +
  
  geom_point(
    aes(y = fitted_values_9, shape = Replicate),
    colour = "black",
    size   = 3,      # a bit bigger
    stroke = 1.1,
    alpha  = 1,
    show.legend = FALSE   # avoid duplicate legend
  ) +
  
  geom_point(
    aes(y = fitted_values_9, colour = BSH_Enzyme, shape = Replicate),
    size   = 2,
    alpha  = 1,
    show.legend = FALSE
  ) +
  
  facet_wrap(~Inhibitor, scales = "free") +
  scale_colour_manual(values = enzyme_colors) +
  labs(
    x = expression(Log[2]*" Ratio (Ct+NI / Enzyme+NI) TUDCA Concentration (No Inhibitor – NI)"),
    y = expression(Log[2]*" Ratio (Enzyme+Inhibitor / Enzyme+NI) TUDCA Concentration (With Inhibitor)"),
    title = expression(Log[2]~"TUDCA Ratios: Effect of Inhibitors on BSH Enzyme Activity")
  )


############EXTRACT COEFFICIENTS VALUES###############
coef_df <- data.frame(
  term = names(fixef(model_9)),
  estimate = as.numeric(fixef(model_9))
)
coef_df

str(coef_df)

library(lme4)
library(broom.mixed)

set.seed(6090)
K = 5
n = nrow(df_combined)
folds = sample(rep(1:K, length.out = n))
cv_coefs = list()

for (k in 1:K) {
  i.train = which(folds != k)
  dat.train = df_combined[i.train, ]
  dat.test  = df_combined[-i.train, ]
  
  model_k = lmer(log2_concentration ~ 
                   BSH_Enzyme * Inhibitor * Bile_Acid + 
                   (1|Run) + (1|Run:Replicate),
                 data = dat.train)
  
  coefs_k = broom.mixed::tidy(model_k, effects = "fixed")
  coefs_k$Fold = k
  cv_coefs[[k]] = coefs_k
}
cv_results = do.call(rbind, cv_coefs)

library(dplyr)

cv_summary <- cv_results %>%
  group_by(term) %>%
  summarise(mean_estimate = mean(estimate), .groups = "drop")
##########extract TUDCA value#########
library(emmeans)

emm <- emmeans(model_9, ~ BSH_Enzyme * Inhibitor * Bile_Acid,pbkrtest.limit = 5000)
tudca_emm <- emmeans(
  model_9,
  ~ BSH_Enzyme * Inhibitor * Bile_Acid,
  at = list(Bile_Acid = "TUDCA"),
  mode = "latent",            # speeds up prediction
  lmer.df = "asymptotic"      # avoids slow df calc
) %>%
  as.data.frame()

tudca_df <- as.data.frame(tudca_emm)

tudca_df_clean <- tudca_df %>%
  filter(!is.na(emmean))  
tudca_final <- tudca_df %>%
  filter(!is.na(emmean)) %>%
  mutate(term = ifelse(Inhibitor == "NI",
                       paste0("BSH_Enzyme", BSH_Enzyme,
                              ":Bile_Acid", Bile_Acid),
                       paste0("BSH_Enzyme", BSH_Enzyme,
                              ":Inhibitor", Inhibitor,
                              ":Bile_Acid", Bile_Acid))) %>%
  dplyr::select(term, mean_estimate = emmean)

head(tudca_final)

combined_summary <- bind_rows(cv_summary, tudca_final)
intercept_value <- combined_summary$mean_estimate[combined_summary$term == "(Intercept)"]

get_coef <- function(term, coef_table) {
  if (is.na(term)) return(0)
  val <- coef_table$mean_estimate[coef_table$term == term]
  if (length(val) == 0) return(0)
  return(val)
}

df_combined <- df_combined %>%
  mutate(
    enzyme_term = if_else(BSH_Enzyme != "Ct", paste0("BSH_Enzyme", BSH_Enzyme), NA_character_),
    inhibitor_term = if_else(Inhibitor != "NI", paste0("Inhibitor", Inhibitor), NA_character_),
    bileacid_term = paste0("Bile_Acid", Bile_Acid),
    
    #Interaction term names
    ei_term = if_else(!is.na(enzyme_term) & !is.na(inhibitor_term),
                      paste0(enzyme_term, ":", inhibitor_term), NA_character_),
    
    eb_term = if_else(!is.na(enzyme_term),
                      paste0(enzyme_term, ":", bileacid_term), NA_character_),
    
    eib_term = if_else(!is.na(enzyme_term) & !is.na(inhibitor_term),
                       paste0(enzyme_term, ":", inhibitor_term, ":", bileacid_term), NA_character_)
  ) %>%
  rowwise() %>%
  mutate(
    enzyme_baseline_9 = if_else(BSH_Enzyme == "Ct", 0,
                                get_coef(enzyme_term, combined_summary)),
    
    #inhibitor_baseline_9 = if_else(Inhibitor == "NI", 0,
                                   #get_coef(paste0("Inhibitor", Inhibitor), combined_summary)),
    
    #include enzyme:bileacid interaction (unless BSH_Enzyme is Ct)
    enzyme_bileacid_interaction_9 = if_else(BSH_Enzyme == "Ct", 0,
                                            get_coef(eb_term, combined_summary)),
    
    #if Inhibitor is not NI
    enzyme_inhibitor_interaction_9 = if_else(Inhibitor == "NI", 0,
                                             get_coef(ei_term, combined_summary)),
    
    enzyme_inhibitor_bileacid_interaction_9 = if_else(Inhibitor == "NI", 0,
                                                      get_coef(eib_term, combined_summary)),
    
    inhibitor_potency_9 = if_else(Inhibitor != "NI",
                                    (enzyme_baseline_9 + enzyme_bileacid_interaction_9 + enzyme_inhibitor_bileacid_interaction_9)-(enzyme_baseline_9 + enzyme_bileacid_interaction_9),
                                  0) ) %>%
  ungroup() %>%
  dplyr::select(-enzyme_term, -inhibitor_term, -bileacid_term, -ei_term, -eb_term, -eib_term)


nrow(df_combined)  # Should be 4950


write.csv(df_combined, "/Users/keerthim/Documents/ST4060/ST6090/R CODE/df_combined_with_components.csv", row.names = FALSE)
write.csv(combined_summary, "/Users/keerthim/Documents/ST4060/ST6090/R CODE/cv_summary_coefficients.csv", row.names = FALSE)


################INSIGHTS################################

library(ggplot2)
library(dplyr)


df_combined %>%
  filter(Inhibitor == "NI", BSH_Enzyme != "Ct") %>%
  group_by(BSH_Enzyme, Bile_Acid) %>%
  summarise(mean_fitted = mean(fitted_values_9, na.rm = TRUE), .groups = "drop") %>%
  group_by(Bile_Acid) %>%
  mutate(highlight = BSH_Enzyme == BSH_Enzyme[which.max(mean_fitted)]) %>%
  ggplot(aes(x = reorder(BSH_Enzyme, -mean_fitted), y = mean_fitted, fill = highlight)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "steelblue"), guide = "none") +
  labs(
    title = "Top BSH Enzyme per Bile Acid (Inhibitor = NI)",
    x = "BSH Enzyme", y = "Mean log₂ concentration"
  ) +
  facet_wrap(~ Bile_Acid, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#############INHIBITOR POTENCY PLOT#########
library(ggplot2)
library(dplyr)


df_combined %>%
  filter(Inhibitor != "NI") %>%  # Exclude reference inhibitor
  group_by(BSH_Enzyme, Inhibitor, Bile_Acid) %>%
  summarise(mean_potency = mean(inhibitor_potency_9, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = Inhibitor, y = BSH_Enzyme, fill = mean_potency)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue", mid = "#FFE5E5", high = "red", midpoint = 0,
    name = "Potency"
  ) +
  facet_wrap(~ Bile_Acid, ncol = 3) +
  labs(
    title = "Inhibitor Potency per Enzyme Across Bile Acids (Excluding NI)",
    x = "Inhibitor",
    y = "BSH Enzyme"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )





library(dplyr)
library(ggplot2)

df_combined %>%
  filter(Inhibitor != "NI") %>%  # Exclude the reference inhibitor
  group_by(BSH_Enzyme, Inhibitor, Bile_Acid) %>%
  summarise(mean_potency = mean(inhibitor_potency_9, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = reorder(BSH_Enzyme, -mean_potency), y = mean_potency)) +
  geom_col(fill = "darkorange") +
  facet_grid(Bile_Acid ~ Inhibitor, scales = "free_y") +
  labs(
    title = "Inhibitor Potency per Enzyme, Faceted by Bile Acid and Inhibitor",
    x = "BSH Enzyme",
    y = "Mean Inhibitor Potency"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 14)
  )

library(dplyr)
library(ggplot2)

df_combined %>%
  filter(BSH_Enzyme != "Ct") %>%  # ← Remove control enzyme
  group_by(BSH_Enzyme, Inhibitor, Bile_Acid) %>%
  summarise(mean_concentration = mean(fitted_values_9, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = reorder(BSH_Enzyme, -mean_concentration), y = mean_concentration)) +
  geom_col(fill = "steelblue") +
  facet_grid(Bile_Acid ~ Inhibitor, scales = "free_y") +
  labs(
    title = "Mean log2 Bile Acid Concentration per Enzyme\nFaceted by Inhibitor and Bile Acid",
    x = "BSH Enzyme",
    y = "Mean log2 Concentration"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 14)
  )


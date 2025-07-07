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

d=d[d$Replicate=="Average",]
d$Replicate=droplevels(d$Replicate)
str(d)



##################correlation####################
num_vars <- names(d)[sapply(d, is.numeric)]
exclude_vars <- c("BSH_IP_Conc", "Inhibitor_IP_Conc")
num_vars_filtered <- setdiff(num_vars, exclude_vars)

numeric_data_filtered <- d[, num_vars_filtered]


cor_matrix_filtered <- cor(numeric_data_filtered, use = "pairwise.complete.obs")


print(cor_matrix_filtered)


high_cor <- which(abs(cor_matrix_filtered) > 0.7 & abs(cor_matrix_filtered) < 1, arr.ind = TRUE)
# Keep only unique pairs (i < j)
high_cor <- high_cor[high_cor[,1] < high_cor[,2],]

cat("Highly correlated variable pairs (|correlation| > 0.7), excluding BSH_IP_Conc and Inhibitor_IP_Conc:\n")
apply(high_cor, 1, function(idx) {
  cat(rownames(cor_matrix_filtered)[idx[1]], "-", colnames(cor_matrix_filtered)[idx[2]], 
      ": ", round(cor_matrix_filtered[idx[1], idx[2]], 3), "\n")
})


dev.new()
library(corrplot)
corrplot(cor_matrix_filtered, method = "color", type = "full", 
         tl.cex = 0.7, tl.col = "black", diag = FALSE)

##############MULTIVARIATE RANDOM FOREST #############
if (!require(randomForestSRC)) {
  install.packages("randomForestSRC")
  library(randomForestSRC)
} else {
  library(randomForestSRC)
}

bile_acids <- c("TUDCA", "TCDCA", "TDCA", "TCA", "TLCA",
                "GUDCA", "GCDCA", "GDCA", "GCA", "GLCA")

inhibitor_vars <- c("inhibitor_compound_canonicalized", "Inhibitor_Complexity",
                    "Inhibitor_H_Bond_Acceptor", "Inhibitor_H_Bond_Donor",
                    "Inhibitor_Rotatable_Bonds", "Inhibitor_Molecular_Weight",
                    "Inhibitor_Polar_Surface_Area")

protein_vars <- grep("^Protein_", names(d), value = TRUE)

predictor_vars <- c(protein_vars, inhibitor_vars)


model_data <- d[, c(predictor_vars, bile_acids)]
model_data <- na.omit(model_data)

#multivariate random forest model
rf_model <- rfsrc(Multivar(TUDCA, TCDCA, TDCA, TCA, TLCA, 
                           GUDCA, GCDCA, GDCA, GCA, GLCA) ~ ., 
                  data = model_data, 
                  importance = TRUE)


print(rf_model)



###VARIABLE IMPORTANCE

importance_list <- lapply(bile_acids, function(ba) {
  rf_model$regrOutput[[ba]]$importance
})
imp_matrix <- do.call(cbind, importance_list)
colnames(imp_matrix) <- bile_acids

#MEAN IMPORTANCE
mean_imp <- rowMeans(imp_matrix, na.rm = TRUE)

##SORT THE DATA
imp_df <- data.frame(Feature = names(mean_imp),
                     MeanImportance = mean_imp)
imp_df <- imp_df[order(-imp_df$MeanImportance), ]


library(ggplot2)

ggplot(head(imp_df, 20), aes(x = reorder(Feature, MeanImportance), y = MeanImportance)) +
  geom_col(fill = "darkgreen") +
  coord_flip() +
  labs(title = "Top 20 Important Features (Mean VIMP across Bile Acids)",
       x = "Feature", y = "Mean Variable Importance") +
  theme_minimal()



###TOP 10 FOR EACH BA#####
top_n <- 10 

top_features_per_ba <- lapply(bile_acids, function(ba) {
  importance_vec <- rf_model$regrOutput[[ba]]$importance
  sorted <- sort(importance_vec, decreasing = TRUE)
  head(sorted, top_n)
})


names(top_features_per_ba) <- bile_acids


for (ba in bile_acids) {
  cat("\nTop features for", ba, ":\n")
  print(top_features_per_ba[[ba]])
}

library(ggplot2)

ba <- "TUDCA"
df <- data.frame(Feature = names(top_features_per_ba[[ba]]),
                 Importance = top_features_per_ba[[ba]])

ggplot(df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = paste("Top Features for", ba),
       x = "Feature", y = "Variable Importance") +
  theme_minimal()


#################CHECK THE relation between concentration to h_bond_acceptor and h_bond_donor


bile_acids <- c("TUDCA", "TCDCA", "TDCA", "TCA", "TLCA",
                "GUDCA", "GCDCA", "GDCA", "GCA", "GLCA")

control_means <- d %>%
  filter(BSH_Enzyme == "JL885", Inhibitor == "NI") %>%
  summarise(across(all_of(bile_acids), \(x) mean(x, na.rm = TRUE))) %>%
  pivot_longer(cols = everything(), names_to = "Bile_Acid", values_to = "Control_Mean")

ggplot(df_long, aes(x = Inhibitor_H_Bond_Acceptor, y = Concentration, color = Inhibitor)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_hline(data = control_means, aes(yintercept = Control_Mean), 
             linetype = "dotted", color = "red", size = 0.7) +
  facet_wrap(~ Bile_Acid, scales = "free_y") +
  labs(title = "JL885: Bile Acid Concentration vs Inhibitor H Bond Acceptor",
       x = "Inhibitor H Bond Acceptor",
       y = "Bile Acid Concentration") +
  theme_minimal() +
  theme(legend.position = "bottom")




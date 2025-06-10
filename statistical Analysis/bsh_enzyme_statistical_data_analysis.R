
rm(list=ls())
d = read.csv("/Users/keerthim/Documents/ST4060/ST6090/CLUSTERING DATA/BSH_compound_protein_bile_acid_detailed.csv",header=TRUE)

head(d)
dim(d)
str(d)

d$Replicate= NULL
d$Sample_ID=NULL
d$Run=NULL
d$Inhibitor=as.factor(d$Inhibitor)
d$BSH_Enzyme=as.factor(d$BSH_Enzyme)
d$Protein_Organism=NULL

head(d)

d <- d[
  !(d$BSH_Enzyme == "Ct" & 
      d$Inhibitor == "NI"), ]

head(d)
dim(d)
str(d)

d$BSH_Enzyme = droplevels(d$BSH_Enzyme)
d$Inhibitor = droplevels(d$Inhibitor)


numeric_cols <- sapply(d, is.numeric)
d_numeric <- d[, numeric_cols]

pairs(d_numeric[, 1:6], pch = 20)

###INITIAL ANALYSES

d_enzyme = subset(d, Inhibitor == "NI")


bile_acids = c("TUDCA", "TCDCA", "TDCA", "TCA", "TLCA", 
                "GUDCA", "GCDCA", "GDCA", "GCA", "GLCA")


library(reshape2)
df_long <- melt(d_enzyme[, c("BSH_Enzyme", bile_acids)], 
                id.vars = "BSH_Enzyme", 
                variable.name = "Bile_Acid", 
                value.name = "Concentration")

# Aggregate mean bile acid concentration per enzyme per bile acid
library(dplyr)
df_summary <- df_long %>%
  group_by(BSH_Enzyme, Bile_Acid) %>%
  summarise(Mean_Concentration = mean(Concentration, na.rm = TRUE))

# Plotting
library(ggplot2)
ggplot(df_summary, aes(x = Bile_Acid, y = Mean_Concentration, fill = BSH_Enzyme)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Enzyme Activity on Bile Acids",
       x = "Bile Acid",
       y = "Mean conjugated BA Concentration (µg/ml)",
       fill = "BSH Enzyme") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Define new labels for plotting only
enzyme_labels <- c(
  "JL885" = "BSH-JL885 (positive control)",
  "A"     = "BSH-A",
  "B"     = "BSH-B",
  "F"     = "BSH-F",
  "1011c" = "BSH-1011c",
  "T2"    = "BSH-T2"
)


###########enzyme BA activity###############
# Plot with relabeled enzymes
ggplot(df_summary, aes(x = Bile_Acid, y = Mean_Concentration, fill = BSH_Enzyme)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "BSH Activity on Bile Acids(No inhibitor)",
       x = "Bile Acid", y = "Mean conjugated BA Concentration (µg/ml)", fill = "BSH Enzyme") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(6, "Set2"),
                    labels = enzyme_labels) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

df_summary[df_summary$Mean_Concentration > 5, ]


############Enzyme Inhibitor Activity on BA###############


# Filter out rows with Inhibitor == "NI"
df_inhibitor = subset(d, Inhibitor != "NI")

df_inhibitor = df_inhibitor[, c("BSH_Enzyme", "Inhibitor", bile_acids)]

# Computing the mean concentration for each Enzyme-Inhibitor pair
df_summary = aggregate(. ~ BSH_Enzyme + Inhibitor, data = df_inhibitor, FUN = mean, na.rm = TRUE)

# Melt the summary for plotting
df_melted = melt(df_summary, id.vars = c("BSH_Enzyme", "Inhibitor"), 
                  variable.name = "Bile_Acid", value.name = "Mean_Concentration")


ggplot(df_melted, aes(x = BSH_Enzyme, y = Mean_Concentration, fill = Inhibitor)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ Bile_Acid, scales = "free_y") +
  labs(title = "Effect of Inhibitors on Enzyme Activity Across Bile Acids",
       x = "BSH Enzyme", y = "Mean BA Concentration((µg/ml)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#####for each BA and enzyme################

  
doc_path = "/Users/keerthim/Documents/ST4060/ST6090/plots"

if (!dir.exists(doc_path)) {
  dir.create(doc_path, recursive = TRUE)
}

for (enzyme in unique(df_melted$BSH_Enzyme)) {
  df_enzyme = subset(df_melted, BSH_Enzyme == enzyme)
  
  for (ba in unique(df_enzyme$Bile_Acid)) {
    df_plot = subset(df_enzyme, Bile_Acid == ba)
    
    p = ggplot(df_plot, aes(x = Inhibitor, y = Mean_Concentration, fill = Inhibitor)) +
      geom_bar(stat = "identity") +
      labs(title = paste("Effect of Inhibitors on", enzyme, "-", ba),
           x = "Inhibitor", y = "Mean BA Concentration (µg/ml)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
    
    file_name = paste0("Effect_", gsub("[^A-Za-z0-9]", "_", enzyme), "_", ba, ".pdf")
    
    ggsave(filename = file.path(doc_path, file_name), plot = p, width = 8, height = 5)
  }
}


#########top 5 inhibitor for each bile acid and enzyme
library(dplyr)

top5_inhibitors_table = df_melted %>%
  group_by(BSH_Enzyme, Bile_Acid) %>%
  arrange(desc(Mean_Concentration)) %>%
  slice_head(n = 1) %>%
  ungroup()

print(top5_inhibitors_table)

########visualize it


ggplot(top5_inhibitors_table, aes(x = BSH_Enzyme, y = Inhibitor, fill = Mean_Concentration)) +
  geom_tile(color = "white") +
  facet_wrap(~ Bile_Acid, scales = "free") +
  scale_fill_gradientn(
    colours = c("#33a02c", "#ffffbf", "#e31a1c"),  # Green -> Yellow -> Red
    name = "BA Conc (µg/ml)"
  ) +
  labs(title = "Top 5 Inhibitors per Enzyme and Bile Acid (Heatmap)",
       x = "BSH Enzyme", y = "Inhibitor") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

####FIND A COMMON INHIBITOR FOR ENZYME FOR BOTH TAURINE AND GLYCINE HEADED

library(dplyr)
library(stringr)

# Taurine-headed BA
taurine_top_inhibitors = df_melted %>%
  filter(str_starts(Bile_Acid, "T")) %>%
  group_by(BSH_Enzyme) %>%
  slice_max(order_by = Mean_Concentration, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Bile_Acid_Headed = "Taurine")

# Glycine-headed BA
glycine_top_inhibitors = df_melted %>%
  filter(str_starts(Bile_Acid, "G")) %>%
  group_by(BSH_Enzyme) %>%
  slice_max(order_by = Mean_Concentration, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Bile_Acid_Headed = "Glycine")

# Combine results
top_inhibitors_by_conc = bind_rows(taurine_top_inhibitors, glycine_top_inhibitors)

top_inhibitors_by_conc

#########################Correlation####################

df = subset(d, select = -c(Inhibitor, BSH_Enzyme))

cor(df)

df =df[, apply(df, 2, sd, na.rm = TRUE) != 0]

cor_matrix =cor(df, use = "pairwise.complete.obs")

diag(cor_matrix)=NA
print(cor_matrix)


high_corr = which(abs(cor_matrix) > 0.7, arr.ind = TRUE)

result = data.frame(
  Var1 = rownames(cor_matrix)[high_corr[,1]],
  Var2 = colnames(cor_matrix)[high_corr[,2]],
  Correlation = cor_matrix[high_corr]
)

result = result[result$Var1 < result$Var2, ]

print(result)

############PCA to find the important features##############
df_pca = d[, sapply(d, is.numeric)]

#columns with zero variance
df_pca  = df_pca[, apply(df_pca, 2, sd, na.rm = TRUE) != 0]

pca = prcomp(df_pca)
pca.s = prcomp(df_pca,scale=TRUE)

par(mfrow=c(2,1))
plot(pca, main="Proportion of variance explained unscaled", 
     xlab="Principle Component")

plot(pca.s, main="Proportion of variance explained scaled ", 
     xlab="Principle Component", col='pink')


# Proportion of variance explained / of information captured
names(pca)
summary(pca)
# reconstruct each of these lines:
pca$sdev
pca$sdev^2/sum(pca$sdev^2)
cumsum(pca$sdev^2/sum(pca$sdev^2))

# explore loadings:
print(pca)
pca$rotation # same matrix :)
apply(abs(pca$rotation),2,which.max)
row.names(pca$rotation)[apply(abs(pca$rotation),2,which.max)]
colnames(pca$rotation)[apply(abs(pca$rotation),2,which.max)]

#scaled data
summary(pca.s)
apply(abs(pca.s$rotation),2,which.max)
row.names(pca.s$rotation)[apply(abs(pca.s$rotation),2,which.max)]

# plot the data on its first 2 dimensions in each space: 
par(mfrow=c(4,1))
plot(df_pca[,3:6], cex=2, pch=20, col="blue", main="Data in original space") 
biplot(pca, main="Data in PCA space") ####plotting the data in PC1, PC2 SPACE
###and it is moving along the first principle component,PCA analysis the structure of correlations in data
biplot(pca.s, main="Data in PCA space")
abline(h=0, col='orange')
abline(v=0, col='orange')


##simply build a classifier using the sign of the coordinates along the PC1
par(mfrow=c(2,2), mar=c(3,3,3,3))
plot(pca) # scree plot 
biplot(pca) # biplot 
plot(pca.s) # scree plot 
biplot(pca.s) # biplot 
# now analyse this plot :)
par(mfrow=c(1,2))
biplot(pca) # biplot 
biplot(pca.s) # biplot 




install.packages("factoextra")
library(ggplot2)
library(factoextra)
library(FactoMineR)

# ----------------------------
# 1. Simulate Data (Replace with YOUR data)
# ----------------------------
# Rows = Samples (e.g., enzyme mutants), Columns = Variables (e.g., bile salt hydrolysis rates)
set.seed(123)


library(ggplot2)
library(factoextra)

# Assuming you have a PCA object named 'pca'
biplot <- fviz_pca_biplot(pca, 
                          col.ind = "black", 
                          pointshape = 19, 
                          pointsize = 3)

# Now apply your theme adjustments
biplot + theme(plot.title = element_text(hjust = 0.5, face = "bold"))


# Get top 10 variables contributing to PC1
contrib <- abs(pca$rotation[, 1:2])
contrib_sum <- rowSums(contrib)
top_contributors <- names(sort(contrib_sum, decreasing = TRUE))[1:5]

fviz_pca_var(pca,
             select.var = list(name = top_contributors),
             col.var = "red",
             repel = TRUE) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

##scaled data

# Get top 10 variables contributing to PC1
contrib.s = abs(pca.s$rotation[, 1:2])
contrib_sum.s = rowSums(contrib.s)
top_contributors.s = names(sort(contrib_sum.s, decreasing = TRUE))[1:2]

fviz_pca_var(pca.s,
             select.var = list(name = top_contributors.s),
             col.var = "red",
             repel = TRUE) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

fviz_pca_biplot(pca)
fviz_pca_biplot(pca.s)



## Script to Understand Side Effects from Yellow Card MHRA with Gene Essentiality and Drug Targets data ##


# Load the necessary libraries ------------------------


library(tidyverse)
library(heatmaply)


# Import the required Datasets ------------------------


#Import SE Dataset

mhra_se_data <- read_csv("C:/Users/HP-ssd/Desktop/Short term project2/mhra/mhra_se_data_with_soc.csv", , 
                         col_types = cols(...1 = col_skip()))

mhra_se_data <- mhra_se_data %>%
  na.omit()

#Import Parlgoy data

human_gene_paralogues <- read_csv("C:/Users/HP-ssd/Desktop/Short term project2/paralogues/human_gene_paralogues.csv", 
                                  col_types = cols(...1 = col_skip()))

#Import FUSIL data

fusil <- read_csv("C:/Users/HP-ssd/Desktop/Short term project2/fusil.csv")


#Import open target data

opentarget_data <- read_csv("C:/Users/HP-ssd/Desktop/Short term project2/Open Target drug data/Opentarget_data.csv", 
                            col_types = cols(...1 = col_skip()))

opentarget_data <- opentarget_data %>%
  replace_na( list(status = 'Unknown status'))



# Filter for Approved drugs only and Integrate FUSIL data----------------------

opentarget_approved <- opentarget_data %>%
  filter(status == "Completed")



# Integrate fusil data

drug_gene <- opentarget_approved %>%
  select(approvedSymbol, prefName) %>%
  unique()

drug_fusil <- fusil %>%
  full_join(drug_gene, by = c( "gene_symbol" = "approvedSymbol")) %>%
  na.omit() %>%
  select(3,4,5) 

# Merge open target and SE data

# filter for drugs in opentarget


mhra_df <- mhra_se_data %>%
  left_join(drug_gene, by = c ("DRUG" = "prefName")) %>%
  na.omit() %>%
  filter(seriousness != "FALSE")


# Merge the FUSIL Drugs Data with Side Effects --------


mhra_fusil <- mhra_se_data %>%
  full_join(drug_fusil, by = c ("DRUG" = "prefName")) %>%
  na.omit() %>%
  filter(seriousness != "FALSE")

length(unique(mhra_fusil$DRUG)) #489 drugs
length(unique(mhra_fusil$gene_symbol)) #308 genes





# count of serious SE for each FUSIL bin --------------

count_fusil_se <- mhra_fusil %>%
  count(fusil)

se_serious <- mhra_fusil %>%
  group_by(fusil, seriousness) %>%
  count(fusil) %>%
  left_join(count_fusil_se, by = "fusil") %>%
  mutate(percentage = (n.x*100)/n.y)


#write.csv(se_serious, "C:/Users/HP-ssd/Desktop/se_serious.csv") # Very Difficult to see trends
se_serious$fusil <- factor(se_serious$fusil, 
                           levels = c("CL", "DL", "SV", "VP", "VnP"  ))

ggplot(se_serious, aes(x = fusil, y = percentage, fill = seriousness)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Percentage of Serious Side Effects in each FUSIL",
       x = "FUSIL bin",
       y = "Percentage") +
  theme_minimal()



# #Compute Enrichment for Each HLGT–FUSIL Pair --------

#SE group for each FUSIL

count_fusil_se <- mhra_fusil %>%
  count(fusil)

se_group <- mhra_fusil %>%
  group_by(fusil, HLGT) %>%
  count(fusil) %>%
  left_join(count_fusil_se, by = "fusil") %>%
  mutate(percentage = (n.x*100)/n.y)


# Step 1: Group by FUSIL and HLGT and sum counts
group_counts <- se_group %>%
  group_by(fusil, HLGT) %>%
  summarise(n_x = sum(n.x), .groups = 'drop')

# Step 2: Total counts per FUSIL
fusil_totals <- group_counts %>%
  group_by(fusil) %>%
  summarise(fusil_total = sum(n_x))

# Step 3: Merge and calculate proportion within FUSIL
group_counts <- group_counts %>%
  left_join(fusil_totals, by = "fusil") %>%
  mutate(fusil_prop = n_x / fusil_total)

# Step 4: Overall HLGT totals and proportions
overall_totals <- se_group %>%
  group_by(HLGT) %>%
  summarise(overall_total = sum(n.x)) %>%
  mutate(overall_prop = overall_total / sum(overall_total))

# Step 5: Merge and compute relative enrichment
enrichment_df <- group_counts %>%
  left_join(overall_totals, by = "HLGT") %>%
  mutate(relative_enrichment = fusil_prop / overall_prop)




# Visualsie outliers

enrichment_qc <- enrichment_df %>%
  mutate(log_enrichment = log2(relative_enrichment))

enrichment_qc$fusil <- factor(enrichment_qc$fusil, 
                              levels = c("CL", "DL", "SV", "VP", "VnP"  ))

ggplot(enrichment_qc, aes(x = "", y = log_enrichment)) +
  geom_boxplot() +
  labs(title = "Distribution of log2(Relative Enrichment)")


ggplot(enrichment_qc, aes(x = fusil, y = log_enrichment)) +
  geom_boxplot() +
  labs(title = "log2(Relative Enrichment) by FUSIL")




# Remove outliers

cleaned_enrichement <- enrichment_qc %>%
  filter(n_x >10) %>%
  filter(
    (fusil == "CL"  & log_enrichment >= -2   & log_enrichment <= 2)   |
      (fusil == "DL"  & log_enrichment >= -1   & log_enrichment <= 1)   |
      (fusil == "SV"  & log_enrichment >= -1   & log_enrichment <= 1)   |
      (fusil == "VP"  & log_enrichment >= -0.5 & log_enrichment <= 0.5) |
      (fusil == "VnP" & log_enrichment >= -1   & log_enrichment <= 1)
  )


# plot afetr cleaning 

ggplot(cleaned_enrichement, aes(x = fusil, y = log_enrichment)) +
  geom_boxplot() +
  labs(title = "log2(Relative Enrichment) by FUSIL")

# Visualize with heatmaply



# Reshape into matrix
heatmap_matrix <- cleaned_enrichement %>%
  select(fusil, HLGT, relative_enrichment) %>%
  pivot_wider(names_from = fusil, values_from = relative_enrichment, values_fill = 1) %>%
  column_to_rownames("HLGT")%>%
  as.matrix()

# Plot interactive heatmap
heatmaply(heatmap_matrix,
          xlab = "FUSIL", ylab = "HLGT",
          colors = colorRampPalette(c("blue", "white", "red"))(100),
          k_col = 0, k_row = 0)  



# SE group enrichement for each Each HLGT–Gene Pair compared to the relative overall side effects -----------------------------


count_gene_se2 <- mhra_df %>%
  count(approvedSymbol)

se_group2 <- mhra_df %>%
  group_by(approvedSymbol, HLGT) %>%
  count(approvedSymbol) %>%
  left_join(count_gene_se2, by = "approvedSymbol") %>%
  mutate(percentage = (n.x*100)/n.y)



#Step 1 : Group by gene and HLGT and sum counts

group_counts2 <- se_group2 %>%
  group_by(approvedSymbol, HLGT) %>%
  summarise(n_x = sum(n.x), .groups = 'drop')

# Step 2: Total counts per Gene
gene_totals2 <- group_counts2 %>%
  group_by(approvedSymbol) %>%
  summarise(gene_total = sum(n_x))

# Step 3: Merge and calculate proportion within Gene
group_counts_merged2 <- group_counts2 %>%
  left_join(gene_totals2, by = "approvedSymbol") %>%
  mutate(gene_prop = n_x / gene_total)

# Step 4: Overall HLGT totals and proportions
overall_totals2 <- se_group2 %>%
  group_by(HLGT) %>%
  summarise(overall_total = sum(n.x)) 

# Step 5: Merge and compute relative enrichment
enrichment_df2 <- group_counts_merged2 %>%
  left_join(overall_totals2, by = "HLGT") %>%
  mutate(overall_prop = n_x /overall_total) %>%
  mutate(relative_enrichment = gene_prop / overall_prop ,
         log2_enrichment = log2(relative_enrichment)) %>%
  filter(n_x >= 50)





# Interactve Plot for HLGT-GENE enrichement 
enrichement_arranged <- enrichment_df2 %>%
  arrange(desc(log2_enrichment))

p <- ggplot(enrichement_arranged, aes(
  x = approvedSymbol,
  y = HLGT,
  size = gene_prop,
  color = log2_enrichment,
  text = paste(
    "Gene:", approvedSymbol,
    "<br>HLGT:", HLGT,
    "<br>log2 enrichment:", round(log2_enrichment, 2),
    "<br>Gene prop:", gene_prop
  )
)) +
  geom_point(alpha = 0.8) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    title = "Gene Enrichment by HLGT",
    x = "Gene",
    y = "HLGT"
  )

# Make it interactive
ggplotly(p, tooltip = "text")


# "Parallelized Fisher's Exact Test and Volcano Plot for Gene–Side Effect Enrichment Analysis" --------


# Load required libraries
library(dplyr)
library(purrr)
library(progressr)
library(furrr)

# Set up parallel backend and progress bar handler
plan(multisession)              # Use all available cores (Windows/macOS)
handlers("txtprogressbar")      # Enable terminal progress bar

# Compute total side effects once
total_side_effects <- sum(overall_totals$overall_total)

# Precompute contingency table values
enrichment_df2 <- enrichment_df2 %>%
  mutate(
    a = n_x,
    b = gene_total - n_x,
    c = overall_total - n_x,
    d = total_side_effects - gene_total - c
  )

# Run Fisher's test in parallel with progress bar
enrichment_df_fisher <- with_progress({
  p <- progressor(steps = nrow(enrichment_df2))
  
  enrichment_df2 %>%
    mutate(
      fisher_p = future_pmap_dbl(
        list(a, b, c, d),
        function(a, b, c, d) {
          p()  # update progress bar
          if (all(c(a, b, c, d) >= 0)) {
            fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
          } else {
            NA_real_
          }
        }
      )
    )
})


enrichment_df_fisher <- enrichment_df_fisher %>%
  mutate(fdr = p.adjust(fisher_p, method = "BH")) %>%
  filter(fdr <= 0.05)





# Add significance and labeling columns
enrichment_df_fisher_cleaned <- enrichment_df_fisher %>%
  mutate(
    neg_log10_fdr = -log10(fdr),
    significant = ifelse(fdr < 0.05 & abs(log2_enrichment) >= 1, "Yes", "No")
  ) %>%
  filter(neg_log10_fdr < 317)

# Plot
ggplot(enrichment_df_fisher_cleaned, aes(x = log2_enrichment, y = neg_log10_fdr)) +
  geom_point(aes(color = significant), alpha = 0.6) +
  scale_color_manual(values = c("No" = "grey70", "Yes" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot of Gene–Side Effect Enrichment",
    x = "Log2 Enrichment",
    y = "-Log10 FDR",
    color = "Significant"
  ) +
  theme_minimal()









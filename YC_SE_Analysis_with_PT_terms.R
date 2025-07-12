

## Script to Understand Side Effects from Yellow Card MHRA with Gene Essentiality and Drug Targets data ##


# Load the necessary libraries ------------------------


library(tidyverse)
library(heatmaply)


# Import the required Datasets ------------------------


#Import SE Dataset

mhra_se_data <- read_csv("C:/Users/HP-ssd/Desktop/Short term project2/mhra/mhra_se_data_with_soc.csv", , 
                         col_types = cols(...1 = col_skip()))

mhra_se_data <- mhra_se_data %>%
  na.omit() %>%
  filter(seriousness != "FALSE")

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





# Count of Individual PT side effects term for each FUSIL and observe their seriousness distribution --------

se_count_pt_serious <- mhra_se_data %>%
  filter(DRUG %in% drug_fusil$prefName) %>%
  select(1,2,3,4,5,9) %>%
  group_by(DRUG, seriousness) %>%
  tally()

se_count_serious_fusil <- se_count_pt_serious %>%
  left_join(drug_fusil, by = c("DRUG"= "prefName" )) %>%
  ungroup() %>%
  select(-4) %>%
  unique()
  
se_serious <- se_count_serious_fusil %>%
  group_by(fusil) %>%
  mutate(total_per_fusil = sum(n)) %>%
  group_by(fusil, seriousness) %>%
  summarise(n = sum(n), total = first(total_per_fusil)) %>%
  mutate(percentage = (n / total) * 100)

# Plot

se_serious$fusil <- factor(se_serious$fusil, 
                           levels = c("CL", "DL", "SV", "VP", "VnP"  ))

ggplot(se_serious, aes(x = fusil, y = percentage, fill = seriousness)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Percentage of Serious Side Effects in each FUSIL",
       x = "FUSIL bin",
       y = "Percentage") +
  theme_minimal()



# #Compute Enrichment for Each PT–FUSIL Pair --------


se_count_pt <- mhra_se_data %>%
  filter(DRUG %in% drug_fusil$prefName) %>%
  select(1,2,3,4,5,9) %>%
  group_by(DRUG, PT) %>%
  tally()

se_count_pt_fusil <- se_count_pt %>%
  left_join(drug_fusil, by = c("DRUG"= "prefName" )) %>%
  ungroup() %>%
  select(-4) %>%
  unique()


# Step 1: Group by FUSIL and PT and sum counts

group_counts <- se_count_pt_fusil %>%
  group_by(fusil, PT) %>%
  summarise(n_x = sum(n), .groups = 'drop')


# Step 2: Total counts per FUSIL
fusil_totals <- group_counts %>%
  group_by(fusil) %>%
  summarise(fusil_total = sum(n_x))

# Step 3: Merge and calculate proportion within FUSIL
group_counts <- group_counts %>%
  left_join(fusil_totals, by = "fusil") %>%
  mutate(fusil_prop = n_x / fusil_total)

# Step 4: Overall HLGT totals and proportions
overall_totals <- se_count_pt_fusil %>%
  group_by(PT) %>%
  summarise(overall_total = sum(n)) %>%
  mutate(overall_prop = overall_total / sum(overall_total))

# Step 5: Merge and compute relative enrichment
enrichment_df <- group_counts %>%
  left_join(overall_totals, by = "PT") %>%
  mutate(relative_enrichment = fusil_prop / overall_prop)


# Visualsie outliers

enrichment_qc <- enrichment_df %>%
  mutate(log_enrichment = log2(relative_enrichment)) %>%
  filter(n_x >= 10)

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
    (fusil == "CL"  & log_enrichment >= -2.5   & log_enrichment <= 2.5)   |
      (fusil == "DL"  & log_enrichment >= -1   & log_enrichment <= 1.3)   |
      (fusil == "SV"  & log_enrichment >= -1.25   & log_enrichment <= 1.25)   |
      (fusil == "VP"  & log_enrichment >= -0.5 & log_enrichment <= 0.5) |
      (fusil == "VnP" & log_enrichment >= -1.25   & log_enrichment <= 1.25)
  )


# plot after cleaning 

ggplot(cleaned_enrichement, aes(x = fusil, y = log_enrichment)) +
  geom_boxplot() +
  labs(title = "log2(Relative Enrichment) by FUSIL")



# Reshape into matrix
heatmap_matrix <- cleaned_enrichement %>%
  select(fusil, PT, relative_enrichment) %>%
  pivot_wider(names_from = fusil, values_from = relative_enrichment, values_fill = 1) %>%
  column_to_rownames("PT")%>%
  as.matrix()

# Plot interactive heatmap
heatmaply(heatmap_matrix,
          xlab = "FUSIL", ylab = "PT",
          colors = colorRampPalette(c("blue", "white", "red"))(100),
          k_col = 0, k_row = 0)  



# Include Higher Biological Term for Better Segreggation  ---------------------



higher_term <- mhra_se_data %>%
  select(PT, SOC_ABBREV) %>%
  unique()


enrichement_soc_df <- cleaned_enrichement %>%
  left_join(higher_term, by = "PT")

##For LOOP for each SOC Term

# Get unique SOC_ABBREV values
soc_list <- unique(enrichement_soc_df$SOC_ABBREV)

# Loop over each SOC_ABBREV
for (soc in soc_list) {
  
  # Filter the data for the current SOC_ABBREV
  subset_data <- enrichement_soc_df %>%
    filter(SOC_ABBREV == soc)
  
  # Reshape into matrix
  heatmap_matrix2 <- subset_data %>%
    select(fusil, PT, relative_enrichment) %>%
    pivot_wider(names_from = fusil, values_from = relative_enrichment, values_fill = 1) %>%
    column_to_rownames("PT") %>%
    as.matrix()
  
  # Plot interactive heatmap
  heatmaply(heatmap_matrix2,
            xlab = "FUSIL", ylab = "PT",
            colors = colorRampPalette(c("blue", "white", "red"))(100),
            k_col = 0, k_row = 0,
            main = paste("SOC:", soc))
}


























# #Compute Enrichment for Each HLGT–FUSIL Pair --------


se_count_hlgt <- mhra_se_data %>%
  filter(DRUG %in% drug_fusil$prefName) %>%
  select(1,2,3,4,7,9) %>%
  group_by(DRUG, HLGT) %>%
  tally()

se_count_hlgt_fusil <- se_count_hlgt %>%
  left_join(drug_fusil, by = c("DRUG"= "prefName" )) %>%
  ungroup() %>%
  select(-4) %>%
  unique()


# Step 1: Group by FUSIL and HLGT and sum counts

group_counts_hlgt <- se_count_hlgt_fusil %>%
  group_by(fusil, HLGT) %>%
  summarise(n_x = sum(n), .groups = 'drop')


# Step 2: Total counts per FUSIL
fusil_totals_hlgt <- group_counts_hlgt %>%
  group_by(fusil) %>%
  summarise(fusil_total = sum(n_x))

# Step 3: Merge and calculate proportion within FUSIL
group_counts_hlgt <- group_counts_hlgt %>%
  left_join(fusil_totals_hlgt, by = "fusil") %>%
  mutate(fusil_prop = n_x / fusil_total)

# Step 4: Overall HLGT totals and proportions
overall_totals_hlgt <- se_count_hlgt_fusil %>%
  group_by(HLGT) %>%
  summarise(overall_total = sum(n)) %>%
  mutate(overall_prop = overall_total / sum(overall_total))

# Step 5: Merge and compute relative enrichment
enrichment_df_hlgt <- group_counts_hlgt %>%
  left_join(overall_totals_hlgt, by = "HLGT") %>%
  mutate(relative_enrichment = fusil_prop / overall_prop)


# Visualsie outliers

enrichment_qc_hlgt <- enrichment_df_hlgt %>%
  mutate(log_enrichment = log2(relative_enrichment)) %>%
  filter(n_x >= 10)

enrichment_qc_hlgt$fusil <- factor(enrichment_qc_hlgt$fusil, 
                              levels = c("CL", "DL", "SV", "VP", "VnP"  ))

ggplot(enrichment_qc_hlgt, aes(x = "", y = log_enrichment)) +
  geom_boxplot() +
  labs(title = "Distribution of log2(Relative Enrichment)")


ggplot(enrichment_qc_hlgt, aes(x = fusil, y = log_enrichment)) +
  geom_boxplot() +
  labs(title = "log2(Relative Enrichment) by FUSIL")


# Remove outliers

cleaned_enrichement_hlgt <- enrichment_qc_hlgt %>%
  filter(n_x >10) %>%
  filter(
    (fusil == "CL"  & log_enrichment >= -2.5   & log_enrichment <= 2.5)   |
      (fusil == "DL"  & log_enrichment >= -1   & log_enrichment <= 1)   |
      (fusil == "SV"  & log_enrichment >= -1.25   & log_enrichment <= 1.25)   |
      (fusil == "VP"  & log_enrichment >= -0.4 & log_enrichment <= 0.4) |
      (fusil == "VnP" & log_enrichment >= -1.25   & log_enrichment <= 1.25)
  )


# plot after cleaning 

ggplot(cleaned_enrichement_hlgt, aes(x = fusil, y = log_enrichment)) +
  geom_boxplot() +
  labs(title = "log2(Relative Enrichment) by FUSIL")



# Reshape into matrix
heatmap_matrix_hlgt <- cleaned_enrichement_hlgt %>%
  select(fusil, HLGT, relative_enrichment) %>%
  pivot_wider(names_from = fusil, values_from = relative_enrichment, values_fill = 1) %>%
  column_to_rownames("HLGT")%>%
  as.matrix()

# Plot interactive heatmap
heatmaply(heatmap_matrix_hlgt,
          xlab = "FUSIL", ylab = "HLGT",
          colors = colorRampPalette(c("blue", "white", "red"))(100),
          k_col = 0, k_row = 0)  


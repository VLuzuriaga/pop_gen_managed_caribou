###############################################################################
### Load sample table, extracted from GenomeStudio 
  # This dataset has added information for sample type, sex and ecotype. 
samples <- read.table("SummaryTable_ALL_samples.txt", header=T, sep="\t")

### Visually examine call rates by sample type

ggplot(samples, aes(
  # This sorts samples by call rate
  x = reorder(Sample.ID,Call.Rate),
  y = Call.Rate,
  fill = Sample_type
)) +
  geom_jitter(
      shape = 21,     
      alpha = 0.8,
      stroke=0,
      size = 3
  ) +
# Highlight low quality thresholds
  # Samples with call rates <0.15 failed 
  geom_hline(yintercept = 0.15, linetype = "dashed", col="red3") +
# Samples with call rates >0.9 are of high quality
  geom_hline(yintercept = 0.90, linetype = "dashed", col="forestgreen") +
# Differentiate genotype data by sample type: blood vs feces
  scale_fill_manual(
    name = "Sample type", 
    values = c("Blood" = "orange3", "Feces" = "mediumpurple4")) +
  ylab("Call rate") + xlab("Sample ID") +
# This increases the size of the symbol in the legend
  guides(
    fill = guide_legend(override.aes = list(size = 5))
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=7))
# ----

###############################################################################
### Examine individuals with paired genotype data for blood and feces

# Sort samples by ID
paired <- samples[order(samples$Sample.ID), ]
# Find pairs with the sample GAN ID 
prefix <- sub("_.*$", "", paired$Sample.ID)
  # Extract index 
  idx <- which(duplicated(prefix) | duplicated(prefix, fromLast = TRUE))
  # Remove rows that are not repeated
  paired <- paired[idx,]

  # Extract paired individuals
paired2 <- paired %>%
  mutate(Individual = sub("_.*$", "", Sample.ID))

## Plot call rates by sample type per individual
ggplot(paired2, aes(
  x = Individual,
  y = Call.Rate,
  fill = Sample_type
)) +
  geom_col(position = position_dodge(width = 0.8), alpha=0.8) +
  labs(
    x = "Individual",
    y = "SNP call rate",
    fill = "Sample type"
  ) +
  scale_fill_manual(
    name = "Sample type", 
    values = c("Blood" = "orange3", "Feces" = "mediumpurple4")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
# ----

## Plot the total proportion of call rates per sample per individual

paired3 <- paired %>%
  mutate(Individual = sub("_.*$", "", Sample.ID)) %>%
  group_by(Individual) %>%
  filter(n() > 1)

ggplot(paired3, aes(
  x = Individual,
  y = Call.Rate,
  fill = Sample_type,
  # This highlights samples that failed to pass minimum QC
  alpha = Call.Rate >= 0.15
)) +
  geom_col(position = "fill") +
  scale_alpha_manual(
    values = c("TRUE" = 0.8, "FALSE" = 0.2),
    guide = "none"
  ) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Individual ID",
    y = "Proportion of total call rate",
    fill = "Sample type"
  ) +
  scale_fill_manual(
    name = "Sample type",
    values = c("Blood" = "orange3", "Feces" = "mediumpurple4")
  ) +
  geom_hline(yintercept = 0.15, linetype = "dashed", color = "black") +
  theme_bw() +
  coord_flip()
# ----

###############################################################################
### Extract failed samples to remove from downstream analyses 
failed_samples <- subset(samples, samples$Call.Rate<0.15)

###############################################################################
### Plot proportions of SNPs in failed samples

qc_long_counts <- qc_data %>%
  pivot_longer(
    cols = c(Num.good.SNPs, Num.nocall.SNPs),
    names_to = "CallType",
    values_to = "Count"
  ) %>%
  mutate(
    CallType = recode(CallType,
                      Num.good.SNPs   = "High quality",
                      Num.nocall.SNPs = "Missing calls")
  )
## sort by quality

qc_long_counts$Sample <- factor(
  qc_long_counts$Sample,
  levels = qc_data$Sample[order(qc_data$Num.good.SNPs)]
)

## Plot 
ggplot(qc_long_counts, aes(x = Sample, y = Count, fill = CallType)) +
  geom_col(show.legend = FALSE, alpha=0.7) +
  facet_wrap(~ CallType, ncol = 2, scales = "free_x") +
  scale_fill_manual(values = c(
    "High quality" = "forestgreen",
    "Missing calls"  = "firebrick"
  )) +
  labs(
    x = NULL,
    y = "Number of SNPs"
  ) +
  coord_flip() +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 9)
  )

###############################################################################
### Pull out proportions of SNPs of high versus low quality from failed samples 
  # Load full table per SNP array from GenomeStudio
full_data_c1 <- read.table("Full_Data_Table_chip1.txt", header = T, sep="\t")
full_data_c23 <- read.table("Full_Data_Table_chip2-3.txt", header = T, sep="\t")

## Function
calc_snp_stats <- function(data, score_cols, sample_ids, threshold = 0.90) {
  results <- lapply(score_cols, function(col) {
    scores <- data[[col]]
    good_snps <- sum(scores > threshold, na.rm = TRUE)
    nocall_snps <- sum(is.nan(scores))
    return(c(good_snps, nocall_snps))
  })
  results <- do.call(rbind, results)
  # Store results in a dataframe
  prop_failed <- data.frame(
    Sample = sample_ids,
    Num.good.SNPs = results[,1],
    Num.nocall.SNPs = results[,2]
  )
  return(prop_failed)
}
# ----

## For samples in array 1
score_cols_1 <- c(
"X24588Lola_MHZ204729310007R03C01.Score",
"X24588Rosa_TZ204729310007R07C01.Score")

failed_prop1 <- calc_snp_stats(
  data = full_data_c1,
  score_cols = score_cols_1,
  sample_ids = failed_samples$Sample.ID[1:2]
)

## For samples in arrays 2 and 3
score_cols_2 <- c(
  "X24588Suri_APZ_F25204729310013R12C01.Score",
  "X24588Neige_PO_F34204729310013R09C02.Score",
  "X24588VRM24.00100_BCWP_F21204729310013R08C01.Score",
  "X24588Victor_APZ_F27204729310013R02C02.Score",
  "X24588Theodore_APZ_F24204729310013R11C01.Score",
  "X24588Karma_PO_F33204729310013R08C02.Score",
  "X24588F20_BWP_VRM24.00703204729310008R02C01.Score",
  "X24588Raina_GVZ_F31204729310013R06C02.Score",
  "X24588Whitney_APZ_F28204729310013R03C02.Score",
  "X24588Charlotte_PO_F32204729310013R07C02.Score",
  "X24588Freddy_APZ_F29204729310013R04C02.Score",
  "X24588Avalanche_APZ_F22204729310013R09C01.Score",
  "X24588LBJ21.11767_Primrose204729310013R11C02.Score"
)

failed_prop2 <- calc_snp_stats(
  data = full_data_c23,
  score_cols = score_cols_2,
  sample_ids = failed_samples$Sample.ID[3:15]
)

failedProps <- rbind(failed_prop1, failed_prop2)
# ----

## Plot Plot props of SNPs of failed samples
qc_long_counts <- failedProps %>%
  pivot_longer(
    cols = c(Num.good.SNPs, Num.nocall.SNPs),
    names_to = "CallType",
    values_to = "Count"
  ) %>%
  mutate(
    CallType = recode(CallType,
                      Num.good.SNPs   = "High quality",
                      Num.nocall.SNPs = "Missing calls")
  )

## Sort by quality
qc_long_counts$Sample <- factor(
  qc_long_counts$Sample,
  levels = failedProps$Sample[order(failedProps$Num.good.SNPs)]
)

## Plot 
ggplot(qc_long_counts, aes(x = Sample, y = Count, fill = CallType)) +
  geom_col(show.legend = FALSE, alpha=0.7) +
  facet_wrap(~ CallType, ncol = 2, scales = "free_x") +
  scale_fill_manual(values = c(
    "High quality" = "forestgreen",
    "Missing calls"  = "firebrick"
  )) +
  labs(
    x = NULL,
    y = "Number of SNPs"
  ) +
  coord_flip() +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 9)
  )
# ----

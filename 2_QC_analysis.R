###############################################################################
### Load sample table, extracted from GenomeStudio 
  # This dataset has added information for sample type and sex.
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
### Use this function to extract proportions of SNPs of high versus low quality 

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

## For example, for failed samples Lola and Rosa
score_cols_1 <- c(
"X24588Lola_MHZ204729310007R03C01.Score",
"X24588Rosa_TZ204729310007R07C01.Score")

failed_prop <- calc_snp_stats(
  data = full_data_c1, # Load full table per SNP array from GenomeStudio
  score_cols = score_cols_1, # Vector containing the column names of sample scores
  sample_ids = failed_samples$Sample.ID[1:2]
)

# ----

## Plot Plot props of SNPs of failed samples
qc_long_counts <- failed_prop %>%
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
  levels = failed_prop$Sample[order(failed_prop$Num.good.SNPs)]
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

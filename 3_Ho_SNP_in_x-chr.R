####################################################################
### Calculate Ho of SNPs in x-chr to determine sex 

## Load vcf file, obtained from Genome Quebe's Nanuq repository
vcf <- read.vcfR("ALL_samples.vcf")

## Remove failed samples from vcf file
 ## Extract IDs from vcf file
  vcf_samples <- colnames(vcf@gt)[-1]
## Remove failed samples
  ## Delete blank spaces from failed samples
  failed_samples$Sample.ID <- gsub("\\s+", "", failed_samples$Sample.ID)
  ## match IDs in vcf file with IDs of failed samples
  idx <- match(vcf_samples, failed_samples$Sample.ID)
  ## Remove samples 
  indx_to_remove <- which(!is.na(idx))
  ## add 1 to the index to keep the "format" column
  idx <- indx_to_remove + 1
## Filter vcf file 
vcf_filt <- vcf[, -idx]

## -----

####################################################################
## Load probe ID for SNPs in x-chr, available from Carrier et al. 2022
x_chr <- read.table("x-chr_SNPs.txt", header=F, sep="\t")

## Extract SNP IDs from vcf file 
vcf_snps <- vcf_filt@fix[, "ID"]
## Match VCF SNP IDs with IDs in x-chr
ind_x <- match(vcf_snps, x_chr$V1)
## Locate SNPs in x-chr 
ind_x_only <- which(!is.na(ind_x))
## Subset vcf file
vcf_X <- vcf_filt[ind_x_only,]
## Extract genotypes (exclude "format" column)
gt <- extract.gt(vcf_X, element = "GT")  # returns a matrix
rownames(gt)

## -----

### Function to compute Ho for each SNP in x-chr (n=1,626)
ho_per_snp <- function(x) {
  # remove NAs
  x <- x[!is.na(x)]
  if(length(x) == 0) return(NA)   # no data
  mean(x %in% c("0/1", "1/0"))
}
# apply row-wise
ho <- apply(gt, 1, ho_per_snp)

### Calculate Ho per individual 
ho_ind <- apply(gt, 2, function(x) {
  x <- x[!is.na(x)]
  mean(x %in% c("0/1", "1/0"))
})
## -----

### Extract actual individual sex 
index_sex <- match(colnames(gt), samples$Sample.ID)
index_sex[22] <- 7

sex <- samples[index_sex,"Sex"]

## Create a data frame
ho_df <- data.frame(
  Individual = colnames(gt),
  Ho = ho_ind,
  Sex = sex 
)
## -----

### Calculate % of missing genotypes
miss_g <- function(x) {sum(is.na(x)) / length(x)}
## per individual
prop_miss <- apply(gt, 2, miss_g)
## add to data frame
ho_df$Miss <- prop_miss
## -----

## Plot Ho against prop of missing genotypes
ggplot(ho_df, aes(
  x = Ho,
  y = Miss,
  fill = Sex
)) +
  geom_jitter(
    shape = 21, 
    stroke=0,
    alpha = 0.8,
    size = 4
  ) +
  ## Homozygous individuals are males 
  geom_vline(xintercept = 0.05, linetype = "dashed", col="black") +
  labs(
    x = "X-chr Ho",
    y = "Missing genotypes (%)",
    fill = "Sex"
  ) +
  scale_fill_manual(
    name = "Sex", 
    values = c("female" = "magenta3", "male" = "royalblue"), 
    labels = c("Female", "Male")) +
  theme_bw() +
  theme()

## -----


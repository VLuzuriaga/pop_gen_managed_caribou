## The following libraries must be loaded: 
library(vcfR)
library(SNPRelate)
library(adegenet)
library(HardyWeinberg)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(mclust)                                  
library(pegas)
library(ade4)
library(ape)
library(ggtree)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(Rsamtools)
library(Biostrings)
library(GenomicRanges)

## Load input files
in_files <- list(
  samples  = "samples.csv",                 # sample tables extracted from the array data in GenomeStudio
  vcf      = "caribou_snps_raw.vcf",        # Raw genotype data converted into vcf format with vcfR
  bam      = "snps.bam",                    # Output BAM file after mapped array SNPs to the reference genome 
  fna      = "GCA_019903745.2_ULRtarCaribou_2v2_genomic.fna",  # chromosome-level reference genome
  no_call  = "No_call_SNPs.txt",            # SNPs excluded in Carrier et al. 2022 
  mono     = "Monomorphic_SNPs.txt",        # SNPs excluded in Carrier et al. 2022 
  x_probes = "x-chr_SNPs.txt",              # SNPs mapping to the X chromosome (Carrier et al. 2022)
  ecotype  = "hq_uniques.txt"               # classifier output (co-author)
)

## Parameters
P <- list(
  cr_fail  = 0.15,   # call rate < this = failed, excluded from everything
  cr_high  = 0.90,   # call rate >= this = high quality (Dataset 2)
  miss_max = 0.05,   # per-SNP missingness, computed in the selection samples only
  mac_min  = 3,      # minor allele count
  ld_win   = 50,     # LD window (number of SNPs)
  ld_r2    = 0.5,    # r^2 threshold (PLINK --indep-pairwise 50 5 0.5 equivalent)
  ho_x_male = 0.05,  # X-chromosome Ho below this = male
  n_expected = c(DS1_blood = 44, DS1_feces = 44, DS2 = 34)  # edit if final numbers differ
)


## Kinship categories (KING cut-offs)
kin_breaks <- c(-Inf, 0.0442, 0.0884, 0.177, 0.354, Inf)
kin_labels <- c("Unrelated", "Third degree", "Second degree", "First degree", "Duplicate")
kin_cat    <- function(phi) cut(phi, kin_breaks, labels = kin_labels, right = FALSE)

## Palettes (same colours in every figure)
inst_levels <- c("YWP","BCWP","GVZ","EVZ","CZ","APZ","SN","TZ","RPZ","ZSSF","PO","MHZ")
pal_inst <- setNames(brewer.pal(12, "Paired"), inst_levels)
pal_type <- c(Blood = "orange3", Feces = "mediumpurple4")
pal_qual <- c("High quality" = "grey85", "Lower confidence" = "black")
pal_clu  <- c("1" = "#1b9e77", "2" = "#d95f02", "3" = "#7570b3")

## Helpers
save_obj <- function(obj, name) saveRDS(obj, file.path(dir_derived, paste0(name, ".rds")))
load_obj <- function(name)      readRDS(file.path(dir_derived, paste0(name, ".rds")))
open_gds <- function()          snpgdsOpen(file.path(dir_derived, "caribou_all.gds"))

## Individual ID = first "_"-separated field of the sample name
ind_of <- function(s) sub("_.*$", "", s)

## Drop samples whose kinship row is mostly NaN; needed before hclust().
kin_clean <- function(K) {
  bad <- rownames(K)[rowMeans(is.na(K)) > 0.5]
  if (length(bad)) message("Dropped from kinship matrix: ", paste(bad, collapse = ", "))
  K[!rownames(K) %in% bad, !colnames(K) %in% bad]
}

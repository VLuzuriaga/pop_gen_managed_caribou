#########################################################################################
## I only kept high quality genotypes of invidivuals with blood and fecal samples.
## Removed all failed samples (n=15).
#########################################################################################

### Load filtered data with good quality individuals 
vcf_filt_QC <-  read.vcfR("HighQ_Individuals.vcf.gz")
  # 44 individuals

### Convert vcf file to Genomic Data Structures (GDS)
snpgdsVCF2GDS(
  "HighQ_Individuals.vcf.gz",
  "genotypes_HQ.gds",
  method = "biallelic.only"
)

### Open genofile
genofile <- snpgdsOpen("genotypes_HQ.gds")

# Remove SNPs with >20% missingnes
snpset <- snpgdsSelectSNP(
  genofile,
  sample.id = read.gdsn(index.gdsn(genofile, "sample.id")),
  autosome.only = FALSE,
  missing.rate = 0.20
)   # Excluding 5,130 SNPs 

## Check number of SNPs
length(snpset)  # 37K

### Calculate kinship using identity-by-descent (IBD)
kinship <- snpgdsIBDKING(
  genofile,
  snp.id = snpset, 
  autosome.only = FALSE)

  ## Extract kinship matrix
  kin_mat <- kinship$kinship

  ## Assign row and column names
  rownames(kin_mat) <- kinship$sample.id
  colnames(kin_mat) <- kinship$sample.id

  # Fill diagonals with self-kinship = 0.5
  diag(kin_mat) <- 0.5

  write.table(kin_mat, "Kinship_matrix.txt", sep="\t", col.names = T, row.names = T)

##-------
### PLOTS 

## Heatmap of pairwise kisnhip estimates 

pheatmap(
  kin_mat,
  clustering_distance_rows = as.dist(1 - kin_mat),
  clustering_distance_cols = as.dist(1 - kin_mat),
  border_color = NA, 
  fontsize=8,
  color= viridis(44, option = "C"))

## Phylogeny with genetic distances

  ## Calculate genetic distance
  dist_mat <- as.dist(0.5 - kin_mat)
    ## distance = 0.5−kinship
    #duplicates (kinship ≈ 0.5) → distance ≈ 0
    #unrelated (kinship ≈ 0) → distance ≈ 0.5

  ## Hierarchical clustering
  hc <- hclust(dist_mat, method = "average")  # UPGMA

## Plot 
tree <- as.phylo(hc)

plot_tree <- ggtree(tree) +
  geom_tiplab(size = 3, offset = 0.01) +
  theme_tree2() +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5.5, 95, 5.5, 5.5))
 
## reverse x-axis
  revts(plot_tree) + scale_x_continuous(labels=abs)

##-------







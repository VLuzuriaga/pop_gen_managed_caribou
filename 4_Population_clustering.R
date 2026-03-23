##############################################################
### Discriminant Analysis of Principal Components

### Convert VCF file into genind 
gen_ind <- vcfR2genind(vcf_filt, return.alleles = TRUE) 
  #markers with no scored alleles have been removed

# remove loci with >20% missing
fil_genind <- gen_ind[ , colMeans(is.na(tab(gen_ind))) < 0.2]

### Find clusters
  # I am limiting the number of clusters to 6, corresponding to 5 populations used in the SNP chip: 
    # migratory, mountain, western boreal, center boreal, eastern boreal 
    # plus Eurasian reindeer, present in the samples

groups <- find.clusters(fil_genind, max.n.clust = 6)
  # Retain all PCs. The plot is showing that the cummulative variation platoes at 50 PCs so will put *100* to keep them all 
  # The number of clusters that minimizes the BIC (Bayesian Information Criterion) is *4*
## Check groups size
groups$size

### Perform DAPC
dapc1 <- dapc(fil_genind, groups$grp)
  #Based on the cummulative variance figure, after 35 PCs little variation is added. Retained *35* PCs
  #Retained all *3* discriminant functions. 
## Check information
dapc1

## Extract table with clustering groups
clusters <- dapc1$grp
write.table(clusters, "Clusters_55_samples_35PCs.txt", sep="\t")

## Check % of contribution of each eigen value
ld_contrib <- dapc1$eig / sum(dapc1$eig) * 100

## Check clusters
df_dapc[order(df_dapc$cluster), ]

##------

### Plot DAPC

## Extract scores
df_dapc <- data.frame(
  ind = rownames(dapc1$ind.coord),
  LD1 = dapc1$ind.coord[, 1],
  LD2 = dapc1$ind.coord[, 2],
  cluster = groups$grp
)

## Create start/end for each cluster connection
segments <- df_dapc %>%
  group_by(cluster) %>%
  mutate(
    centroid_LD1 = mean(LD1),
    centroid_LD2 = mean(LD2)
  ) %>%
  ungroup()

## Calculate centroids for cluster labels
centroids <- df_dapc %>%
  group_by(cluster) %>%
  summarise(
    LD1 = mean(LD1),
    LD2 = mean(LD2)
  )

### Plot
ggplot(df_dapc, aes(x = LD1, y = LD2, colour = factor(cluster))) +
  #Points
  geom_point(size = 3, alpha = 0.8) +
  #Lines connecting points to their cluster centroid
  geom_segment(
    data = segments,
    aes(x = LD1, y = LD2, xend = centroid_LD1, yend = centroid_LD2),
    alpha = 0.8
  ) +
  #Ellipses
  stat_ellipse(type = "norm", linewidth = 0.3) +
  #Cluster labels
  geom_label(
    data = centroids,
    aes(x = LD1, y = LD2, label = cluster),
    fill = "white",
    size = 4, 
    fontface= "plain"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #Add labels
  labs(
    x = "Discriminant 1 (92.1%)",
    y = "Discriminant 2 (7.84%)",
    colour = "Cluster"
  ) +
  theme_bw() +
  theme(legend.position = "none")

##--------

## Optional: plot all loadings in 3D

library(plotly)

# Extract DAPC coordinates
dapc_df <- as.data.frame(dapc1$ind.coord)

# Add group/cluster info
dapc_df$cluster <- as.factor(dapc1$grp)

plot_ly(
  dapc_df,
  x = ~LD1,
  y = ~LD2,
  z = ~LD3,
  color = ~cluster,
  colors = "Set1",
  type = "scatter3d",
  mode = "markers"
) %>%
  layout(
    scene = list(
      xaxis = list(title = "D1 (73%)"),
      yaxis = list(title = "D2 (20.3%)"),
      zaxis = list(title = "D3 (6.7%)")
    ))

##--------





source("1_Libraries_to_load.R")
list2env(load_obj("sample_sets"), envir = .GlobalEnv)
ids_kin  <- load_obj("snpsets_ids")$Kinship
genofile <- open_gds()

king <- lapply(datasets, function(s)
  snpgdsIBDKING(genofile, sample.id = s, snp.id = ids_kin,
                autosome.only = FALSE, verbose = FALSE))
kin_mat <- function(k) {
  m <- k$kinship; dimnames(m) <- list(k$sample.id, k$sample.id)
  diag(m) <- 0.5; m
}
K_list <- lapply(king, kin_mat)

## Number of SNPs used and any samples with undefined (NaN) kinship
sapply(K_list, function(K) c(n_samples = nrow(K), n_snps = length(ids_kin)))
lapply(K_list, function(K) rownames(K)[rowMeans(is.na(K)) > 0.5])   \

## Pairs by relationship category (cut-offs defined in 1_Libraries_to_load.R)
pair_table <- function(K) {
  idx <- which(upper.tri(K), arr.ind = TRUE)
  kin_cat(K[idx])
}
lapply(K_list, function(K) table(pair_table(K)))

## Removal of one member of each close pair (kinship > 0.177)
drop_close <- function(K, thr = 0.177) {
  K <- kin_clean(K)
  K0 <- K; diag(K0) <- 0; out_s <- character()
  while (any(K0 > thr, na.rm = TRUE)) {
    i <- which.max(rowSums(K0 > thr, na.rm = TRUE))
    out_s <- c(out_s, rownames(K0)[i]); K0 <- K0[-i, -i]
  }
  out_s
}
dropped        <- lapply(K_list, drop_close)
datasets_noRel <- Map(setdiff, datasets, dropped)   # use for sensitivity runs
sapply(datasets_noRel, length)

for (d in names(K_list))
  write.csv(K_list[[d]], file.path(out, paste0("kinship_matrix_", d, ".csv")))

snpgdsClose(genofile)
save_obj(list(K_list = K_list, dropped = dropped, datasets_noRel = datasets_noRel), "kinship")

##############################################################################
## Figures

list2env(load_obj("sample_sets"), envir = .GlobalEnv)   # info, lowconf, paired_ids
list2env(load_obj("kinship"),     envir = .GlobalEnv)   # K_list

## Samples with undefined (NaN) kinship (e.g. the one with almost no overlapping
## genotypes in Dataset 1-feces) cannot be clustered: drop them for plotting only.
K_plot <- lapply(K_list, kin_clean)
ds_names <- c(DS2 = "Dataset 2", DS1_blood = "Dataset 1-blood", DS1_feces = "Dataset 1-feces")
star <- function(s) ifelse(s %in% lowconf, paste0(s, " *"), s)

## Hierarchical clustering (shared by heatmaps and trees)
hc_list <- lapply(K_plot, function(K) hclust(as.dist(0.5 - K), method = "average"))

## Pheatmaps (one per dataset). Scale from -0.5 to 0.5.
bk       <- seq(-0.5, 0.5, length.out = 101)
heat_col <- colorRampPalette(c("navy", "white", "firebrick"))(100)
for (d in names(K_plot)) {
  K  <- K_plot[[d]]
  ph <- pheatmap(K, cluster_rows = hc_list[[d]], cluster_cols = hc_list[[d]],
                 labels_row = star(rownames(K)), breaks = bk, color = heat_col,
                 show_colnames = FALSE, fontsize_row = 7, main = ds_names[[d]],
                 silent = TRUE)
  pdf(file.path(out, paste0("Fig_kinship_heatmap_", d, ".pdf")), width = 9, height = 8)
  grid::grid.newpage(); grid::grid.draw(ph$gtable); dev.off()
}

## ComplexHeatmap (lower triangle). 
## Paired individuals (blood + feces from the same animal) have coloured sample names.
col_fun <- colorRamp2(c(-0.5, 0, 0.5), c("navy", "white", "firebrick"))

kin_ht <- function(K, title, pair_col, show_legend = TRUE, keep_diag = TRUE) {
  hc  <- hclust(as.dist(0.5 - K), method = "average")
  K   <- K[hc$order, hc$order]
  lab_col <- ifelse(ind_of(rownames(K)) %in% paired_ids, pair_col, "black")
  rownames(K) <- colnames(K) <- star(rownames(K))
  K[upper.tri(K, diag = !keep_diag)] <- NA            # keep lower triangle only

  ComplexHeatmap::Heatmap(
    K, name = "KING\nkinship", col = col_fun, na_col = "white",
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_names_side = "left", column_names_side = "bottom",
    row_names_gp    = gpar(fontsize = 7, col = lab_col),
    column_names_gp = gpar(fontsize = 7, col = lab_col), column_names_rot = 90,
    column_title = title, column_title_side = "top",
    column_title_gp = gpar(fontsize = 12, fontface = "plain"),
    rect_gp = gpar(col = NA), show_heatmap_legend = show_legend,
    heatmap_legend_param = list(title_gp  = gpar(fontface = "plain", fontsize = 9),
                                labels_gp = gpar(fontface = "plain", fontsize = 8)))
}

## Dataset 2 on its own
ht_DS2 <- kin_ht(K_plot$DS2, "Dataset 2", pair_col = "black")
pdf(file.path(out, "kinship_DS2.pdf"), width = 9, height = 8)
ComplexHeatmap::draw(ht_DS2); dev.off()

## Dataset 1: blood and feces side by side
ht_b <- kin_ht(K_plot$DS1_blood, "Dataset 1-blood", pair_col = "goldenrod2",    show_legend = FALSE)
ht_f <- kin_ht(K_plot$DS1_feces, "Dataset 1-feces", pair_col = "mediumpurple4")
g_b  <- grid.grabExpr(ComplexHeatmap::draw(ht_b, padding = unit(c(2, 2, 2, 2), "mm")))
g_f  <- grid.grabExpr(ComplexHeatmap::draw(ht_f, padding = unit(c(2, 2, 2, 2), "mm")))
p_kin <- wrap_elements(g_b) + wrap_elements(g_f) +
  plot_annotation(caption = paste("Paired individuals (blood and feces from the same animal):",
                                  "yellow labels = blood, purple labels = feces.",
                                  "* lower-confidence sample"))
p_kin
ggsave(file.path(out, "kinship_DS1_blood_feces.png"), p_kin, width = 16, height = 8, dpi = 300)

##  UPGMA trees. Tips coloured by institution, shape = sample type.
make_tree <- function(K, hc, title, mirror = FALSE, extra = 0.6, off = 0.05) {
  m  <- match(rownames(K), info$Sample)
  md <- data.frame(label       = rownames(K),
                   Institution = factor(info$Institution[m], levels = inst_levels),
                   Sample      = info$Type[m],
                   tiplab      = star(rownames(K)))
  tr    <- as.phylo(hc); tr$edge.length <- tr$edge.length * 2
  xmax  <- max(node.depth.edgelength(tr))
  delta <- off * xmax

  p    <- ggtree(tr) %<+% md
  tips <- p$data[p$data$isTip, ]
  p <- p +
    geom_text(data = tips, aes(x = x + delta, y = y, label = tiplab, colour = Institution),
              hjust = if (mirror) 1 else 0, size = 3, show.legend = FALSE) +
    geom_tippoint(aes(colour = Institution, shape = Sample), size = 2) +
    scale_colour_manual(values = pal_inst, drop = TRUE) +
    expand_limits(x = xmax * (1 + extra)) +             # room for the full tip labels
    theme_tree2() +
    labs(x = "Distance (0.5 - kinship)", title = title) +
    theme(plot.title = element_text(hjust = 0.5))
  if (mirror) p <- p + scale_x_reverse()                # root on the right, labels on the left
  p
}
cap <- "* lower-confidence sample (call rate 0.15-0.90)"

p_tree_DS2  <- make_tree(K_plot$DS2,       hc_list$DS2,       "Dataset 2") + labs(caption = cap)
p_tree_DS1b <- make_tree(K_plot$DS1_blood, hc_list$DS1_blood, "Dataset 1-blood")
p_tree_DS1f <- make_tree(K_plot$DS1_feces, hc_list$DS1_feces, "Dataset 1-feces", mirror = TRUE)

## Dataset 1: blood normal, feces mirrored 
p_mirror <- (p_tree_DS1b + p_tree_DS1f) +
  plot_layout(guides = "collect") + plot_annotation(caption = cap)
p_mirror

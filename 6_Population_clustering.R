## DAPC following Thia (2023) and the Miller et al. (2020) checklist.

source("1_Libraries_to_load.R")
list2env(load_obj("sample_sets"), envir = .GlobalEnv)
list2env(load_obj("snp_data"),    envir = .GlobalEnv)   # snp, G
list2env(load_obj("kinship"),     envir = .GlobalEnv)   # K_list, datasets_noRel
vcf        <- load_obj("vcf")
ids_struct <- load_obj("snpsets_ids")$Structure

## Helpers
## genind for a set of samples, restricted to a SNP-id list (sid)
make_genind <- function(samples, ids) {
  rows <- snp$row[match(ids, snp$sid)]                # row in the ORIGINAL vcf
  cols <- c(1, match(samples, colnames(vcf@gt)))      # column 1 = FORMAT
  vcfR2genind(vcf[rows, cols], return.alleles = TRUE)
}

dapc_thia <- function(gi, k_max = 6, pcs = c(10, 20, 30), k = NULL, seed = 1) {
  n   <- nInd(gi)
  pcs <- sort(unique(c(pcs[pcs < n - 1], n - 1)))
  bic <- sapply(pcs, function(p) {
    set.seed(seed)
    find.clusters(gi, max.n.clust = k_max, n.pca = p, choose.n.clust = FALSE,
                  criterion = "min", n.start = 100)$Kstat
  })
  dimnames(bic) <- list(paste0("k", seq_len(k_max)), paste0("PC", pcs))
  k_bic <- which.min(bic[, ncol(bic)])           
  if (is.null(k)) k <- k_bic
  stopifnot(k >= 3)                               # need >= 2 discriminant functions
  set.seed(seed)
  fc <- find.clusters(gi, n.pca = n - 1, n.clust = k, n.start = 100)
  da <- dapc(gi, fc$grp, n.pca = k - 1, n.da = k - 1)   # k - 1 PCs, k - 1 DFs
  grp <- fc$grp; names(grp) <- indNames(gi)
  list(bic = bic, pcs = pcs, k = k, k_bic = k_bic, grp = grp, dapc = da,
       post_assigned = da$posterior[cbind(seq_len(n), as.integer(fc$grp))],
       n_reassigned  = sum(as.character(da$assign) != as.character(fc$grp)))
}

align_labels <- function(ref, new) {
  sh  <- intersect(names(ref), names(new))
  t2  <- table(ref[sh], new[sh]); map <- setNames(rep(NA_character_, ncol(t2)), colnames(t2))
  while (nrow(t2) > 0 && ncol(t2) > 0 && any(t2 > 0)) {
    ix <- which(t2 == max(t2), arr.ind = TRUE)[1, ]
    map[colnames(t2)[ix[2]]] <- rownames(t2)[ix[1]]
    t2 <- t2[-ix[1], -ix[2], drop = FALSE]
  }
  un <- is.na(map); map[un] <- paste0("new", seq_len(sum(un)))
  out <- factor(unname(map[as.character(new)]), levels = unique(c(levels(ref), map[un])))
  names(out) <- names(new); out
}

## DAPC on every dataset (k = 3)
gi_list  <- lapply(datasets, make_genind, ids = ids_struct)
dapc_res <- lapply(gi_list, dapc_thia, k = 3)
sapply(dapc_res, function(r) c(k = r$k, k_bic = r$k_bic, n_reassigned = r$n_reassigned))
sapply(dapc_res, function(r) ncol(r$dapc$ind.coord))     # CHECK: 2 for all

## Align clusters to Dataset 2
grp_ind <- lapply(dapc_res, function(r) { g <- r$grp; names(g) <- ind_of(names(g)); g })
ref     <- grp_ind$DS2
grp_al  <- lapply(grp_ind, function(g) align_labels(ref, g))

## Robustness to the number of PCs (Thia 2023)
pc_robust <- function(d) {
  r <- dapc_res[[d]]
  bind_rows(lapply(r$pcs, function(p) {
    set.seed(1)
    g <- find.clusters(gi_list[[d]], n.pca = p, n.clust = r$k, n.start = 100)$grp
    names(g) <- ind_of(indNames(gi_list[[d]]))
    g <- align_labels(grp_al[[d]], g)
    data.frame(dataset = d, n_pcs = p, ARI = adjustedRandIndex(g, grp_al[[d]][names(g)]),
               n_changed = sum(g != grp_al[[d]][names(g)]))
  }))
}
pc_check <- bind_rows(lapply(names(dapc_res), pc_robust)); pc_check


## ---------------------------------------------------------------------------
## D. BIC tables (Table S2)
## ---------------------------------------------------------------------------
lapply(dapc_res, function(r) round(r$bic, 2))      # rows = k, columns = PCs used

bic_long <- do.call(rbind, lapply(names(dapc_res), function(d) {
  b <- dapc_res[[d]]$bic
  data.frame(dataset = d,
             k       = rep(rownames(b), times = ncol(b)),
             PCs     = rep(colnames(b), each  = nrow(b)),
             BIC     = round(as.vector(b), 2),
             row.names = NULL)
}))
## difference from the best k at the same PC setting
bic_long$delta <- ave(bic_long$BIC, bic_long$dataset, bic_long$PCs, FUN = function(x) x - min(x))
subset(bic_long, k == "k6")                        # k = 6 alone, with delta from best k
write.csv(bic_long, file.path(out, "TableS2_BIC.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## E. Sensitivity: repeat without close relatives
## ---------------------------------------------------------------------------
res_noRel <- bind_rows(lapply(names(datasets_noRel), function(d) {
  r <- dapc_thia(make_genind(datasets_noRel[[d]], ids_struct), k = dapc_res[[d]]$k)
  g <- r$grp; names(g) <- ind_of(names(g)); g <- align_labels(grp_al[[d]], g)
  data.frame(dataset = d, n_individuals = length(g),
             ARI = adjustedRandIndex(g, grp_al[[d]][names(g)]),
             n_changed = sum(g != grp_al[[d]][names(g)]))
}))
res_noRel
write.csv(res_noRel, file.path(out, "dapc_sensitivity_no_relatives.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## F. Cluster composition and revised Table S1
## ---------------------------------------------------------------------------
dapc_df <- bind_rows(lapply(names(dapc_res), function(d) {
  co <- dapc_res[[d]]$dapc$ind.coord
  data.frame(Dataset = d,
             Sample  = rownames(co),
             LD1     = co[, 1],
             LD2     = if (ncol(co) >= 2) co[, 2] else NA_real_,
             PostAssigned = dapc_res[[d]]$post_assigned,
             row.names = NULL)
})) %>% mutate(Individual  = ind_of(Sample),
               Cluster     = unname(mapply(function(d, i) as.character(grp_al[[d]][i]), Dataset, Individual)),
               Institution = factor(info$Institution[match(Sample, info$Sample)], levels = inst_levels),
               Type        = info$Type[match(Sample, info$Sample)],
               LowConf     = Sample %in% lowconf)

cluster_comp <- dapc_df %>% dplyr::count(Dataset, Cluster, Institution, name = "n_samples")
cluster_comp
write.csv(cluster_comp, file.path(out, "dapc_cluster_composition.csv"), row.names = FALSE)

## Revised Table S1: one row per sample, call-rate class, dataset membership,
## cluster and discriminant scores per dataset. Define D1/D2 and the category
## labels in the caption.
wide <- dapc_df %>% select(Sample, Dataset, Cluster, LD1, LD2) %>%
  pivot_wider(names_from = Dataset, values_from = c(Cluster, LD1, LD2))
tabS1 <- info %>% filter(Category != "Failed") %>%
  mutate(in_DS1_blood = Sample %in% datasets$DS1_blood,
         in_DS1_feces = Sample %in% datasets$DS1_feces,
         in_DS2       = Sample %in% datasets$DS2) %>%
  left_join(wide, by = "Sample") %>% arrange(Institution, Individual, Type)
write.csv(tabS1, file.path(out, "TableS1_revised.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## G. Figures. Axes are arbitrary per dataset (free scales); report % variance
##    per discriminant from dapc_res[[d]]$dapc$eig in the caption.
## ---------------------------------------------------------------------------
p_clu <- ggplot(dapc_df, aes(LD1, LD2)) +
  stat_ellipse(aes(group = Cluster, colour = Cluster),
               type = "norm", linewidth = 0.3, show.legend = FALSE) +
  geom_point(aes(colour = Cluster, shape = Type), size = 2.5, alpha = 0.85) +
  geom_point(data = dplyr::filter(dapc_df, LowConf), shape = 8, size = 1.5, colour = "black") +
  scale_colour_manual(values = pal_clu) +
  labs(caption = "* lower-confidence sample (call rate 0.15-0.90)") +
  facet_wrap(~Dataset, scales = "free") + theme_bw()
p_clu

cent <- aggregate(cbind(LD1, LD2) ~ Dataset + Cluster, data = dapc_df, FUN = mean)
p_inst <- ggplot(dapc_df, aes(LD1, LD2)) +
  stat_ellipse(aes(group = Cluster, linetype = Cluster), colour = "grey40",
               type = "norm", linewidth = 0.4, show.legend = FALSE) +
  geom_point(aes(colour = Institution, shape = Type), size = 2.5, alpha = 0.85) +
  geom_point(data = dplyr::filter(dapc_df, LowConf), shape = 8, size = 1.5, colour = "black") +
  geom_label(data = cent, aes(label = Cluster), size = 3.5, fontface = "plain",
             colour = "grey20", fill = "white", alpha = 0.8,
             label.size = 0.2, label.padding = unit(0.15, "lines")) +
  scale_colour_manual(values = pal_inst) +
  scale_linetype_manual(values = c("1" = "dashed", "2" = "dashed", "3" = "dashed")) +
  guides(colour = guide_legend(ncol = 2, byrow = TRUE, order = 1,
                               override.aes = list(size = 3.5)),
         shape  = guide_legend(order = 2, override.aes = list(size = 3.5))) +
  labs(caption = "Numbers = DAPC cluster; * lower-confidence sample (call rate 0.15-0.90)") +
  facet_wrap(~Dataset, scales = "free") +
  theme_bw() +
  theme(legend.position = "right", legend.box = "vertical",
        legend.title = element_text(size = 9), legend.text = element_text(size = 8))
p_inst

ggsave(file.path(out, "FigS_DAPC_by_cluster.pdf"),     p_clu,  width = 12, height = 4.5)
ggsave(file.path(out, "FigS_DAPC_by_institution.pdf"), p_inst, width = 12, height = 4.5)

## % variance per discriminant function (for the captions)
lapply(dapc_res, function(r) round(100 * r$dapc$eig / sum(r$dapc$eig), 1))

## ---------------------------------------------------------------------------
## H. Save for 5_Genetic_relatedness.R and 6_Paired_blood_feces.R
## ---------------------------------------------------------------------------
save_obj(list(dapc_res = dapc_res, grp_al = grp_al, dapc_df = dapc_df), "dapc")

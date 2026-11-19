source("1_Libraries_to_load.R")
list2env(load_obj("sample_sets"), envir = .GlobalEnv)
list2env(load_obj("kinship"),     envir = .GlobalEnv)   # K_list
list2env(load_obj("dapc"),        envir = .GlobalEnv)   # grp_al

## Pairwise kinship, categories, within vs among institutions
pairs_df <- function(K, d) {
  idx <- which(upper.tri(K), arr.ind = TRUE)
  data.frame(a = rownames(K)[idx[, 1]], b = colnames(K)[idx[, 2]], phi = K[idx]) %>%
    mutate(ia = ind_of(a), ib = ind_of(b),
           inst_a = info$Institution[match(a, info$Sample)],
           inst_b = info$Institution[match(b, info$Sample)],
           cat    = kin_cat(phi),
           where  = ifelse(inst_a == inst_b, "Within institution", "Among institutions"),
           same_cluster = as.character(grp_al[[d]][ia]) == as.character(grp_al[[d]][ib]))
}
pairs_list <- Map(pairs_df, K_list, names(K_list))

## Long table of every pair, all datasets
kin_all <- bind_rows(lapply(names(pairs_list), function(d)
  pairs_list[[d]] %>%
    transmute(Dataset = d, Sample1 = a, Sample2 = b, Individual1 = ia, Individual2 = ib,
              Kinship = phi, Category = cat)))

## Related pairs only (second degree or closer)
kin_rel <- kin_all[!is.na(kin_all$Kinship) & kin_all$Kinship >= 0.0884, ]
kin_rel <- kin_rel[order(kin_rel$Dataset, -kin_rel$Kinship), ]

## Relationship categories x within/among, per dataset

table1 <- lapply(pairs_list, function(p) table(Category = p$cat, p$where))
table1

## Numbers for the text ("X of N within-institution pairs were first/second degree")
close_summary <- bind_rows(lapply(names(pairs_list), function(d)
  pairs_list[[d]] %>% group_by(where) %>%
    summarise(dataset = d, n_pairs = n(),
              n_first_or_second = sum(cat %in% c("First degree", "Second degree"), na.rm = TRUE),
              .groups = "drop")))
close_summary
write.csv(close_summary, file.path(out, "kinship_close_pairs_summary.csv"), row.names = FALSE)

## Are first-degree pairs in the same DAPC cluster?
first_deg_clusters <- bind_rows(lapply(names(pairs_list), function(d)
  pairs_list[[d]] %>% filter(cat == "First degree") %>%
    summarise(dataset = d, n_first = n(), n_same_cluster = sum(same_cluster))))
first_deg_clusters

## AMOVA on Dataset 2, institutions with >= 3 individuals
inst2 <- info$Institution[match(rownames(K_list$DS2), info$Sample)]
keep  <- inst2 %in% names(which(table(inst2) >= 3))
K2    <- K_list$DS2[keep, keep]
meta2 <- data.frame(Institution = factor(inst2[keep]))

D2 <- 1 - 2 * K2; diag(D2) <- 0; D2[D2 < 0] <- 0; D2 <- (D2 + t(D2)) / 2
D2d <- as.dist(D2)                       

ade4::is.euclid(as.dist(sqrt(D2)))      

set.seed(1)
am <- pegas::amova(D2d ~ Institution, data = meta2, nperm = 999, is.squared = TRUE)
am

vc <- as.data.frame(am$varcomp)          
s2 <- pmax(vc$sigma2, 0)                 # report negative sigma2 as ~0
phi_ST <- s2[1] / sum(s2)
amova_res <- c(percent_among = 100 * phi_ST, phi_ST = phi_ST, P = vc$P.value[1])
amova_res

## Concordance of DAPC clusters with institution
clu2 <- as.character(grp_al$DS2[ind_of(rownames(K2))])
fisher.test(table(meta2$Institution, clu2), simulate.p.value = TRUE, B = 1e5)
adjustedRandIndex(as.character(meta2$Institution), clu2)

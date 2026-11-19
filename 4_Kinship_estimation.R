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

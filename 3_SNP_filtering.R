source("1_Libraries_to_load.R")
list2env(load_obj("sample_sets"), envir = .GlobalEnv)   

## Load vcf
vcf <- read.vcfR(in_files$vcf, verbose = FALSE)

colnames(vcf@gt)[-1] <- gsub("\\s+", "", colnames(vcf@gt)[-1])
stopifnot(all(unlist(datasets) %in% colnames(vcf@gt)))

fx  <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)
snp <- data.frame(row = seq_len(nrow(fx)), id = fx$ID, chr = fx$CHROM,
                  pos = suppressWarnings(as.integer(fx$POS)),
                  ref = fx$REF, alt = fx$ALT, stringsAsFactors = FALSE)

## Remap SNPs to the chromosome-level assembly (GCA_019903745.2)
## Needs: snps.bam (bowtie2 --sensitive, MAPQ >= 10, sorted)
## The VCF ID column must equal the first field (the SNP name).  
    # CHECK head(snp$id)
fna <- in_files$fna
if (file.exists(in_files$bam)) {
  bam <- scanBam(in_files$bam, param = ScanBamParam(
    what = c("qname", "flag", "rname", "pos", "mapq", "cigar", "qwidth"),
    flag = scanBamFlag(isUnmappedQuery = FALSE)))[[1]]
  aln <- data.frame(qname = bam$qname, flag = bam$flag, chr = as.character(bam$rname),
                    pos0 = bam$pos, mapq = bam$mapq, cigar = bam$cigar,
                    qlen = bam$qwidth, stringsAsFactors = FALSE)
  aln <- aln[bitwAnd(aln$flag, 256L) == 0L & aln$mapq >= 10, ]   # primary, MAPQ >= 10
  aln$strand <- ifelse(bitwAnd(aln$flag, 16L) > 0L, "-", "+")
  hp <- strsplit(aln$qname, "|", fixed = TRUE)
  aln$id    <- vapply(hp, `[`, "", 1)
  aln$flank <- as.integer(vapply(hp, `[`, "", 2))
  aln$a1    <- vapply(hp, `[`, "", 3)
  aln$a2    <- vapply(hp, `[`, "", 4)
  aln$offset <- aln$flank + 1L                                   # 1-based SNP position in the read
  aln$qpos   <- ifelse(aln$strand == "+", aln$offset, aln$qlen - aln$offset + 1L)

  ## query position -> reference position through the CIGAR string
  q2r <- function(pos, cigar, qpos) {
    ops <- regmatches(cigar, gregexpr("[0-9]+[MIDNSHP=X]", cigar))[[1]]
    len <- as.integer(sub("[A-Z=]$", "", ops)); op <- sub("^[0-9]+", "", ops)
    q <- 0L; r <- pos - 1L
    for (i in seq_along(op)) {
      if (op[i] %in% c("M", "=", "X")) {
        if (qpos <= q + len[i]) return(r + (qpos - q))
        q <- q + len[i]; r <- r + len[i]
      } else if (op[i] %in% c("I", "S")) {
        if (qpos <= q + len[i]) return(NA_integer_)
        q <- q + len[i]
      } else if (op[i] %in% c("D", "N")) r <- r + len[i]
    }
    NA_integer_
  }
  aln$rpos <- mapply(q2r, aln$pos0, aln$cigar, aln$qpos, USE.NAMES = FALSE)

  ## reference base at each SNP position (FASTA index; genome is not loaded)
  if (!file.exists(paste0(fna, ".fai"))) indexFa(fna)
  fa <- FaFile(fna)
  ok <- !is.na(aln$rpos)
  aln$ref_base <- NA_character_
  aln$ref_base[ok] <- as.character(getSeq(fa, GRanges(aln$chr[ok], IRanges(aln$rpos[ok], width = 1))))
  comp <- c(A = "T", C = "G", G = "C", T = "A")
  b1 <- ifelse(aln$strand == "+", aln$a1, unname(comp[aln$a1]))   # alleles on the reference strand
  b2 <- ifelse(aln$strand == "+", aln$a2, unname(comp[aln$a2]))
  aln$ref_ok <- aln$ref_base == b1 | aln$ref_base == b2

  ## assembled chromosomes = headers that say "chromosome" and not scaffold/unlocalized
  hdr <- system(paste("grep '>'", shQuote(fna)), intern = TRUE)
  keep_h <- grepl("chromosome", hdr, ignore.case = TRUE) &
            !grepl("scaffold|unlocalized|unplaced", hdr, ignore.case = TRUE)
  chr_assembled <- sub("^>([^ ]+) .*$", "\\1", hdr[keep_h])          

  remap <- aln[!is.na(aln$rpos) & aln$ref_ok %in% TRUE, c("id", "chr", "rpos")]
  names(remap) <- c("id", "chr", "pos"); remap <- remap[!duplicated(remap$id), ]
  m <- match(snp$id, remap$id); snp$chr <- remap$chr[m]; snp$pos <- remap$pos[m]

  ## Numbers for the Methods: mutually exclusive groups
  remap_counts <- c(vcf_snps = nrow(snp), mapped = sum(snp$id %in% remap$id),
    no_alignment_or_MAPQ_below_10 = sum(!snp$id %in% aln$id),
    reference_base_mismatch_or_indel = sum(snp$id %in% aln$id & !snp$id %in% remap$id),
    aligned_to_unplaced_scaffold = sum(!is.na(snp$chr) & !snp$chr %in% chr_assembled))
  print(remap_counts)     # the last group is a subset of 'mapped'; it is excluded below
} else warning(in_files$bam, " not found: using the positions already in the VCF")

## Inspect chromosome names; restrict to assembled chromosomes if needed
#table(snp$chr, useNA = "ifany")                       
if (!exists("chr_assembled")) chr_assembled <- NULL   
snp$mapped <- !is.na(snp$chr) & !snp$chr %in% c(".", "0", "") &
              !is.na(snp$pos) & snp$pos > 0 &
              (is.null(chr_assembled) | snp$chr %in% chr_assembled)

## SNP annotation: Carrier et al. 2022 exclusions, X-linked, biallelic, duplicates
no_call <- read.table(in_files$no_call) 
mono <- read.table(in_files$mono)
x_chr   <- read.table(in_files$x_probes)
snp$carrier_drop <- snp$id %in% c(no_call$V1, mono$V1)
c(no_call = sum(snp$id %in% no_call$V1), monomorphic = sum(snp$id %in% mono$V1),
  removed_total = sum(snp$carrier_drop))        # numbers for the Methods; compare with the text
snp$is_X <- snp$id %in% x_chr$V1

x_chr_name <- names(which.max(table(snp$chr[snp$is_X & !is.na(snp$chr)])))
snp$is_X   <- snp$is_X | (!is.na(snp$chr) & snp$chr == x_chr_name)
snp$biallelic <- !grepl(",", snp$alt)
key <- paste(snp$chr, snp$pos)
snp$dup <- duplicated(snp$id) | duplicated(snp$id, fromLast = TRUE) |
           duplicated(key)    | duplicated(key,    fromLast = TRUE)

## Genotype matrix 
gt <- extract.gt(vcf, element = "GT")
gt <- gsub("\\|", "/", gt)
G  <- matrix(NA_integer_, nrow(gt), ncol(gt), dimnames = list(NULL, colnames(gt)))
G[which(gt == "0/0")]             <- 2L
G[which(gt %in% c("0/1", "1/0"))] <- 1L
G[which(gt == "1/1")]             <- 0L

## Sort by chromosome and position (needed for LD pruning)
snp$chr_int <- match(snp$chr, unique(snp$chr[snp$mapped]))
ord <- order(snp$chr_int, snp$pos, na.last = TRUE)
snp <- snp[ord, ]; G <- G[ord, ]
snp$sid <- seq_len(nrow(snp))

## GDS for SNPRelate (LD pruning, KING). 
Ggds <- G; Ggds[is.na(Ggds)] <- 3L
gds_fn <- file.path(dir_derived, "caribou_all.gds")
if (file.exists(gds_fn)) file.remove(gds_fn)
snpgdsCreateGeno(gds_fn, genmat = Ggds, sample.id = colnames(G), snp.id = snp$sid,
                 snp.chromosome = ifelse(is.na(snp$chr_int), 0L, snp$chr_int),
                 snp.position   = ifelse(is.na(snp$pos), 0L, snp$pos),
                 snp.allele     = paste(snp$ref, snp$alt, sep = "/"),
                 snpfirstdim = TRUE)
genofile <- snpgdsOpen(gds_fn)

## SNP filtering function 
build_snpset <- function(sel, miss_max, mac_min, hwe_p = NULL,
                         ld = TRUE, ld_r2 = 0.5) {
  L <- list()
  rec <- function(step, ids) L[[length(L) + 1]] <<- data.frame(step = step, n_snps = length(ids))

  ids <- snp$sid;                                                rec("1 Variants in VCF", ids)
  ids <- ids[!snp$carrier_drop[ids]];                            rec("2 After Carrier et al. 2022 exclusions", ids)
  ids <- ids[snp$mapped[ids] & snp$biallelic[ids] & !snp$dup[ids]]
                                                                 rec("3 Mapped, biallelic, unique ID/position", ids)
  ids <- ids[!snp$is_X[ids]];                                    rec("4 Autosomal (X-linked set aside)", ids)

  Gs <- G[ids, sel, drop = FALSE]

  ok <- rowMeans(is.na(Gs)) <= miss_max
  ids <- ids[ok]; Gs <- Gs[ok, , drop = FALSE];                  rec("5 Per-SNP missingness", ids)

  n_ref <- rowSums(Gs, na.rm = TRUE); n_all <- 2 * rowSums(!is.na(Gs))
  ok <- pmin(n_ref, n_all - n_ref) >= mac_min                    # also drops monomorphic
  ids <- ids[ok]; Gs <- Gs[ok, , drop = FALSE];                  rec("6 Minor allele count", ids)

  if (!is.null(hwe_p)) {
    cnt <- cbind(AA = rowSums(Gs == 2, na.rm = TRUE),
                 AB = rowSums(Gs == 1, na.rm = TRUE),
                 BB = rowSums(Gs == 0, na.rm = TRUE))
    p <- HWExactMat(cnt, verbose = FALSE)$pvalvec                
    ok <- p >= hwe_p
    ids <- ids[ok]; Gs <- Gs[ok, , drop = FALSE]
  };                                                             rec("7 HWE exact test (not used by default)", ids)

  if (ld) {
    ## SNPRelate thresholds on |r|, so convert r^2 -> r. 
    pr  <- snpgdsLDpruning(genofile, sample.id = sel, snp.id = ids,
                           autosome.only = FALSE, slide.max.n = P$ld_win,
                           slide.max.bp = 1e8, ld.threshold = sqrt(ld_r2),
                           method = "corr", start.pos = "first",   # fixed start = reproducible
                           verbose = FALSE)
    ids <- sort(unlist(pr, use.names = FALSE))
  };                                                             rec("8 LD pruning", ids)

  list(ids = ids, log = do.call(rbind, L))
}

## Filter configurations: main + sensitivity. 
cfg <- list(
  Structure = list(miss_max = P$miss_max, mac_min = P$mac_min, hwe_p = NULL, ld = TRUE,  ld_r2 = P$ld_r2),
  Kinship   = list(miss_max = P$miss_max, mac_min = 1,         hwe_p = NULL, ld = TRUE,  ld_r2 = P$ld_r2),
  Lenient   = list(miss_max = 0.10,       mac_min = 2,         hwe_p = NULL, ld = FALSE),
  StrictLD  = list(miss_max = P$miss_max, mac_min = P$mac_min, hwe_p = NULL, ld = TRUE,  ld_r2 = 0.2)
)
snpsets <- lapply(cfg, function(a) do.call(build_snpset, c(list(sel = datasets$DS2), a)))

## Reproducibility table for Methods/SI (counts after each step, per configuration)
filter_log <- cbind(Step = snpsets[[1]]$log$step,
                    sapply(snpsets, function(s) s$log$n_snps))
params <- t(sapply(cfg, function(a)
  sapply(c("miss_max", "mac_min", "hwe_p", "ld", "ld_r2"),
         function(k) if (is.null(a[[k]])) NA else a[[k]])))
filter_log; params

## Per-sample missingness on the final (Structure) SNP set
ids_struct <- snpsets$Structure$ids
miss_tab <- bind_rows(lapply(names(datasets), function(d)
  data.frame(dataset = d, sample = datasets[[d]],
             missing = colMeans(is.na(G[ids_struct, datasets[[d]], drop = FALSE])))))
miss_tab %>% group_by(dataset) %>%
  summarise(n = n(), median_miss = median(missing), n_over_50 = sum(missing > 0.5))
write.csv(miss_tab, file.path(out, "per_sample_missingness.csv"), row.names = FALSE)

## Save for downstream scripts
snpgdsClose(genofile)
save_obj(vcf, "vcf")
save_obj(list(snp = snp, G = G), "snp_data")
save_obj(lapply(snpsets, `[[`, "ids"), "snpsets_ids")

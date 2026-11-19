
source("1_Libraries_to_load.R")

## Load sample table
samples <- read.csv(in_files$samples)

smp <- samples %>%
  transmute(Sample     = gsub("\\s+", "", Sample.ID),
            Individual = sub("_.*$", "", gsub("\\s+", "", Sample.ID)),
            Type       = Sample_type,       # "Blood" / "Feces"
            CallRate   = Call.Rate)

## Full metadata table used by all downstream scripts
info <- samples %>%
  transmute(Sample        = gsub("\\s+", "", Sample.ID),
            Individual    = sub("_.*$", "", Sample),
            Institution   = sub("^[^_]+_([^_]+)_.*$", "\\1", Sample),
            Type          = Sample_type,                  # "Blood" / "Feces"
            Sex           = tolower(Sex),                 # "male" / "female"
            CallRate      = Call.Rate,
            Ecotype       = Ecotype,                      # as subspecies recorded in ZIMS
            ZIMS_reindeer = Ecotype %in% c("Reindeer", "Eurasian tundra reindeer"),
            Category      = cut(CallRate, c(-Inf, P$cr_fail, P$cr_high, Inf),
                                labels = c("Failed", "Lower confidence", "High quality"),
                                right = FALSE))

## Lower-confidence samples (call rate 0.15-0.90) are flagged with * in every figure
lowconf <- info$Sample[info$Category == "Lower confidence"]

## Define sample sets
## Drop failed samples; within individual x sample type keep the sample with the
## higher call rate (resolves the transfer duplicates DKJ23-05769 / -05773).
smp1 <- smp %>%
  filter(CallRate >= P$cr_fail) %>%
  group_by(Individual, Type) %>%
  slice_max(CallRate, n = 1, with_ties = FALSE) %>%
  ungroup()

## One sample per individual; for paired individuals choose which type to keep
pick <- function(smp, type) {
  smp <- as.data.frame(smp)
  smp <- smp[order(-smp$CallRate), ]            # best call rate first
  has_type <- as.logical(ave(smp$Type == type, smp$Individual, FUN = any))
  keep <- smp[ifelse(has_type, smp$Type == type, TRUE), ]
  keep <- keep[!duplicated(keep$Individual), ]  # one row per individual
  keep$Sample
}
ds1_blood <- pick(smp1, "Blood")
ds1_feces <- pick(smp1, "Feces")

## Dataset 2: high quality only, one sample per individual (highest call rate)
ds2 <- smp %>%
  filter(CallRate >= P$cr_high) %>%
  group_by(Individual) %>%
  slice_max(CallRate, n = 1, with_ties = FALSE) %>%
  ungroup() %>% pull(Sample)

datasets <- list(DS1_blood = ds1_blood, DS1_feces = ds1_feces, DS2 = ds2)
sapply(datasets, length)
stopifnot(all(sapply(datasets, length) == P$n_expected[names(datasets)]))

## Individuals with both a blood and a feces sample (n = 9)
paired <- smp1 %>% group_by(Individual) %>%
  filter(all(c("Blood", "Feces") %in% Type)) %>% ungroup()
paired_ids <- unique(paired$Individual)


## Save for downstream scripts
save_obj(list(smp = smp, smp1 = smp1, info = info, lowconf = lowconf,
              datasets = datasets, paired = paired, paired_ids = paired_ids),
         "sample_sets")

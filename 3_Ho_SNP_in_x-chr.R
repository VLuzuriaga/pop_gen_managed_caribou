source("1_Libraries_to_load.R")
list2env(load_obj("sample_sets"), envir = .GlobalEnv)
list2env(load_obj("snp_data"),    envir = .GlobalEnv)   

## X-linked probes: not subject to the autosomal MAC / LD filters
x_ids <- snp$sid[snp$is_X & snp$mapped & !snp$carrier_drop & !snp$dup & snp$biallelic]

ho_x <- function(s) colMeans(G[x_ids, s, drop = FALSE] == 1, na.rm = TRUE)
ho_tab <- bind_rows(lapply(names(datasets), function(d)
  data.frame(dataset = d, sample = datasets[[d]], Ho = ho_x(datasets[[d]]))))
ho_tab$n_called <- sapply(ho_tab$sample, function(s) sum(!is.na(G[x_ids, s])))
save_obj(ho_tab, "ho_tab")

## Accuracy of the sex assignment

with(filter(info, Category != "Failed"), table(Type, Sex, Category))

sex_tab <- ho_tab %>%
  left_join(select(info, Sample, Individual, Type, Sex, CallRate, Category),
            by = c("sample" = "Sample")) %>%
  mutate(Predicted = ifelse(Ho < P$ho_x_male, "male", "female"),   
         Correct   = Predicted == Sex)

sex_summary <- sex_tab %>%
  group_by(dataset, Type, Sex) %>%
  summarise(n = n(), n_correct = sum(Correct, na.rm = TRUE),
            n_no_Ho = sum(is.na(Ho)), .groups = "drop")
sex_summary
filter(sex_tab, !Correct | is.na(Correct))             

## Figure: X-chromosome Ho vs missing genotypes; * = lower-confidence sample
sex_tab$Pair <- ifelse(sex_tab$dataset %in% c("DS1_blood", "DS1_feces") &
                         sex_tab$Individual %in% paired_ids,
                       paste(sex_tab$Type, "(paired)"), "Unpaired")

p_ho <- ggplot(sex_tab, aes(Ho, 100 * (1 - CallRate),
                            colour = Sex, shape = Type, fill = Pair)) +
  geom_point(size = 3, alpha = 0.9, stroke = 0.7) +
  geom_text(aes(label = ifelse(Category == "Lower confidence", "*", "")),
            colour = "black", size = 5, vjust = -0.2, show.legend = FALSE) +
  geom_vline(xintercept = P$ho_x_male, linetype = "dashed") +
  scale_colour_manual(values = c(female = "magenta3", male = "royalblue")) +
  scale_shape_manual(values = c(Blood = 21, Feces = 24)) +
  scale_fill_manual(
    name   = "Paired individuals",
    values = c("Blood (paired)" = "orange2",
               "Feces (paired)" = "mediumpurple4",
               "Unpaired"       = "white")) +
  guides(fill   = guide_legend(override.aes = list(shape = 21, colour = "grey30", size = 4)),
         colour = guide_legend(override.aes = list(shape = 21, fill = "white", size = 4))) +
  facet_wrap(~dataset) +
  labs(x = expression("X-chromosome "~italic(Ho)), y = "Missing genotypes (%)",
       caption = "* lower-confidence sample (call rate 0.15-0.90)") +
  theme_bw()
p_ho

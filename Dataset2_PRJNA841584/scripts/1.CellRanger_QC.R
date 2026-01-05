# Script ran on high performance Linux machine - access provided by Prof. Pieter De Maayer

# Load libraries
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(reshape2)
library(tibble)

# Names of samples (same name as folders stored in data)
samples <- c("SAMN28600744", "SAMN28600743", "SAMN28600742", "SAMN28600741", "SAMN28600740", "SAMN28600739", "SAMN28600738", "SAMN28600737", "SAMN28600736", "SAMN28600735", "SAMN28600734", "SAMN28600733", "SAMN28600732", "SAMN28600731", "SAMN28600730", "SAMN28600729", "SAMN28600728", "SAMN28600727", "SAMN28600726", "SAMN28600725", "SAMN28600724", "SAMN28600723", "SAMN28600722")

# Loop over each sample and read the metrics summary in
metrics <- list()
for (sample in samples) {
  path_csv <- paste0("/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA841584/CellRanger/", sample, "/outs/metrics_summary.csv")
  df <- read.csv(path_csv)
  rownames(df) <- sample
  metrics[[sample]] <- df
}

# Concatenate each sample metrics together
metrics <- ldply(metrics, rbind)

# Remove periods and percentages to make the values numeric
metrics <- metrics %>%
  column_to_rownames(".id") %>%
  mutate_all(funs(parse_number(str_replace(., ",", "")))) %>%
  mutate_all(funs(parse_number(str_replace(., "%", ""))))
metrics$sample <- rownames(metrics)

# Columns of interest
cols <- c("Reads.Mapped.Confidently.to.Intergenic.Regions",
          "Reads.Mapped.Confidently.to.Intronic.Regions",
          "Reads.Mapped.Confidently.to.Exonic.Regions",
          "sample")

# Data wrangling to sculpt dataframe in a ggplot friendly manner
df <- metrics %>%
  select(cols) %>%
  melt(id.vars = "sample") %>%  # Preserve 'sample' during melting
  mutate(variable = str_replace_all(variable, "Reads.Mapped.Confidently.to.", "")) %>%
  mutate(variable = str_replace_all(variable, ".Regions", ""))

# ggplot code to make a barplot
each_region_reads <- df %>% ggplot() +
  geom_bar(
    aes(x = sample, y = value, fill = variable),
    position = "stack",
    stat = "identity") +
  coord_flip() +
  labs(
    x = "Sample",
    y = "Percentage of Reads",
    title = "Percent of Reads Mapped to Each Region",
    fill = "Region")

# Save the plot
ggsave("/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/Data/PRJNA841584/QC/figures/CellRanger/each_region_reads.png", plot = each_region_reads)

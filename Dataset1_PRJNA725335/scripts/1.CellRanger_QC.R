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
samples <- c("SAMN18876198", "SAMN18876199", "SAMN18876200", "SAMN18876201", "SAMN18876202", "SAMN18876203", "SAMN18876204", "SAMN18876205", "SAMN18876206", "SAMN18876207", "SAMN18876208", "SAMN18876209", "SAMN18876210", "SAMN18876211", "SAMN18876212")

# Loop over each sample and read the metrics summary in
metrics <- list()
for (sample in samples) {
  path_csv <- paste0("/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA725335/CellRanger/", sample, "/outs/metrics_summary.csv")
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
ggsave("/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/Data/PRJNA725335/QC/figures/CellRanger/each_region_reads.png", plot = each_region_reads)

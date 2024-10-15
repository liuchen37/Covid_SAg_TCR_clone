rm(list = ls())

library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)

# Import files
X1D0 <- read_excel("1D0.xlsx")
X1D6 <- read_excel("1D6.xlsx")
X2D0 <- read_excel("2D0.xlsx")
X2D6 <- read_excel("2D6.xlsx")
X3D0 <- read_excel("3D0.xlsx")
X3D6 <- read_excel("3D6.xlsx")

# Combine all the datasets
all_data <- rbind(
  X1D0 %>% mutate(Sample = "X1D0"),
  X1D6 %>% mutate(Sample = "X1D6"),
  X2D0 %>% mutate(Sample = "X2D0"),
  X2D6 %>% mutate(Sample = "X2D6"),
  X3D0 %>% mutate(Sample = "X3D0"),
  X3D6 %>% mutate(Sample = "X3D6")
)

# Keep the entry with the higher fraction value for each allVHitsWithScore
all_data_filtered <- all_data %>%
  group_by(Sample, allVHitsWithScore) %>%
  summarise(cloneFraction = max(cloneFraction))

# Create a new table with allVHitsWithScores and cloneFractions
new_table <- all_data_filtered %>%
  pivot_wider(names_from = Sample, values_from = cloneFraction) %>%
  drop_na()

#Split the fraction values based on their group
new_table_split <- new_table %>%
  mutate(
    X1_Resting = X1D0,
    X2_Resting = X2D0,
    X3_Resting = X3D0,
    X1_Peptide = X1D6,
    X2_Peptide = X2D6,
    X3_Peptide = X3D6
  ) %>%
  select(allVHitsWithScore, X1_Resting, X2_Resting, 
         X3_Resting, X1_Peptide, X2_Peptide, X3_Peptide)

# Calculate mean and standard deviation for Resting and Peptide groups
summary_table <- new_table_split %>%
  rowwise() %>%
  mutate(
    Resting_Mean = mean(c(X1_Resting, X2_Resting, X3_Resting)),
    Resting_SD = sd(c(X1_Resting, X2_Resting, X3_Resting)),
    Peptide_Mean = mean(c(X1_Peptide, X2_Peptide, X3_Peptide)),
    Peptide_SD = sd(c(X1_Peptide, X2_Peptide, X3_Peptide)),
    Mean_Diff = Peptide_Mean - Resting_Mean
  )

# Select the top 10 rows based on the absolute value of Mean_Diff
top_10_table <- summary_table %>%
  arrange(desc(abs(Mean_Diff))) %>%
  head(10)

# Convert the top_10_table to a long format for plotting
plot_data <- top_10_table %>%
  select(allVHitsWithScore, X1_Resting, X2_Resting, X3_Resting, 
         X1_Peptide, X2_Peptide, X3_Peptide) %>%
  pivot_longer(cols = -allVHitsWithScore, 
               names_to = c("Sample", "Group"),
               names_sep = "_", 
               values_to = "cloneFraction")

# Reorder the levels of the Group factor
plot_data$Group <- factor(plot_data$Group, levels = c("Resting", "Peptide"))

# Create the box plot with log10 scale
ggplot(plot_data, aes(x = allVHitsWithScore, y = cloneFraction, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.5) +
  stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.75), width = 0.2) +
  labs(x = " ", y = "Clonal Fraction", fill = "Group") +
  scale_y_log10(labels = scales::scientific) +
  scale_fill_manual(values = c("azure3", "deepskyblue1"), labels = c("CTL", "P3")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )

## Seperate figures by chain
# Split the data based on gene type
trav_data <- new_table_split[grepl("TRAV", new_table_split$allVHitsWithScore), ]
trbv_data <- new_table_split[grepl("TRBV", new_table_split$allVHitsWithScore), ]
trdv_data <- new_table_split[grepl("TRDV", new_table_split$allVHitsWithScore), ]

# Create a data frame with chain types and their counts
chain_counts <- data.frame(
  Chain = c("α", "β", "δ"),
  Count = c(603, 800, 8)
)

# Create the histogram
ggplot(chain_counts, aes(x = Chain, y = Count, fill = Chain)) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(x = "TCR Chain", y = "Number of Observations", title = "Events by Chains") +
  scale_fill_manual(values = c("chocolate1", "cyan3", "darkseagreen2")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        legend.position = "none")

# Calculate mean and standard deviation for Resting and Peptide groups based on A/B/D
summary_table_trav <- trav_data %>%
  rowwise() %>%
  mutate(
    Resting_Mean = mean(c(X1_Resting, X2_Resting, X3_Resting)),
    Resting_SD = sd(c(X1_Resting, X2_Resting, X3_Resting)),
    Peptide_Mean = mean(c(X1_Peptide, X2_Peptide, X3_Peptide)),
    Peptide_SD = sd(c(X1_Peptide, X2_Peptide, X3_Peptide)),
    Mean_Diff = Peptide_Mean - Resting_Mean
  )

summary_table_trbv <- trbv_data %>%
  rowwise() %>%
  mutate(
    Resting_Mean = mean(c(X1_Resting, X2_Resting, X3_Resting)),
    Resting_SD = sd(c(X1_Resting, X2_Resting, X3_Resting)),
    Peptide_Mean = mean(c(X1_Peptide, X2_Peptide, X3_Peptide)),
    Peptide_SD = sd(c(X1_Peptide, X2_Peptide, X3_Peptide)),
    Mean_Diff = Peptide_Mean - Resting_Mean
  )

# Select the top 10 (8 for delta) rows based on the absolute value of Mean_Diff
top_10_table_trav <- summary_table_trav %>%
  arrange(desc(abs(Mean_Diff))) %>%
  head(10)

top_10_table_trbv <- summary_table_trbv %>%
  arrange(desc(abs(Mean_Diff))) %>%
  head(10)

top_8_table_trdv <- summary_table_trdv %>%
  arrange(desc(abs(Mean_Diff))) %>%
  head(8)

# Convert the top_10(8)_tables to a long format for plotting
plot_data_trav <- top_10_table_trav %>%
  select(allVHitsWithScore, X1_Resting, X2_Resting, X3_Resting, 
         X1_Peptide, X2_Peptide, X3_Peptide) %>%
  pivot_longer(cols = -allVHitsWithScore, 
               names_to = c("Sample", "Group"),
               names_sep = "_", 
               values_to = "cloneFraction")

plot_data_trbv <- top_10_table_trbv %>%
  select(allVHitsWithScore, X1_Resting, X2_Resting, X3_Resting, 
         X1_Peptide, X2_Peptide, X3_Peptide) %>%
  pivot_longer(cols = -allVHitsWithScore, 
               names_to = c("Sample", "Group"),
               names_sep = "_", 
               values_to = "cloneFraction")

plot_data_trdv <- top_8_table_trdv %>%
  select(allVHitsWithScore, X1_Resting, X2_Resting, X3_Resting, 
         X1_Peptide, X2_Peptide, X3_Peptide) %>%
  pivot_longer(cols = -allVHitsWithScore, 
               names_to = c("Sample", "Group"),
               names_sep = "_", 
               values_to = "cloneFraction")

# Reorder the levels of the Group factor
plot_data_trav$Group <- factor(plot_data_trav$Group, levels = c("Resting", "Peptide"))
plot_data_trbv$Group <- factor(plot_data_trbv$Group, levels = c("Resting", "Peptide"))
plot_data_trdv$Group <- factor(plot_data_trdv$Group, levels = c("Resting", "Peptide"))

# Create the box plot with log10 scale
ggplot(plot_data_trav, aes(x = allVHitsWithScore, y = cloneFraction, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.5) +
  stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.75), width = 0.2) +
  labs(x = " ", y = "Clonal Fraction", fill = "Group") +
  scale_y_log10(labels = scales::scientific) +
  scale_fill_manual(values = c("azure3", "deepskyblue1"), labels = c("CTL", "P3")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )

ggplot(plot_data_trbv, aes(x = allVHitsWithScore, y = cloneFraction, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.5) +
  stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.75), width = 0.2) +
  labs(x = " ", y = "Clonal Fraction", fill = "Group") +
  scale_y_log10(labels = scales::scientific) +
  scale_fill_manual(values = c("azure3", "deepskyblue1"), labels = c("CTL", "P3")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )


# Statistices
# Reshape the data from wide to long format
trav_data_long <- trav_data %>%
  pivot_longer(cols = starts_with(c("X1", "X2", "X3")),
               names_to = c("Sample", "Condition"),
               names_sep = "_",
               values_to = "Value")

trbv_data_long <- trbv_data %>%
  pivot_longer(cols = starts_with(c("X1", "X2", "X3")),
               names_to = c("Sample", "Condition"),
               names_sep = "_",
               values_to = "Value")

# Summary
top_10_trav_summary <- top_10_table_trav %>%
  pivot_longer(cols = starts_with(c("X1", "X2", "X3")),
               names_to = c("Sample", "Condition"),
               names_sep = "_",
               values_to = "Value") %>%
  group_by(allVHitsWithScore) %>%
  summarise(
    Resting_Mean = mean(Value[Condition == "Resting"]),
    Resting_SD = sd(Value[Condition == "Resting"]),
    Peptide_Mean = mean(Value[Condition == "Peptide"]),
    Peptide_SD = sd(Value[Condition == "Peptide"]),
    t_value = t.test(Value[Condition == "Resting"], Value[Condition == "Peptide"])$statistic,
    p_value = t.test(Value[Condition == "Resting"], Value[Condition == "Peptide"])$p.value
  )

top_10_trbv_summary <- top_10_table_trbv %>%
  pivot_longer(cols = starts_with(c("X1", "X2", "X3")),
               names_to = c("Sample", "Condition"),
               names_sep = "_",
               values_to = "Value") %>%
  group_by(allVHitsWithScore) %>%
  summarise(
    Resting_Mean = mean(Value[Condition == "Resting"]),
    Resting_SD = sd(Value[Condition == "Resting"]),
    Peptide_Mean = mean(Value[Condition == "Peptide"]),
    Peptide_SD = sd(Value[Condition == "Peptide"]),
    t_value = t.test(Value[Condition == "Resting"], Value[Condition == "Peptide"])$statistic,
    p_value = t.test(Value[Condition == "Resting"], Value[Condition == "Peptide"])$p.value
  )

write.csv(top_10_trav_summary, "top_10_trav_summary.csv", row.names = FALSE)
write.csv(top_10_trbv_summary, "top_10_trbv_summary.csv", row.names = FALSE)

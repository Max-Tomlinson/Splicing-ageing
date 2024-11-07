library(ggplot2)
library(ggsci)
library(readxl)

theme_new <-  theme_classic(base_size = 24) +
  theme(text = element_text(family = "Helvetica"))
theme_set(new = theme_new)

pax_selected_samples <- read.table(file = "~/Documents/KCL/PhD/TwinsUK/Data/Nanopore/pax_selected_samples.txt", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
pax_selected_samples <- pax_selected_samples[!(pax_selected_samples$Study_No %in% c(371, 431, 1832, 2621, 32331, 371, 431, 32861, 51901)), ]

master_table <- read.table(file = "~/Documents/KCL/PhD/TwinsUK/Output/Master_table.txt", header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

bmi_data <- master_table[master_table$Description == "Body Mass Index", ]
bmi_data <- merge(x = pax_selected_samples, y = bmi_data, by = 1)
bmi_data$TimeDiff1 <- difftime(as.Date(bmi_data$VisitDate), as.Date(bmi_data$SampleDate1, format = "%d/%m/%Y"), units = "days")
bmi_data$TimeDiff2 <- difftime(as.Date(bmi_data$VisitDate), as.Date(bmi_data$SampleDate2, format = "%d/%m/%Y"), units = "days")

BMI_1 <- sapply(unique(x = bmi_data$Study_No), function(i) {
  id <- bmi_data[bmi_data$Study_No == i, ]
  BMI_date1 <- id[which.min(abs(id$TimeDiff1)), ]$Result
  return(BMI_date1)
})

BMI_2 <- sapply(unique(x = bmi_data$Study_No), function(i) {
  id <- bmi_data[bmi_data$Study_No == i, ]
  BMI_date2 <- id[which.min(abs(id$TimeDiff2)), ]$Result
  return(BMI_date2)
})

pax_selected_samples <- cbind(pax_selected_samples[order(pax_selected_samples$Study_No), ], BMI_1, BMI_2)
pax_selected_samples$Group <- ifelse(test = pax_selected_samples$Age1 < 60, yes = "Middle-aged", no = "Older-aged")

middle_aged <- pax_selected_samples[pax_selected_samples$Age1 < 60, ]
older_aged <- pax_selected_samples[pax_selected_samples$Age1 >= 60, ]

n <- nrow(x = pax_selected_samples)
SampleDate1_range <- paste0(min(as.Date(pax_selected_samples$SampleDate1, format = "%d/%m/%Y")), "-", max(as.Date(pax_selected_samples$SampleDate1, format = "%d/%m/%Y")))
SampleDate2_range <- paste0(min(as.Date(pax_selected_samples$SampleDate2, format = "%d/%m/%Y")), "-", max(as.Date(pax_selected_samples$SampleDate2, format = "%d/%m/%Y")))
male <- sum(pax_selected_samples$Sex)
female <- n - male
mz <- sum(pax_selected_samples$Zygosity == "MZ")
dz <- sum(pax_selected_samples$Zygosity == "DZ")

study_characteristics <- data.frame(
  Statistic = c(
    n, male, female, mz, dz,
    paste0(round(x = mean(x = middle_aged$Age1), digits = 2), "??", round(x = sd(x = middle_aged$Age1), digits = 2)),
    paste0(round(x = mean(x = middle_aged$Age2), digits = 2), "??", round(x = sd(x = middle_aged$Age2), digits = 2)),
    paste0(round(x = mean(x = older_aged$Age1), digits = 2), "??", round(x = sd(x = older_aged$Age1), digits = 2)),
    paste0(round(x = mean(x = older_aged$Age2), digits = 2), "??", round(x = sd(x = older_aged$Age2), digits = 2)),
    paste0(round(x = mean(x = pax_selected_samples$Agediff), digits = 2), "??", round(x = sd(x = pax_selected_samples$Agediff), digits = 2)),
    paste0(round(x = mean(x = middle_aged$BMI_1), digits = 2), "??", round(x = sd(x = middle_aged$BMI_1), digits = 2)),
    paste0(round(x = mean(x = middle_aged$BMI_2), digits = 2), "??", round(x = sd(x = middle_aged$BMI_2), digits = 2)),
    paste0(round(x = mean(x = older_aged$BMI_1), digits = 2), "??", round(x = sd(x = older_aged$BMI_1), digits = 2)),
    paste0(round(x = mean(x = older_aged$BMI_2), digits = 2), "??", round(x = sd(x = older_aged$BMI_2), digits = 2)),
    paste0(round(x = mean(BMI_1 - BMI_2), digits = 2), "??", round(x = sd(x = BMI_1 - BMI_2), digits = 2)),
    round(t.test(x = pax_selected_samples$BMI_1, y = pax_selected_samples$BMI_2)$p.value, digits = 2),
    SampleDate1_range, SampleDate2_range
  ),
  row.names = c(
    "n", "Male", "Female", "MZ", "DZ",
    "Middle-aged age 1 mean", "Middle-aged age 2 mean", "Older-aged age 1 mean", "Older-aged age 2 mean",
    "Time difference mean", "Middle-aged BMI 1 mean", "Middle-aged BMI 2 mean",
    "Older-aged BMI 1 mean", "Older-aged BMI 2 mean", "BMI difference mean", "BMI difference p-value",
    "Sample 1 date range", "Sample 2 date range"
  )
)

write.table(x = pax_selected_samples, file = "~/Documents/KCL/PhD/TwinsUK/Data/Nanopore/pax_selected_samples_processed.txt", quote = FALSE, sep = "\t")
write.csv(x = study_characteristics, file = "~/Desktop/Nanopore/Analysis/Study characteristics/Study characteristics.csv")

plot_boxplot <- function(data, x_col, group_col, file_name) {
  
  boxplot_plot <- ggplot(data = data, mapping = aes(x = .data[[group_col]], y = .data[[x_col]], fill = .data[[group_col]])) +
    geom_boxplot() +
    geom_point(position = position_dodge(width = 1)) +
    labs(x = NULL, y = gsub(pattern = "_", replacement = " ", x = x_col)) +
    scale_fill_futurama() +
    theme(legend.position = "none")
  
  ggsave(filename = file_name, plot = boxplot_plot, width = 11, height = 8.5, dpi = 300)
}

data <- data.frame(
  Study_No = c(pax_selected_samples$Study_No, pax_selected_samples$Study_No),
  BMI = c(pax_selected_samples$BMI_1, pax_selected_samples$BMI_2),
  Time = factor(x = c(rep(x = "Time 1", 18), rep(x = "Time 2", 18)))
)

plot_boxplot(data = data, x_col = "BMI", group_col = "Time", file_name = "~/Desktop/Nanopore/Analysis/Study characteristics/BMI.pdf")

RNA_extraction_data <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Nanopore/Extraction data/RNA_extraction_data_updated1[1].xlsx", skip = 1)

data <- read.csv(file = "~/Desktop/Nanopore/sample_info.csv")
data <- merge(x = data, y = RNA_extraction_data[, c("...1", "RIN...13")], by.x = "No", by.y = "...1", sort = FALSE)
data <- merge(x = pax_selected_samples, y = data[, c("Study_No", "Time", "RIN...13")], by = "Study_No", all.x = TRUE)

data <- data.frame(
  Study_No = data$Study_No,
  RIN = data$RIN...13,
  Time = factor(x = paste("Time", data$Time))
)

plot_boxplot(data = data, x_col = "RIN", group_col = "Time", file_name = "~/Desktop/Nanopore/Analysis/Study characteristics/RIN.pdf")

data <- data.frame(
  Study_No = pax_selected_samples$Study_No,
  Age = c(pax_selected_samples$Age1, pax_selected_samples$Age2),
  Group = paste(c(pax_selected_samples$Group, pax_selected_samples$Group), "Time", factor(rep(1:2, each = nrow(pax_selected_samples))))
)

plot <- ggplot(data = data, mapping = aes(x = Age, fill = Group)) +
    geom_density(adjust = 2) +
    labs(y = "Density") +
    scale_fill_futurama() +
    theme(legend.position = "none")

ggsave(filename = "~/Desktop/Nanopore/Analysis/Study characteristics/Age combined.pdf", plot = plot, width = 11, height = 8.5, dpi = 300)

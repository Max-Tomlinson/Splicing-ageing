library(cowplot)
library(dplyr)
library(ggplot2)
library(ggsci)
library(tidyr)

load(file = "~/Documents/KCL/PhD/TwinsUK/Output/Master_table_appended.RData")

master2 <- read.table(file = "~/Documents/KCL/PhD/TwinsUK/Output/Master_table.txt", header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
master2 <- split(x = master2, f = master2$Description)

master <- c(master, master2)

RI <- read.table(file = "~/Desktop/Nanopore/Analysis/glm/bambu_RI_fitted.txt", header = TRUE, sep = '\t', row.names = 1)
RI <- sqrt(x = apply(X = RI[, 19:36] - RI[, 1:18], MARGIN = 2, FUN = function(x) IQR(x = x, na.rm = TRUE)))

laplace <- read.table(file = "~/Desktop/Nanopore/Analysis/glm/laplace_entropy_fitted.txt", header = TRUE, sep = '\t', row.names = 1)
laplace <- sqrt(x = apply(X = laplace[, 19:36] - laplace[, 1:18], MARGIN = 2, FUN = function(x) IQR(x = x, na.rm = TRUE)))

gini <- read.table(file = "~/Desktop/Nanopore/Analysis/glm/gini_index_fitted.txt", header = TRUE, sep = '\t', row.names = 1)
gini <- sqrt(x = apply(X = gini[, 19:36] - gini[, 1:18], MARGIN = 2, FUN = function(x) IQR(x = x, na.rm = TRUE)))

extraction_data <- readxl::read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Nanopore/Extraction data/RNA_extraction_data_updated1[1].xlsx", skip = 1)[-c(12, 14, 32, 34), c("RIN...13")]

covar <- read.table(file = "~/Desktop/Nanopore/Analysis/Differential expression/dCell.csv", header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
covar <- cbind.data.frame(covar[1:18, ], Group = factor(x = c(rep("Middle_aged", 8), rep("Older_aged", 10))), RIN = extraction_data[1:18, ]$RIN...13, RI = RI, Laplace = laplace, Gini = gini)

results <- list()

for (i in 1:length(x = master)){
  
  trait <- master[[i]]
  trait <- trait[which(x = complete.cases(trait)), ]
  
  data <- merge(x = covar, y = trait, by.x = "Study_No", by.y = "ParticipantID", sort = FALSE)
  data <- data %>%
    arrange(Study_No, VisitDate) %>%
    distinct(Study_No, VisitDate, .keep_all = TRUE) %>%
    group_by(Study_No) %>%
    filter(row_number() == 1 | row_number() == n()) %>%
    mutate(Index = if_else(condition = row_number() == 1, true = "1", false = "2")) %>%
    pivot_wider(names_from = Index, values_from = c(VisitDate, Result), names_sep = "_")
  data <- data[which(x = complete.cases(data)), ]
  
  if ("VisitDate_2" %in% names(x = data)){
    if (!is.character(x = data$Result_1) && nrow(x = data) >= 9 && var(x = data$Result_1) != 0){
      
      data$diffres  <- data$Result_2 - data$Result_1
      
      if (var(x = data$diffres) != 0) {
        models <- list(
          model1 <- lm(formula = diffres ~ RI:Group, data = data),
          model2 <- lm(formula = diffres ~ Laplace:Group, data = data),
          model3 <- lm(formula = diffres ~ Gini:Group, data = data))
      }
      
      model_results <- data.frame(model = paste0("model", rep(x = 1:length(x = models))), 
                                  beta1 = numeric(length = length(x = models)), 
                                  beta2 = numeric(length = length(x = models)), 
                                  pval1 = numeric(length = length(x = models)), 
                                  pval2 = numeric(length = length(x = models)))
      
      for (i in seq_along(along.with = models)) {
        model <- summary(object = models[[i]])
        beta1 <- coef(object = model)[2, "Estimate"]
        beta2 <- coef(object = model)[3, "Estimate"]
        pval1 <- coef(object = model)[2, "Pr(>|t|)"]
        pval2 <- coef(object = model)[3, "Pr(>|t|)"]
        model_results$beta1[i] <- beta1
        model_results$beta2[i] <- beta2
        model_results$pval1[i] <- pval1
        model_results$pval2[i] <- pval2
      }
      
      model_results <- model_results %>% pivot_wider(names_from = model, values_from = c(beta1, beta2, pval1, pval2), names_glue = "{model}_{.value}")
        
      results[[trait$Description[1]]] <- cbind(n = nrow(x = data), model_results)
    }
  }
}

lm_molecular <- bind_rows(results, .id = "Phenotype")

write.csv(x = lm_molecular, file = "~/Desktop/Nanopore/Analysis/glm/lm_molecular.csv", quote = FALSE, row.names = FALSE)

theme_new <- theme_classic(base_size = 12) +
  theme(legend.position = "none", text = element_text(family = "Helvetica"))
theme_set(new = theme_new)

process <- function(trait, covar) {

  trait <- trait[which(x = complete.cases(trait)), ]
  
  data <- merge(x = covar, y = trait, by.x = "Study_No", by.y = "ParticipantID", sort = FALSE)
  data <- data %>%
    arrange(Study_No, VisitDate) %>%
    distinct(Study_No, VisitDate, .keep_all = TRUE) %>%
    group_by(Study_No) %>%
    filter(row_number() == 1 | row_number() == n()) %>%
    mutate(Index = if_else(condition = row_number() == 1, true = "1", false = "2")) %>%
    pivot_wider(names_from = Index, values_from = c(VisitDate, Result), names_sep = "_")
  data <- data[which(x = complete.cases(data)), ]
  data$diffres  <- data$Result_2 - data$Result_1

  return(data)
}

data1 <- process(trait = master$`Lean tissue in trunk region`, covar = covar)
data2 <- process(trait = master$`Percentage of fat in largest visceral fat region / Percentage of fat in android region`, covar = covar)
data3 <- process(trait = master$`Total sodium in serum`, covar = covar)
data4 <- process(trait = master$`Ventricular heart rate based on automatic measurements. It is the number of heart beats per minute measured at rest`, covar = covar )
data5 <- process(trait = master$`Lung function measured forced vital capacity`, covar = covar)
data6 <- process(trait = master$`Frailty index`, covar = covar)
data7 <- process(trait = master$`Bone mineral density for the thoracic region`, covar = covar)
data8 <- process(trait = master$`Bone mineral content for the thoracic region`, covar = covar)

plot1 <- ggplot(data = data1, mapping = aes(x = diffres, y = RI, colour = Group, fill = Group)) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_point(size = 1, shape = 20, stroke = 3) +
  geom_point(size = 3, shape = 1, stroke = 1, colour = "black") +
  labs(x = "Lean tissue in trunk region", y = "Intron retention") +
  scale_colour_futurama() +
  scale_fill_futurama()
  
plot2 <- ggplot(data = data2, mapping = aes(x = diffres, y = RI, colour = Group, fill = Group)) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_point(size = 1, shape = 20, stroke = 3) +
  geom_point(size = 3, shape = 1, stroke = 1, colour = "black") +
  labs(x = "Visceral fat: android fat ratio", y = "Intron retention") +
  scale_colour_futurama() +
  scale_fill_futurama()

plot3 <- ggplot(data = data3, mapping = aes(x = diffres, y = Laplace, colour = Group, fill = Group)) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_point(size = 1, shape = 20, stroke = 3) +
  geom_point(size = 3, shape = 1, stroke = 1, colour = "black") +
  labs(x = "Sodium", y = "Laplace entropy") +
  scale_colour_futurama() +
  scale_fill_futurama()

plot4 <- ggplot(data = data4, mapping = aes(x = diffres, y = Laplace, colour = Group, fill = Group)) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_point(size = 1, shape = 20, stroke = 3) +
  geom_point(size = 3, shape = 1, stroke = 1, colour = "black") +
  labs(x = "Ventricular heart rate", y = "Laplace entropy") +
  scale_colour_futurama() +
  scale_fill_futurama()

plot5 <- ggplot(data = data5, mapping = aes(x = diffres, y = Gini, colour = Group, fill = Group)) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_point(size = 1, shape = 20, stroke = 3) +
  geom_point(size = 3, shape = 1, stroke = 1, colour = "black") +
  labs(x = "Forced vital capacity", y = "Gini index") +
  scale_colour_futurama() +
  scale_fill_futurama()

plot6 <- ggplot(data = data6, mapping = aes(x = diffres, y = Gini, colour = Group, fill = Group)) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_point(size = 1, shape = 20, stroke = 3) +
  geom_point(size = 3, shape = 1, stroke = 1, colour = "black") +
  labs(x = "Frailty index", y = "Gini index") +
  scale_colour_futurama() +
  scale_fill_futurama()

plot7 <- ggplot(data = data7, mapping = aes(x = diffres, y = Laplace, colour = Group, fill = Group)) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_point(size = 1, shape = 20, stroke = 3) +
  geom_point(size = 3, shape = 1, stroke = 1, colour = "black") +
  labs(x = "Bone mineral density (thoracic region)", y = "Laplace entropy") +
  scale_colour_futurama() +
  scale_fill_futurama()

plot8 <- ggplot(data = data8, mapping = aes(x = diffres, y = Laplace, colour = Group, fill = Group)) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_point(size = 1, shape = 20, stroke = 3) +
  geom_point(size = 3, shape = 1, stroke = 1, colour = "black") +
  labs(x = "Bone mineral content (thoracic region)", y = "Laplace entropy") +
  scale_colour_futurama() +
  scale_fill_futurama()

plot1 <- plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, ncol = 2, labels = "AUTO", label_size = 8)
plot2 <- plot_grid(plot6, plot2, ncol = 2, labels = "AUTO", label_size = 8)

ggsave(filename = "~/Desktop/Nanopore/Analysis/glm/lm_molecular1.pdf", plot = plot1, width = 11, height = 8.5, useDingbats = FALSE)
ggsave(filename = "~/Desktop/Nanopore/Analysis/glm/lm_molecular2.pdf", plot = plot2, width = 11, height = 8.5, useDingbats = FALSE)
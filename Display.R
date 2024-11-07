#### Figure 1 ####

library(cowplot)
library(ggplot2)
library(ggsci)
library(magick)

theme_new <- theme_classic(base_size = 24) +
  theme(text = element_text(family = "Helvetica"))
theme_set(new = theme_new)

f1a <- ggdraw() + draw_image(image = image_read(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/Study design.png", density = 300))
f1b <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/Age combined.pdf"))
f1c <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/BMI.pdf"))
f1d <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/RIN.pdf"))

f1 <- plot_grid(f1a, plot_grid(f1b, f1c, f1d, labels = c("B", "C", "D"), nrow = 1, label_size = 8), nrow = 2, labels = "AUTO", label_size = 8)

ggsave(filename = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/Figures/1.pdf", plot = f1, width = 11, height = 8.5, useDingbats = FALSE)

#### Figure 2 ####

f2a <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/tmap.pdf"))
f2b <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/chr_gene.pdf"))
f2c <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/coding_probability.pdf"))
f2d <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/length.pdf"))
f2 <- plot_grid(f2a, f2b, f2c, f2d, ncol = 2, labels = "AUTO", label_size = 8)

ggsave(filename = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/Figures/2.pdf", plot = f2, width = 11, height = 8.5, useDingbats = FALSE)

#### Figure 3 ####

f3a <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/RPC.pdf"))
f3b <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/Pheatmap.pdf"))
f3c <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/Volcano_gene.pdf"))
f3d <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/Volcano_transcript.pdf"))
f3e <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/MultiMuTHER.pdf"))
f3f <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/gostplot.pdf"))
f3 <- plot_grid(f3a, f3b, f3c, f3d, f3e, f3f, ncol = 2, labels = "AUTO", label_size = 8)

ggsave(filename = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/Figures/3.pdf", plot = f3, width = 8.5, height = 11, useDingbats = FALSE)

#### Figure 4 ####

f4a <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/DEXSeq_DNAJA3.pdf"))
f4b <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/DEXSeq_MBNL3.pdf"))
f4c <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/PValues.pdf"))
f4d <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/Precision.pdf"))
f4e <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/Proportions2.pdf"))
f4f <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/Proportions3.pdf"))
f4 <- plot_grid(f4a, f4b, f4e, f4f, ncol = 2, labels = "AUTO", label_size = 8)

ggsave(filename = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/Figures/4.pdf", plot = f4, width = 11, height = 8.5, useDingbats = FALSE)

#### Figure 5 ####

f5a <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/SUPPA_events.pdf"))
f5b <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/DSE.pdf"))
f5c <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/bambu_isoform_volcano.pdf"))
f5d <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/Density_laplace.pdf"))
f5e <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/MA_entropy2.pdf"))
f5f <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/MA_gini2.pdf"))
f5g <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/MA_entropy_fitted.pdf"))
f5h <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/MA_gini_fitted.pdf"))

f5 <- plot_grid(f5a, f5b, f5c, f5d, f5g, f5h, ncol = 2, labels = "AUTO", label_size = 8)

ggsave(filename = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/Figures/5.pdf", plot = f5, width = 8.5, height = 11, useDingbats = FALSE)

#### Figure 6 ####

f6a <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/Figures/6/Volcano_bambu.pdf"))
f6b <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/Figures/6/DEXSeq2.pdf"))
f6c <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/Figures/6/Proportions.pdf"))
f6d <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/Figures/6/Enrichment.pdf"))
f6 <- plot_grid(f6a, f6b, f6c, f6d, ncol = 2, labels = "AUTO", label_size = 8)

ggsave(filename = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/Figures/6.pdf", plot = f6, width = 8.5, height = 11, useDingbats = FALSE)

#### Supplementary ####

fsa <- ggdraw() + draw_image(image = image_read_pdf(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/Figures/5/Supplementary/DEGs_summary.pdf"))
fsb <- ggdraw() + draw_image(image = image_read(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/Figures/5/Supplementary/gProfiler_DEXSeq.png"))
fsc <- ggdraw() + draw_image(image = image_read(path = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/Figures/5/Supplementary/gProfiler_DRIMSeq.png"))

fs <- plot_grid(fsa, fsb, fsc, ncol = 1, labels = "AUTO", label_size = 8)

ggsave(filename = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 3/Figures/5/Supplementary.pdf", plot = fs, width = 8.5, height = 11, useDingbats = FALSE)

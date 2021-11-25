library(data.table)
library(dplyr)
library(ggplot2)
library(argparse)
library(ggpubr)


#Command line
parser <- ArgumentParser()
parser <- ArgumentParser(description='This script is for generating PCA plots using results from acnestry PCA and assignment.')
parser$add_argument('-g1000', help='File with 1000G scores')
parser$add_argument('-g1000_info', help='File with 1000G pops')
parser$add_argument('-phenofile', help='Phenotype file used in QC')
parser$add_argument('-pca5', help='PCA file at prob 0.5')
parser$add_argument('-pca8', help='PCA file at prob 0.8')
args <- parser$parse_args()

# Running the script
# Rscript pca.R -g1000 1000G_scores.txt -g1000_info G1000_pops.txt \
# -phenofile cases_controls_phenotypes_ita_be_brazil_swe_ger_27_10_2020 \
# -pca5 cases_controls_pca_sup_pops_0.5_probs.txt \
# -pca8 cases_controls_pca_sup_pops_0.8_probs.txt

g1000 <- fread(args$g1000)
g1000_info <- fread(args$g1000_info)  %>%
  mutate(s=Sample)
g1000 <- g1000 %>%
  left_join(g1000_info, by="s")

cov_info <- fread(args$phenofile) %>%
  select(ID, Cohort, Analysis_C2) %>%
  mutate(s=ID) %>%
  distinct()

d5 <- fread(args$pca5)
d8 <- fread(args$pca8)

d5 <- d5 %>%
  left_join(cov_info, by = "s")
d8 <- d8 %>%
  left_join(cov_info, by = "s")

print("Table showing the population breakdown at p > 0.5")
table(d5$pop)
print("Table showing the population breakdown at p > 0.8")
table(d8$pop)

#d5 <- d5 %>%
  #mutate(c = case_when(is.na(Cohort) ~ "Sweden",
                       #TRUE ~ Cohort))
#d5$case <- d$Analysis_C2
#d5$case[which(is.na(d$case))] <- 0

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

d5_plt <- ggplot() +
  geom_point(data=g1000, aes(x=PC1, y=PC2, color=factor(SuperPop)), size = 1, alpha = 0.1) +
  geom_point(data=d5, aes(x=PC1, y=PC2, color=factor(pop)), size = 1) +
  scale_colour_manual(values=cbPalette) +
  ggtitle("PCA (p > 0.5)") +
  theme_minimal()

d5_plt_two <- ggplot() +
  geom_point(data=g1000, aes(x=PC2, y=PC3, color=factor(SuperPop)), size = 1, alpha = 0.1) +
  geom_point(data=d5, aes(x=PC2, y=PC3, color=factor(pop)), size = 1) +
  scale_colour_manual(values=cbPalette) +
  ggtitle("PCA (p > 0.5)") +
  theme_minimal()

d8_plt <- ggplot() +
  geom_point(data=g1000, aes(x=PC1, y=PC2, color=factor(SuperPop)), size = 1, alpha = 0.1) +
  geom_point(data=d8, aes(x=PC1, y=PC2, color=factor(pop)), size = 1) +
  scale_colour_manual(values=cbPalette) +
  ggtitle("PCA (p > 0.8)") +
  theme_minimal()

d8_plt_two <- ggplot() +
  geom_point(data=g1000, aes(x=PC2, y=PC3, color=factor(SuperPop)), size = 1, alpha = 0.1) +
  geom_point(data=d8, aes(x=PC2, y=PC3, color=factor(pop)), size = 1) +
  scale_colour_manual(values=cbPalette) +
  ggtitle("PCA (p > 0.8)") +
  theme_minimal()

arrange <- ggarrange(d5_plt, d8_plt, d5_plt_two, d8_plt_two, ncol = 2, nrow = 2)
ggsave("PCA_plots.png", arrange, width = 16, height = 12)
#arrange <- ggarrange(d5_plt_two, d8_plt_two, ncol = 2, nrow = 1)
#ggsave("PCA_plotsPC.png", arrange, width = 16, height = 6)





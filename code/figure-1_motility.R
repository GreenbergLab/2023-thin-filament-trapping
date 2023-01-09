library(readxl)
library(here)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggtext)
library(ggpubr)

# read in motility data and transform to long format got ggplot
motility_file <- here("data/motility/WT & Biotin Regulated Motility 200nM Tm&Tn.xlsx")
motility_data <- read_xlsx(motility_file, sheet = 4, skip = 1)
motility_data <- melt(setDT(motility_data),
                      measure.vars = c("regular-actin", "biotin-actin"),
                      variable.name = "id",
                     value.name = "velocity")

motility_data$id <- factor(motility_data$id, levels = c("regular-actin", "biotin-actin"))

motility_summary <- motility_data[, .(velocity_avg = mean(velocity, na.rm = TRUE),
                                      velocity_sd = sd(velocity, na.rm = TRUE)),
                                  by = id]


# check normality
shapiro.test(motility_data$velocity)

# it fails - use non-parametrics
t_test_results <- compare_means(velocity~id, data = motility_data, method = "wilcox.test")

colorz <- c("black", "grey50")

# Motility plot
(
gg_mot <-
ggplot(data = motility_data)+
  geom_boxplot(data = motility_data,
                aes(x = id,
                    y = velocity,
                     fill = id),
               outlier.color = NA,
               size = 1,
               alpha = 0.5,
                width = 0.3)+
  geom_jitter(data = motility_data,
              aes(x = id,
                  y = velocity,
                  color = id),
               width = 0.1,
              size =  1.5,
              alpha = 0.5,
              shape = 16)+
  stat_pvalue_manual(data = t_test_results, label = "p = {p.adj}", y.position=0.62, tip.length = 0, label.size = 7)+
  scale_y_continuous(expand = expansion(c(0, 0.1)))+
  scale_x_discrete(labels = c("regular-actin"="Actin", "biotin-actin"="Actin+Biotin"))+
  scale_color_manual(values = colorz)+
  scale_fill_manual(values = colorz)+
  coord_cartesian(ylim = c(0, NA))+
  xlab("")+
  ylab("Velocity (&#xb5;m/s)")+
  ggtitle("Regulated Motility pCa 4")+
  theme_cowplot(22)+
  theme(
    legend.position = "none",
    axis.title.y = element_markdown()
)
)

## ggsave(here("img", "fig-1_motility-boxplot.png"), dpi = 500, bg = "white")

# pca4 reg trace
gg_pca4 <- readRDS(here("img", "pca4-trace.rds"))
## saveRDS(gg_pca4, here("img/pca4-trace.rds"))
## ggsave(here("img/pca4-trace.png"), dpi = 500)

## after generating for figures, I imported to inkscape to combine with
## an existing .ai file I inherited for this project containing a pCa 9 trace.

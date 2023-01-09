library(readxl)
library(ggplot2)
library(cowplot)
library(tibble)
library(dplyr)
library(tidyr)
library(purrr)
library(ggtext)
library(ggpubr)
library(here)
library(magick)
library(pdftools)

textsize <- 9
basesize <- 11


# colors to use for plots
colz <- unname(palette.colors())[c(1, 4)]

## kslow of ATP binding
## data is from the data/stopped-flow/stopped-flow-data.xlsx sheet 1
kslow_reg <- tribble(
~"conc", ~"kslow", ~"id",
54.34, 22.04, "reg",
104.2275, 70.5, "reg",
268.0275, 67.58, "reg",
555.1, 66.36, "reg",
1067.625, 101.97, "reg",
2092.35, 84.16, "reg",
3646.5, 108.29, "reg",
6015.75, 110.49, "reg",
46.7025, 41.34, "reg",
91, 39.03, "reg",
234.91, 48.19, "reg",
476.775, 48.99, "reg",
946.075, 101.69, "reg",
2087.475, 89.71, "reg",
3357.25, 99.23, "reg",
5414.5, 99.75, "reg"
)

kslow_unreg <- tribble(
~"conc", ~"kslow", ~"id",
59.7025, 25.35, "unreg",
112.5475, 62.72, "unreg",
278.07, 49.14, "unreg",
610.675, 52.56, "unreg",
1240.525, 68.17, "unreg",
2226.9, 96.19, "unreg",
4042.5, 81.74, "unreg",
5885.75, 81.99, "unreg",
90.025, 51.0111, "unreg",
207.675, 52.69477, "unreg",
367.25, 54.3501, "unreg",
874.25, 71.0193, "unreg",
1472.25, 67.81424, "unreg",
2288, 74.80627, "unreg",
2288, 70.6904, "unreg",
4748.25, 62.534, "unreg"
)



kslow_df <- rbind(kslow_reg, kslow_unreg)
kslow_df$id <- factor(kslow_df$id, levels = c("unreg", "reg"))

kslow_nest <- kslow_df |>
  group_by(id) |>
  nest() |>
  mutate(mod = map(data, ~nls(kslow ~ (k2*conc) / (k1+conc),
                              data = .x,
                              start = list(k1 = 150, k2 = 100))),
         pred_df = map(mod, ~tibble(x = 1:6000,
                                    y = predict(.x, newdata = data.frame(conc = 1:6000)))))

kslow_predict_df <- kslow_nest |> select(id, pred_df) |> unnest(cols = c(pred_df))


(gg_kslow <-
ggplot()+
  geom_point(data = kslow_df,
             aes(x = conc,
                 y = kslow,
                 color = id))+
  geom_line(data = kslow_predict_df,
            aes(x = x,
                y = y,
                color = id))+
##   coord_cartesian(ylim = c(0, 1205))+
##   ggtitle("ATP Binding")+
     scale_color_manual(name = "", labels = c("Unregulated", "Regulated"), values = unname(palette.colors()[c(1, 4)]))+
##   xlab("Time (s)")+
   ggtitle("k<sub>+&#945;</sub>")+
   ylab("k<sub>slow</sub> (s<sup>-1</sup>)")+
  xlab("[ATP] &#40;&#xb5;M&#41;")+
  theme_cowplot(basesize)+
  theme(
    legend.position = "none",
     plot.title = element_markdown(),
     axis.title.y = element_markdown(),
     axis.title.x = element_markdown())
)



## calculate Afast:Aslow

ratio_df <- map2_df(c("unreg-ratio-amplitudes", "reg-ratio-amplitudes"),
                     c("unreg", "reg"),
                     ~read_xlsx(here("data/stopped-flow/stopped-flow-data.xlsx"), sheet = .x) |> mutate(id = .y)
            ) |>
  mutate(ratioA = a1/a2)

ratio_df$id <- factor(ratio_df$id, levels = c("unreg", "reg"))

# this value is K_alpha (capital K)
ratio_plataeu <- ratio_df |> filter(atp >= 2000) |> group_by(id) |> summarize(avg = mean(ratioA, na.rm = TRUE))







# from https://stats.stackexchange.com/questions/30394/how-to-perform-two-sample-t-tests-in-r-by-inputting-sample-statistics-rather-tha
# m1, m2: the sample means
# s1, s2: the sample standard deviations
# n1, n2: the same sizes
# m0: the null value for the difference in means to be tested for. Default is 0.
# equal.variance: whether or not to assume equal variance. Default is FALSE.
t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
    if( equal.variance==FALSE )
    {
        se <- sqrt( (s1^2/n1) + (s2^2/n2) )
        # welch-satterthwaite df
        df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
    } else
    {
        # pooled standard deviation, scaled by the sample sizes
        se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) )
        df <- n1+n2-2
    }
    t <- (m1-m2-m0)/se
    dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))
    names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
    return(dat)
}

# ttest from summary stats
adp_unreg_label <- paste0(adp_unreg_mean, "%+-%", round(adp_unreg_sd, 2))
adp_reg_label <- paste0(adp_reg_mean, "%+-%", adp_reg_sd)

# calcualte second order binding rate and propogate error
# k_{+2} / (1/K_1)
## unreg_2nd_order <- 1220/305.8

# k[+alpha] / K[+alpha]
unreg_k_minus_alpha <- 73.2/8.4

#                       sqrt( (((1/B)^2)*sd(A)^2) + (((A/B^2)^2)*sd(B)^2) )
## unreg_2nd_order_sd <- sqrt( (((1/305.8)^2)*(56^2)) + ( ((1220/(305.8^2))^2)*(59^2)) )
unreg_k_minus_alpha_sd <- sqrt( (((1/8.4)^2)*(9.8^2)) + ( ((73.2/(8.4^2))^2)*(0.8^2)) )

## unreg_2nd_order_label <- paste0(round(unreg_2nd_order, 1), "%+-%", round(unreg_2nd_order_sd, 1))
unreg_k_minus_alpha_label <- paste0(round(unreg_k_minus_alpha, 1), "%+-%", round(unreg_k_minus_alpha_sd, 1))


reg_k_minus_alpha <- 103.7/4.9
#                       sqrt( (((1/B)^2)*sd(A)^2) + (((A/B^2)^2)*sd(B)^2) )
## unreg_2nd_order_sd <- sqrt( (((1/305.8)^2)*(56^2)) + ( ((1220/(305.8^2))^2)*(59^2)) )
reg_k_minus_alpha_sd <- sqrt( (((1/4.9)^2)*(14^2)) + ( ((103.7/(4.9^2))^2)*(0.3^2)) )

## unreg_2nd_order_label <- paste0(round(unreg_2nd_order, 1), "%+-%", round(unreg_2nd_order_sd, 1))
reg_k_minus_alpha_label <- paste0(round(reg_k_minus_alpha, 1), "%+-%", round(reg_k_minus_alpha_sd, 1))

k_minus_alpha_p <- t.test2(m1 = unreg_k_minus_alpha, s1 = unreg_k_minus_alpha_sd, n1 = 3,
                        m2 = reg_k_minus_alpha, s2 = reg_k_minus_alpha_sd, n2 = 3)

(pretty_table_supp <-
  tribble(
   ~"Term",          ~"Unregulated",     ~"Regulated",               ~"P-Value",
   "K[+α]",         "8.4%+-%0.8", "4.9%+-%0.3", "<0.00001",
   "k[+α]~(s^-1)",  "73%+-%10", "104%+-%14", "0.01",
   "k[-α]~(s^-1)",  unreg_k_minus_alpha_label, reg_k_minus_alpha_label, as.character(round(k_minus_alpha_p["p-value"], 2))
   )|>
   ggtexttable(rows = NULL,
               theme = ttheme(
                              tbody.style = tbody_style(parse = TRUE,
                                                        fill = "white",
                                                        linecolor = "black"),
                              colnames.style = colnames_style(fill = "Black", color = "white"))) |>
  table_cell_bg(row = 2:4, column = 3, fill = alpha(colz[2], 0.3), color = "black") |>
    table_cell_bg(row = 2:4, column = 2, fill = alpha(colz[1], 0.3), color = "black")
   )


(gg_ratioA <-
   ggplot(ratio_df)+
   geom_point(aes(x = atp, y = ratioA, color = id), show.legend = FALSE)+
   geom_segment(data = ratio_plataeu,
                aes(x = 2000, xend = 6000,
                    y = avg, yend = avg,
                    color = id),
                linetype = "dashed",
                alpha = 0.4,
                show.legend = FALSE)+
   scale_color_manual(values = colz)+
   ggtitle("K<sub>&#945;</sub>")+
   ylab("k<sub>fast</sub>:k<sub>slow</sub>")+
   xlab("[ATP] &#40;&#xb5;M&#41;")+
   theme_cowplot(basesize)+
   theme(
     plot.title = element_markdown(),
     axis.title.y = element_markdown(),
     axis.title.x = element_markdown()
   )
)


# Finalize Plots

# table of values



scheme <- image_ggplot(image_read(here("code", "supplement", "scheme.pdf")))



fig_supp_1_top <- plot_grid(gg_kslow, gg_ratioA, labels = c("b", "c"))

png(filename = "img/supplement-figure-1.png", width = 6.75, height = 6, res = 300, units = "in")
plot_grid(scheme,
          fig_supp_1_top,
          pretty_table_supp,
          nrow = 3,
          rel_heights = c(0.2, 0.7, 0.5),
          labels = c("a", "", "d"))
dev.off()

svg(filename = "img/supplement-figure-1.svg", width = 6.75, height = 6)
plot_grid(scheme,
          fig_supp_1_top,
          pretty_table_supp,
          nrow = 3,
          rel_heights = c(0.2, 0.7, 0.5),
          labels = c("a", "", "d"))
dev.off()

cairo_ps(filename = "img/supplement-figure-1.eps", width = 6.75, height = 6)
plot_grid(scheme,
          fig_supp_1_top,
          pretty_table_supp,
          nrow = 3,
          rel_heights = c(0.2, 0.7, 0.5),
          labels = c("a", "", "d"))
dev.off()

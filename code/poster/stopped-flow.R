library(readxl)
library(ggplot2)
library(cowplot)
library(tibble)
library(tidyverse)
library(ggtext)
library(ggpubr)
dat <- lapply(4:7, \(x) read_xlsx("stopped-flow-data.xlsx", sheet = x))

adp_reg <- dat[[4]] |> mutate(id = "reg", y = `1-((Y1-C)/A)`) |> filter(Time <= 0.25)
adp_unreg <- dat[[1]] |> mutate(id = "1unreg", y = `1-((Y1-C)/A)`) |> filter(Time <= 0.25)
adp_dat <- rbind(adp_reg, adp_unreg)

colz <- c( "#007360","#a51417")


(gg_adp <- ggplot(adp_dat)+
  geom_line(aes(x = Time, y = y, color = id, size = id))+
  coord_cartesian(xlim = c(0, 0.25))+
  scale_size_manual(values = c(1, 0.5), guide = "none")+
  scale_color_manual(name = "", labels = c("Unregulated", "Regulated"), values = colz)+
  scale_x_continuous(expand = expansion(c(0, 0.02)))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), expand = expansion(c(0, 0), c(0, 0.11)))+
  xlab("Time (s)")+
  ylab("Normalized Fluorescence")+
  ggtitle("ADP Release")+
  theme_cowplot()+
  theme(
  legend.position = c(0.6, 0.4)
  )
)
## ggsave("poster/img/gg-adp.pdf")

atp_reg <- tribble(
~"conc", ~"kfast", ~"id",
11.26125, 32.67, "reg",
20.93325, 52.22,"reg",
54.34, 118.85,"reg",
104.2275, 231.92,"reg",
268.0275, 395.44,"reg",
555.1, 555.66,"reg",
1067.625, 759,"reg",
2092.35, 813.64,"reg",
3646.5, 954.26,"reg",
6015.75, 1031.47,"reg",
9.53225, 26.87,"reg",
18.0765, 44.88,"reg",
46.7025, 110.06,"reg",
91, 167.98,"reg",
234.91, 323.86,"reg",
476.775, 437.66,"reg",
946.075, 772.96,"reg",
2087.475, 860.91,"reg",
3357.25, 940.44,"reg",
5414.5, 1008.19, "reg")

atp_unreg <- tribble(
~"conc", ~"kfast", ~"id",
10.22775,38.13, "unreg",
22.568,85.56, "unreg",
59.7025,193.42, "unreg",
112.5475,354.91, "unreg",
278.07,663.38, "unreg",
610.675,788.53, "unreg",
1240.525,981.65, "unreg",
2226.9,1060.09, "unreg",
4042.5,1222.54, "unreg",
5885.75,1257.36, "unreg",
13.715,49.01087, "unreg",
21.4825,74.4418, "unreg",
90.025,287.6394, "unreg",
207.675,493.5051, "unreg",
367.25,640.2071, "unreg",
874.25,866.5379, "unreg",
1472.25,950.9917, "unreg",
2288,1112.106, "unreg",
2288,1076.4, "unreg",
4748.25,1128.845, "unreg"
)


atp_dat <- rbind(atp_reg, atp_unreg)

atp_dat$id <- factor(atp_dat$id, levels = c("unreg", "reg"))

atp_nest <- atp_dat |>
  group_by(id) |>
  nest() |>
  mutate(mod = map(data, ~nls(kfast ~ (k2*conc) / (k1+conc),
                              data = .x,
                              start = list(k1 = 1/305, k2 = 1000))),
         pred_df = map(mod, ~tibble(x = 1:6000,
                                    y = predict(.x, newdata = data.frame(conc = 1:6000)))))

atp_predict_df <- atp_nest |> select(id, pred_df) |> unnest(cols = c(pred_df))


(gg_all <-
ggplot()+
annotate("rect",
         xmin = -100,
         xmax = 305,
         ymin = -100,
         ymax = 320,
         alpha = 0.3,
         fill = "grey20")+
  geom_point(data = atp_dat,
             aes(x = conc,
                 y = kfast,
                 color = id))+
  geom_line(data = atp_predict_df,
            aes(x = x,
                y = y,
                color = id))+
  coord_cartesian(ylim = c(0, 1205))+
  ggtitle("ATP Binding")+
scale_color_manual(name = "", labels = c("Unregulated", "Regulated"), values = colz)+
  xlab("Time (s)")+
  ylab("k<sub>fast</sub> (s<sup>-1</sup>)")+
  theme_cowplot()+
  theme(
    legend.position = "none",
  axis.title.y = element_markdown())
)

(gg_short <-
  gg_all+
  coord_cartesian(xlim=c(0, 105), ylim=c(0, 300))+
  xlab("")+
  ylab("")+
  ggtitle("")+
  theme_cowplot(10)+
  theme(
  legend.position = "none"
  )
  )


gg_atp <-
  ggdraw(gg_all)+
  draw_plot(gg_short,
            x = 0.6,
            y = 0.15,
            width = 0.3,
            height = 0.5)


tab <- tribble(
  ~"Conditions", ~"ADP Release (s^-1)", ~"ATP Binding (M^-1*s^-1)",
  "Unregulated",            78.7,                   "3.99",
  "Regulated",              82.1,                   "2.18*"
  )

(ggtab <-
   ggtexttable(tab,
               theme = ttheme(colnames.style = colnames_style(parse = TRUE)))
)



plot_grid(gg_adp, gg_atp, nrow = 2)

ggsave("poster/img/gg-stopped-flow.pdf")

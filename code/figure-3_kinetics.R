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

set.seed(2022)
#### trap on times ####

colz <- unname(palette.colors())[c(1, 4)]

files <- list(
          Unregulated = here("data", "standard-trap", "unregulated_1uM-ATP_standard", "spasm-unregulated-standard.xlsx"),
          Regulated = here("data", "standard-trap" ,"regulated_1uM-ATP_standard", "spasm-regulated-standard-2.xlsx")
)

read_spasm_workbook <- function(spasm_file, id){
   map_df(excel_sheets(spasm_file),
           ~read_excel(spasm_file,
                       sheet = .x,
                       skip = 19)
        )|>
    mutate(id = id)

}

spasm_data <-
  imap_dfr(files,
           ~read_spasm_workbook(.x, .y)
  )

spasm_data$id <- factor(spasm_data$id, levels = c("Unregulated", "Regulated"))

spasm_data <- mutate(spasm_data,
                     "Attachment Times" = `duration (s)`)

# define a function to use mle to estimate detachment rate and boostrap CI
fit_ton <- function(data){
  #get attachement times
  ton <- data$`Attachment Times`
  fit_ton_pdf_mle <- function(ton, k){
  # pass the pars through to the negative log likelihood function
  # optim will optimize these
  # the variables ton will be inherited from parent function (no need to pass directly)
    nll_fun <- function(k){
      # PDF function from SPASM
      -sum(log(k*exp(-k*ton)))
    }

    fit <- optimize(nll_fun, lower = 0, upper = 100)

    return(fit)
  } #close fit_ton_pdf_mle

  # find k
  mod <- fit_ton_pdf_mle(ton = ton, k = 5)

  # extract fitted value from the model
  k <- mod$minimum[1]

  # function to generate fitted values to the exponential cumulative distribution
  fit_cdf <- function(x, k){ 1-exp(-k*x) }

  #calculate number of missing events
  n_missing <- fit_cdf(min(ton), k)*length(ton)

  #layout the range of x values for the real data bound by upper/lower bounds of data
  x_range <- seq(min(ton), max(ton), by = 1/1000)

  # generate the cdf values for the data
  cdf <- sapply(x_range, \(x) (sum(ton<=x)+n_missing) / (length(ton)+n_missing) )

  real_cdf <- data.frame(x = x_range,
                         y = cdf)



  # predict the fitted values from optimized points
  predict_x_range <- seq(0, max(ton), by = 1/1000)
  predict_y <- fit_cdf(k = k, x = predict_x_range)

  predict_df <- data.frame(x = predict_x_range,
                           y = predict_y)

#### BOOTSTRAP ####
  boostrap_ton_ci <- function(ton){

    boot_ton_fit <- function(ton){
      s <- sample(1:length(ton), replace = TRUE)
      new_ton <- ton[s]
      mod <- fit_ton_pdf_mle(ton = new_ton,
                              k = 5)
    } #close boot_ton_fit

    boot <- replicate(1000, boot_ton_fit(ton), simplify = FALSE)

    boot_df <- data.frame(k = sapply(boot, \(x) x$minimum[1]))

    ks <- sort(boot_df$k)

    k_lower <- ks[25]
    k_upper <- ks[975]

    return(list(boot_df = boot_df,
                k_95 = c(k_lower, k_upper)))

  } #close bootstrap_ton_ci

  ci <- boostrap_ton_ci(ton = ton)

  k_low_diff <- round(mod$minimum[[1]] - ci$k_95[[1]], 2)
  k_high_diff <- round(ci$k_95[[2]] - mod$minimum[[1]], 2)

  html_label <- paste0("k<sub>0</sub> = ",
                     round(mod$minimum[1], 1),
                     "s<sup>-1</sup> (+",
                     k_high_diff,
                     "/-",
                     k_low_diff,
                     ")")

  parse_label <- paste0(round(mod$minimum[1], 1),
                        "~(-",
                        k_low_diff,
                        "/+",
                        k_high_diff,
                        ")")

  list(
    data_df = real_cdf,
    mod = mod,
    predict_df = predict_df,
    boot_df = ci$boot_df,
    boot_ci = list(k1_low = k_low_diff,
                   k1_up = k_high_diff),
    html_label = html_label,
    parse_label = parse_label
  )
}

# fit the data
ton_boot <-
  spasm_data |>
  select(id, `Attachment Times`)|>
  group_by(id)|>
  nest() |>
  mutate(fit = map(data, fit_ton),
         predict = map(fit, `[[`, "predict_df"),
         parse_label  = map_chr(fit, `[[`, "parse_label"),
         boot_df = map(fit, `[[`, "boot_df"),
         cdf_data = map(fit, `[[`, "data_df"))

# get and save the bootstrap data for barrick-greenberg-ttest in matlab
## walk2(ton_boot$id, ton_boot$boot_df, ~readr::write_csv(x = .y,
##                                                 file = here("data",
##                                                             "standard-trap",
##                                                              paste0(.x, "_attachment-time-bootstrap.csv"))))

## ton_trap_pvalue <- kruskal.test(`Attachment Times` ~ id, data = spasm_data)$p.value

# unravel data from the nest back to long data for ggplot
# gets the real emperical CDF
ton_real_df <-
  ton_boot |>
  select(id, cdf_data) |>
  unnest(cols = c(cdf_data))

# unravel data from the nest back to long data for ggplot
# gets the fit prediction line
ton_predict_df <-
  ton_boot |>
  select(id, predict) |>
  unnest(cols = c(predict))

# make the plots
(ton_cdf <-
ggplot()+
  geom_step(data = ton_real_df,
            aes(x = x,
                y = y,
                color = id),
            alpha = 0.6,
            size = 1)+
  geom_line(data = ton_predict_df,
            aes(x = x,
                y = y,
                color = id),
            linetype = "dashed")+
  scale_color_manual(values = colz)+
  ggtitle("Attachment Times")+
  ylab("Cumulative Probability")+
  xlab("Time (s)")+
  theme_cowplot()+
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.title.y = element_text(size = 11)
)
)


## ADP Release
# get the stopped flow data from the excel sheet
adp_df <-
  map2_df(c("tidy_unreg_fluor", "tidy_reg_fluor"),
                  c("unreg", "reg"),
                  ~read_xlsx(here("data", "stopped-flow", "stopped-flow-data.xlsx"), sheet = .x) |>
                    mutate(id = .y, y = `1-((Y1-C)/A)`)) |>
  filter(Time <= 0.15)

adp_df$id <- factor(adp_df$id, levels = c("unreg", "reg"))

colz <- unname(palette.colors()[c(1, 4)])

(gg_adp <- ggplot(adp_df)+
  geom_line(aes(x = Time, y = y, color = id, size = id))+
  coord_cartesian(xlim = c(0, 0.15))+
  scale_size_manual(values = c(1, 0.5), guide = "none")+
  scale_color_manual(name = "", labels = c("Unregulated", "Regulated"), values = unname(palette.colors()[c(1, 4)]))+
  scale_x_continuous(expand = expansion(c(0, 0.02)))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), expand = expansion(c(0, 0), c(0, 0.11)))+
  xlab("Time (s)")+
  ylab("Normalized Fluorescence")+
  ggtitle("ADP Release")+
  theme_cowplot()+
  theme(
  plot.title = element_markdown(hjust = 0.5),
  legend.position = "none",
  axis.title.y = element_text(size = 10))
)

## kfast of ATP binding
## data from data/stopped-flow/stopped-flow-data.xlsx
kfast_reg <- tribble(
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

kfast_unreg <- tribble(
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


kfast_df <- rbind(kfast_reg, kfast_unreg)

kfast_df$id <- factor(kfast_df$id, levels = c("unreg", "reg"))

# fitting michalaes-menten for plotting purposes; values in paper used from data given
kfast_nest <- kfast_df |>
  group_by(id) |>
  nest() |>
  mutate(mod = map(data, ~nls(kfast ~ (k2*conc) / (k1+conc),
                              data = .x,
                              start = list(k1 = 1/305, k2 = 1000))),
         pred_df = map(mod, ~tibble(x = 1:6000,
                                    y = predict(.x, newdata = data.frame(conc = 1:6000)))))

kfast_predict_df <- kfast_nest |> select(id, pred_df) |> unnest(cols = c(pred_df))

unreg_atp_diss <- (1/305.8)*1220
reg_atp_diss <- (1/496.3)*1084

(gg_kfast <-
ggplot()+
annotate("rect",
         xmin = -100,
         xmax = 305,
         ymin = -100,
         ymax = 320,
         alpha = 0.3,
         fill = "grey20")+
  geom_point(data = kfast_df,
             aes(x = conc,
                 y = kfast,
                 color = id))+
  geom_line(data = kfast_predict_df,
            aes(x = x,
                y = y,
                color = id))+
  coord_cartesian(ylim = c(0, 1250))+
scale_color_manual(name = "", labels = c("Unregulated", "Regulated"), values = colz)+
  ggtitle("ATP Induced Dissociation")+
  xlab("[ATP] &#40;&#xb5;M&#41;")+
  ylab("k<sub>fast</sub> (s<sup>-1</sup>)")+
  theme_cowplot()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.position = "none",
  axis.title.y = element_markdown(),
  axis.title.x = element_markdown())
)

# make an inset plot
(gg_kfast_zoom <-
  gg_kfast+
  coord_cartesian(xlim=c(0, 105), ylim=c(0, 300))+
  xlab("")+
  ylab("")+
  ggtitle("")+
  theme_cowplot(8)+
  theme(
  legend.position = "none"
  )
  )

gg_kfast_zoom$layers[[4]] <- NULL

# draw the inset on the parent plot
(gg_kfast2 <- ggdraw(gg_kfast)+draw_plot(gg_kfast_zoom, 0.5, 0.2, 0.4, 0.5))

# adp release stats
adp_reg_mean <- 76
adp_reg_sd <- 5
adp_reg_n <- 3

adp_unreg <- c(79.96, 78.74, 77.43)
adp_unreg_mean <- mean(adp_unreg)
adp_unreg_sd <- sd(adp_unreg)
adp_unreg_n <- length(adp_unreg)

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
adp_p <- t.test2(m1 = adp_unreg_mean, s1 = adp_unreg_sd, n1 = adp_unreg_n,
                 m2 = adp_reg_mean, s2 = adp_reg_sd, n2 = adp_reg_n)

adp_unreg_label <- paste0(adp_unreg_mean, "%+-%", round(adp_unreg_sd, 2))
adp_reg_label <- paste0(adp_reg_mean, "%+-%", adp_reg_sd)

# calcualte second order binding rate and propogate error
# k_{+2} / (1/K_1)
unreg_2nd_order <- 1220/305.8
#                       sqrt( (((1/B)^2)*sd(A)^2) + (((A/B^2)^2)*sd(B)^2) )
unreg_2nd_order_sd <- sqrt( (((1/305.8)^2)*(56^2)) + ( ((1220/(305.8^2))^2)*(59^2)) )
unreg_2nd_order_label <- paste0(round(unreg_2nd_order, 2), "%+-%", round(unreg_2nd_order_sd, 1))

reg_2nd_order <- 1084/496.3
#                       sqrt( (((1/B)^2)*sd(A)^2) + (((A/B^2)^2)*sd(B)^2) )
reg_2nd_order_sd <- sqrt( (((1/496.3)^2)*(49^2)) + ( ((1084/(496.3^2))^2)*(82^2)) )
reg_2nd_order_label <- paste0(round(reg_2nd_order, 2), "%+-%", round(reg_2nd_order_sd, 1))

k_plus_ATP_p <- t.test2(m1 = unreg_2nd_order, s1 = unreg_2nd_order_sd, n1 = 3,
                        m2 = reg_2nd_order, s2 = reg_2nd_order_sd, n2 = 3)

# combine all data into a nice looking table
(pretty_table <-
  tribble(
   ~"Term",           ~"Unregulated",     ~"Regulated",               ~"P-Value",
   "t[on]~(s^-1)",     ton_boot$parse_label[[1]], ton_boot$parse_label[[2]], "0.002",
   "k[-ADP]~(s^-1)",   adp_unreg_label,             adp_reg_label,       as.character(round(adp_p[["p-value"]], 2)),
   "1/K[1]~(μM^-1)",           "305.8%+-%59", "496.3%+-%82",     "0.008",
   "k[+2]~(s^-1)",     "1220%+-%56", "1084%+-%49",         "0.01",
   "k[K[1]%.%k[+2]]~(μM^-1*s^-1)", unreg_2nd_order_label,  reg_2nd_order_label,   as.character(round(k_plus_ATP_p[["p-value"]], 2))
   )|>
   ggtexttable(rows = c("a", "b", "c", "c", "c"),
               theme = ttheme(base_size = 14,
                              tbody.style = tbody_style(size = 14,
                                                        parse = TRUE,
                                                        fill = "white",
                                                        linecolor = "black"),
                              colnames.style = colnames_style(size = 14,
                                                              face = "bold",
                                                              fill = "Black",
                                                              color = "white"))) |>
  table_cell_bg(row = 2:6, column = 4, fill = alpha(colz[2], 0.3), color = "black") |>
    table_cell_bg(row = 2:6, column = 3, fill = alpha(colz[1], 0.3), color = "black")
   )



# plot everything together
fig3_top <- plot_grid(ton_cdf, gg_adp, gg_kfast2, ncol = 3, labels = letters)

plot_grid(fig3_top, pretty_table, ncol = 1, labels = c("", "d"), rel_heights = c(0.5, 0.5))


## ggsave(here("img", "fig-3_kinetics.png"), dpi=500, bg = "white")

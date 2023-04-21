# this is mostly copy paste from figure-2 and figure-3 but includes the low calcium data
library(data.table)
library(dplyr)
library(purrr)
library(readxl)
library(here)
library(cowplot)
library(ggtext)
library(ggpubr)
library(tidyr)

textsize <- 9
basesize <- 11
# get file locations for the SPASM output files
# these files contain all the data for the analyzed trap events
files <- list(
          Unregulated = here("data", "standard-trap", "unregulated_1uM-ATP_standard", "spasm-unregulated-standard.xlsx"),
          "pCa 4" = here("data", "standard-trap" ,"regulated_1uM-ATP_standard", "spasm-regulated-standard.xlsx"),
          "pCa 6.25" = here("data", "standard-trap" ,"pCa-6.25", "spasm-6.25-v7.xlsx")
)

# defines a helper function to help read in SPASM trap data
read_spasm_workbook <- function(spasm_file, id){
   map_df(excel_sheets(spasm_file),
           ~read_excel(spasm_file,
                       sheet = .x,
                       skip = 19)
        )|>
    mutate(id = id)

}

# read in data
spasm_data <-
  imap_dfr(files,
           ~read_spasm_workbook(.x, .y)
  )

spasm_data$id <- factor(spasm_data$id, levels = c("Unregulated", "pCa 4", "pCa 6.25"))

# make new columns based on avg of 2 beads position / rename a few
spasm_data <- mutate(spasm_data,
                     "Displacements" = (`A total step (nm)` + `B total step (nm)`)/2,
                     "Attachment Times" = `duration (s)`,
                     "Substep 1" = (`A step 1 (nm)` + `B step 1 (nm)`) / 2,
                     "Substep 2" = (`A step 2 (nm)` + `B step 2 (nm)`) / 2)


spasm_sum <-
  spasm_data |>
  group_by(id) |>
  summarize(
    total_step_avg = round(mean(Displacements), 1),
    total_step_sd= round(sd(Displacements), 1),
    step_1_avg = round(mean(`Substep 1`), 1),
    step_1_sd= round(sd(`Substep 1`), 1),
    step_2_avg = round(mean(`Substep 2`), 1),
    step_2_sd= round(sd(`Substep 2`), 1),
    n = n())


#### Raw traces ###
colz <- c(unname(palette.colors())[c(1, 4)], "#c4b100ff")


#### Ensemble Averages ####
# define a helper function to read id ensemble average data from SPASM
# and refit the forwards and backwards ensemble to smaller subset of data
spasm_plot_ensemble_average <- function(spasm_file, color, x_shift, title){

  ensemble_data <- read_xlsx(spasm_file,
                             range = cell_cols(1:9))

  numbers <- read_xlsx(spasm_file,
                       range = "K1:L5",
                       col_names =  c("key", "value"))

  new_names <- c("seconds", "nanometers", "id", "direction")

  select_cols <- list(1:2,
                      3:4,
                      6:7,
                      8:9)

  id <- rep(c("data", "fit"), 2)
  direction <- c("f", "f", "b", "b")

  ea_data <-
 pmap(list(select_cols, id, direction),
      ~select(ensemble_data, all_of(..1)) |> mutate(id = ..2,
                                                           direction = ..3))

  ea_data <- map_df(ea_data, ~set_names(.x, new_names))

#### fit to subset of data
  df_subset <-
    ea_data |>
    filter(direction == "f" & id == "data" & between(seconds, 0, 0.2 ))

  ea_forward_fit <-
    nls(nanometers ~ d1 + (d2*(1 - exp(-k1 * seconds))),
        data = df_subset,
        start = list(d1 = 4, d2 = 1, k1 = 25))

  ## print(summary(ea_forward_fit))

  ## print(coef(ea_forward_fit)["k1"])

  f_pred <- df_subset |> mutate(y_fit = predict(ea_forward_fit))

  df_subset_back <-
    ea_data |>
    filter(direction == "b" & id == "data" & between(seconds, -0.5, 0 ))

  ea_backwards_fit <-
    nls(nanometers ~ d1+(d2*exp(seconds*k2)),
        data = df_subset_back,
        start = list(d1 = 4, d2 = 1, k2 = 5))


  b_pred <- df_subset_back |> mutate(y_fit = predict(ea_backwards_fit),
                                     seconds = seconds + x_shift)

  ## print(summary(ea_backwards_fit))

  ea_f <- ea_data |> filter(direction == "f" & id == "data" & between(seconds, -0.1, 0.6))
  ea_b <- ea_data |> filter(direction == "b" & id == "data" & between(seconds, -0.6, 0.1)) |> mutate(seconds = seconds + x_shift)

  d2_y_position <- (numbers$value[[3]] + numbers$value[[5]]) / 2
  d1_y_position <- (numbers$value[[3]] + 0) / 2

  x_position <- (min(ea_b$seconds)+ max(ea_f$seconds)) / 2
  d1_triangle_y <- numbers$value[[3]]*0.9

  k1_y <- max(ea_f$nanometers)+0.15
  k2_y <- max(ea_b$nanometers)+0.15

  mod_f <- ea_forward_fit
  mod_b <- ea_backwards_fit

  gg <-
    ggplot()+
    geom_line(data = ea_f,
              aes(x = seconds,
                  y = nanometers),
              color = color,
              ## alpha = 0.8,
              linewidth = 0.3,
              show.legend = FALSE)+
    geom_line(data = f_pred,
              aes(x = seconds,
                  y = y_fit))+
    geom_line(data = ea_b,
              aes(x = seconds,
                  y = nanometers),
              color = color,
              ## alpha = 0.8,
              linewidth = 0.3,
              show.legend = FALSE)+
    geom_line(data = b_pred,
              aes(x = seconds,
                  y = y_fit))+
    annotate("point",
             x = x_position,
             y = d2_y_position,
             shape = 2)+
    annotate("point",
             x = x_position,
             y = d1_triangle_y,
             shape = 17)+
    annotate("richtext",
             x = x_position-0.04,
             y = d1_y_position,
             label = paste0("δ<sub>1</sub> = ",
                            round(numbers$value[[3]], 1),
                            " nm"),
             fill = NA,
             label.color = NA,
             hjust = 1,
             size = textsize/.pt)+
    annotate("richtext",
             x = x_position+0.02,
             y = d2_y_position+0.2,
             label = paste0("δ<sub>2</sub> = ", round(numbers$value[[4]], 1),
                            " nm"),
             fill = NA,
             label.color = NA,
             hjust = 0,
             size = textsize/.pt)+
    annotate("richtext",
             x = 0,
             y = k1_y,
             label = paste0("k<sub>f</sub> = ",
                            round(coef(mod_f)["k1"], 0),
                            " s<sup>-1</sup>"),
             fill = NA,
             label.color = NA,
             hjust = 0.5,
             size = textsize/.pt)+
    annotate("richtext",
             x = x_shift,
             y = k2_y,
             label = paste0("k<sub>r</sub> = ",
                            round(coef(mod_b)["k2"], 0),
                            " s<sup>-1</sup>"),
             fill = NA,
             label.color = NA,
             hjust = 0.5,
             size = textsize/.pt)+
    annotate("segment",
             x = x_position,
             y = 0 ,
             xend = x_position,
             yend = d1_triangle_y,
             linetype = "dotted")+
    ##scale bars
    draw_line(x = c(0.1, 0.35),
              y = rep(min(ea_f$nanometers), 2),
              linewidth = 0.7)+
    draw_line(x = rep(-0.12, 2),
              y = c(0, 2),
              linewidth = 0.7)+
    annotate("text",
             x = 0.225,
             y = min(ea_f$nanometers),
             label = "0.25 s",
             vjust = 1.25,
             hjust = 0.5,
             color = "black",
             size = textsize/.pt)+
    annotate("text",
             x = -0.12,
             y = 1,
             label = "2 nm",
             angle = 90,
             vjust = -0.5,
             hjust = 0.5,
             color = "black",
             size = textsize/.pt)+
    ggtitle(title)+
    ylab("Position (nm)")+
    xlab("Time (s)")+
    scale_color_manual(values = c(color, "black"))+
    coord_cartesian(ylim = c(min(ea_f$nanometers)-0.25, NA),
                    xlim = c(-0.122, NA))+
    theme_cowplot(basesize)+
    theme(
      plot.title = element_text(hjust = 0.5, color = "transparent"),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank()
    )


  return(list(gg = gg, mod_f  = mod_f, mod_b = mod_b))
}


g1 <-
  spasm_plot_ensemble_average(here("data", "standard-trap", "ensemble-averages", "pCa-625-combinedEnsembleAxesData-v7.xlsx"),
                                 color = colz[3],
                                 x_shift = 1.2,
                                 title = "pCa 6.25")


## g1


#### Step ####

# plot the displacements and substeps
# helper function to expedite the base plot
spasm_plot_ecdf <- function(spasm_data, var, colz, x_lab = "x_label", legend = FALSE){

  ggplot(spasm_data)+
    stat_ecdf(aes(base::get(var),
                  color = id,
                  alpha = id),
              linewidth = 1,
              pad = FALSE,
              show.legend = legend)+
    scale_color_manual(name = "",values = colz)+
    scale_alpha_manual(values = c(0.5, 0.5, 1), guide = FALSE)+
    ylab("Cumulative Distribution")+
    xlab(x_lab)+
    ggtitle(var)+
    theme_cowplot(basesize)+
    theme(
      plot.title = element_text(hjust = 0.5, size = basesize)
    )
}

spasm_low_ca <- filter(spasm_data, id == "pCa 6.25")
conditions_to_compare <- c("Unregulated", "pCa 4") #regulated is pCa 4
var_to_test <- c("Displacements", "Substep 1", "Substep 2")
tt_test_results <- list(Unregulated = list(),
                        Regulated = list())
for(cond in seq_along(conditions_to_compare)){
  data_to_compare <- dplyr::filter(spasm_data, id == conditions_to_compare[[cond]])
 for(var in seq_along(var_to_test)){
   var_name <- var_to_test[[var]]
   tt_test_results[[cond]][[var]] <- t.test(spasm_low_ca[[var_name]], data_to_compare[[var_name]])
 }
  names(tt_test_results[[cond]]) <- var_to_test
}

(total_step_ecdf <- spasm_plot_ecdf(spasm_data, var = "Displacements", colz, x_lab = "nanometers")+
  ggtitle("Total Displacements")+
  ylab("Cumulative Probability")+
  xlab("Displacement (nm)")+
 annotate("richtext",
          x = -Inf,
          y = Inf,
          label = paste0("<span style = 'color:", colz[[1]],
                             "'> x̄ = ", spasm_sum$total_step_avg[[1]],
                             " ± ", spasm_sum$total_step_sd[[1]], " nm",
                             "</span> <br>",
                         "<span style = 'color:", colz[[2]],
                             "'> x̄ = ", spasm_sum$total_step_avg[[2]],
                             " ± ", spasm_sum$total_step_sd[[2]], " nm",
                             "</span> <br>",
                         "<span style = 'color:", colz[[3]],
                             "'> x̄ = ", spasm_sum$total_step_avg[[3]],
                             " ± ", spasm_sum$total_step_sd[[3]], " nm",
                             "</span> <br>"),
          ## color = colz[[3]],
          ## parse = TRUE,
          fill = "transparent",
          color = "transparent",
          hjust = 0,
          vjust = 1,
          size = textsize/.pt)+
 annotate("richtext",
          x = Inf,
          y = -Inf,
          label = paste0("<span style = 'color:", colz[[3]],"'>pCa 6.25</span>",
                         "<span style = 'color: black'> vs </span>",
                         "<span style = 'color:", colz[[2]],"'>pCa 4</span>",
                         "<span style = 'color: black'> (P = ",
                         round(tt_test_results$Regulated$Displacements$p.value, 2), ")</span>",
                         "<br>",
                         "<span style = 'color:", colz[[3]],"'>pCa 6.25</span>",
                         "<span style = 'color: black'> vs </span>",
                         "<span style = 'color:", "#808080","'>Unreg</span>",
                         "<span style = 'color: black'> (P = ",
                         round(tt_test_results$Unregulated$Displacements$p.value, 2), ")</span>"
                         ),
          ## color = colz[[3]],
          ## parse = TRUE,
          fill = "transparent",
          color = "transparent",
          hjust = 1,
          vjust = 0,
          size = 6/.pt)+
   coord_cartesian(xlim = c(-32, NA))+
   theme(axis.title.x = element_text(color = "white")
         )
  )


(one_step_ecdf <- spasm_plot_ecdf(spasm_data, var = "Substep 1", colz, x_lab = "nanometers")+
  ylab("Cumulative Probability")+
  xlab("Displacement (nm)")+
 annotate("richtext",
          x = -Inf,
          y = Inf,
          label = paste0("<span style = 'color:", colz[[1]],
                             "'> x̄ = ", spasm_sum$step_1_avg[[1]],
                             " ± ", spasm_sum$step_1_sd[[1]], " nm",
                             "</span> <br>",
                         "<span style = 'color:", colz[[2]],
                             "'> x̄ = ", spasm_sum$step_1_avg[[2]],
                             " ± ", spasm_sum$step_1_sd[[2]], " nm",
                             "</span> <br>",
                         "<span style = 'color:", colz[[3]],
                             "'> x̄ = ", spasm_sum$step_1_avg[[3]],
                             " ± ", spasm_sum$step_1_sd[[3]], " nm",
                             "</span> <br>"),
          fill = "transparent",
          color = "transparent",
          hjust = 0,
          vjust = 1,
          size = textsize/.pt)+
 annotate("richtext",
          x = Inf,
          y = -Inf,
          label = paste0("<span style = 'color:", colz[[3]],"'>pCa 6.25</span>",
                         "<span style = 'color: black'> vs </span>",
                         "<span style = 'color:", colz[[2]],"'>pCa 4</span>",
                         "<span style = 'color: black'> (P = ",
                         round(tt_test_results$Regulated$`Substep 1`$p.value, 2), ")</span>",
                         "<br>",
                         "<span style = 'color:", colz[[3]],"'>pCa 6.25</span>",
                         "<span style = 'color: black'> vs </span>",
                         "<span style = 'color:", "#808080","'>Unreg</span>",
                         "<span style = 'color: black'> (P = ",
                         round(tt_test_results$Unregulated$`Substep 1`$p.value, 2), ")</span>"
                         ),
          fill = "transparent",
          color = "transparent",
          hjust = 1,
          vjust = 0,
          size = 6/.pt)+
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank())
  )



(two_step_ecdf <- spasm_plot_ecdf(spasm_data, var = "Substep 2", colz, x_lab = "nanometers", legend = TRUE)+
  ## ggtitle("Total Displacements")+
  ylab("Cumulative Probability")+
  xlab("Displacement (nm)")+
 annotate("richtext",
          x = -Inf,
          y = Inf,
          label = paste0("<span style = 'color:", colz[[1]],
                             "'> x̄ = ", spasm_sum$step_2_avg[[1]],
                             " ± ", spasm_sum$step_2_sd[[1]], " nm",
                             "</span> <br>",
                         "<span style = 'color:", colz[[2]],
                             "'> x̄ = ", spasm_sum$step_2_avg[[2]],
                             " ± ", spasm_sum$step_2_sd[[2]], " nm",
                             "</span> <br>",
                         "<span style = 'color:", colz[[3]],
                             "'> x̄ = ", spasm_sum$step_2_avg[[3]],
                             " ± ", spasm_sum$step_2_sd[[3]], " nm",
                             "</span> <br>"),
          fill = "transparent",
          color = "transparent",
          hjust = 0,
          vjust = 1,
          size = textsize/.pt)+
 annotate("richtext",
          x = Inf,
          y = -Inf,
          label = paste0("<span style = 'color:", colz[[3]],"'>pCa 6.25</span>",
                         "<span style = 'color: black'> vs </span>",
                         "<span style = 'color:", colz[[2]],"'>pCa 4</span>",
                         "<span style = 'color: black'> (P = ",
                         round(tt_test_results$Regulated$`Substep 2`$p.value, 2), ")</span>",
                         "<br>",
                         "<span style = 'color:", colz[[3]],"'>pCa 6.25</span>",
                         "<span style = 'color: black'> vs </span>",
                         "<span style = 'color:", "#808080","'>Unreg</span>",
                         "<span style = 'color: black'> (P = ",
                         round(tt_test_results$Unregulated$`Substep 2`$p.value, 2), ")</span>"
                         ),
          ## color = colz[[3]],
          ## parse = TRUE,
          fill = "transparent",
          color = "transparent",
          hjust = 1,
          vjust = 0,
          size = 6/.pt)+
   coord_cartesian(xlim = c(-25, NA))+
  theme(
    legend.position = "none",
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(color = "white"))
  )

# combine step plots, figure 2 bottom row
(gg_steps <- plot_grid(total_step_ecdf,
                       one_step_ecdf,
                       two_step_ecdf,
                       nrow = 1,
                       rel_widths = c(1, 0.9, 0.9),
                       labels = c("a", "b", "c")))



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

  html_label <- paste0("<i>k<sub>det</sub></i> = ",
                     round(mod$minimum[1], 1),
                     " (-",
                     round(k_low_diff, 1),
                     "/+",
                     round(k_high_diff, 1),
                     ") s<sup>-1</sup>")

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
         html_label  = map_chr(fit, `[[`, "html_label"),
         boot_df = map(fit, `[[`, "boot_df"),
         cdf_data = map(fit, `[[`, "data_df"))

# get and save the bootstrap data for barrick-greenberg-ttest in matlab

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
                color = id,
                alpha = id),
            linewidth = 1)+
  geom_line(data = ton_predict_df,
            aes(x = x,
                y = y,
                color = id),
            linetype = "dashed")+
 annotate("richtext",
          x = 1.35,
          y = 0.5,
          label = paste0("<span style = 'color:", colz[[1]], "'>", ton_boot$html_label[[1]], "</span> <br>",
                         "<span style = 'color:", colz[[2]], "'>", ton_boot$html_label[[2]], "</span> <br>",
                         "<span style = 'color:", colz[[3]], "'>", ton_boot$html_label[[3]], "</span> <br>"),
          ## color = colz[[3]],
          ## parse = TRUE,
          fill = "white",
          color = "white",
          hjust = 0.5,
          vjust = 0.5,
          size = textsize/.pt)+
  scale_color_manual(values = colz)+
  scale_alpha_manual(values = c(0.3, 0.3, 0.75))+
  ggtitle("Attachment Durations")+
  ylab("Cumulative Probability")+
  xlab("Time (s)")+
  theme_cowplot(basesize)+
  theme(
    plot.title = element_text(hjust = 0.5, size = basesize),
    legend.position = "none",
    ## axis.title.y = element_text(size = 11)
)
)

bottom <- plot_grid(ton_cdf, g1$gg, labels = c("d", "e"), nrow = 1 )

png(filename="img/figure-s2_pCa-625-trap.png", width = 7, height = 5, units = "in", res = 500)
plot_grid(gg_steps, bottom, nrow = 2)
dev.off()



## get and save the bootstrap data for barrick-greenberg-ttest in matlab
## walk2(ton_boot$id, ton_boot$boot_df, ~readr::write_csv(x = .y,
##                                                 file = here("data",
##                                                             "standard-trap",
##                                                              paste0(.x, "_attachment-time-bootstrap-2023-04-17.csv"))))

here::i_am("code/figure-4_isometric-force-clamp.R")
library(here)
library(readr)
library(purrr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(cowplot)
library(ggtext)
library(data.table)

set.seed(2022)
textsize <- 9
basesize <- 11

## make function to read in analyzed feedback data from SPASM
read_feedback_data <- function(spasm_xlsx_path){
  # determine which bead was motor vs transducer
  # this info is in the trap .txt file, and not in spasm
  # use the spasm.xlsx to look for the trap txt in the same folder
  # and read the header info
  obs_dir <- sub("spasm.xlsx", "", spasm_xlsx_path)

  trap_file <- list.files(obs_dir, pattern = ".txt", full.names = TRUE)

  trap_header <-
    read_tsv(trap_file,
             n_max = 68,
             col_names = c("key", "value"),
             show_col_types = FALSE)

  # line 30 of trap header contains info on which bead was motor
  motor_bead <- trap_header |> slice(30) |> pull(value) |> as.integer()

  # read the the spasm output file skipping the header
  spasm_data <-
    read_xlsx(spasm_xlsx_path,
              sheet = 1,
              skip = 19)

  # extract information from spasm.xlsx based on correct motor bead
  if(motor_bead == 1){
    new_data <-
      tibble(
        time_on_s = spasm_data$`duration (s)`,
        time_on_ms = time_on_s/1000,
        force_before = spasm_data$`A avg force before (pN)`,
        force_during = spasm_data$`A avg force during (pN)`,
        force_avg = force_during - force_before
      )
  } else if(motor_bead == 2){
    new_data <-
      tibble(
        time_on_s = spasm_data$`duration (s)`,
        time_on_ms = time_on_s/1000,
        force_before = spasm_data$`B avg force before (pN)`,
        force_during = spasm_data$`B avg force during (pN)`,
        force_avg = force_during - force_before
      )
  } else {
    stop("motor_bead did not have a matching value")
  }

  return(new_data)
}


## get all feedback SPASM excel files names and read in
feedback_files <- list.files(here("data", "isometric-force-clamp"),
                                                pattern = "spasm.xlsx",
                                                recursive = TRUE,
                                                full.names = FALSE)


feedback_data <- map_df(feedback_files,
                        ~read_feedback_data(here("data", "isometric-force-clamp", .)) |>
                          mutate(path = .))|>
  separate(path, c("conditions",
                   "date",
                   "obs",
                   "file"),
           sep = "/")|>
  unite("date_obs",
        c(date, obs),
        sep = "_",
        remove = FALSE) |>
  filter(force_avg >=  -1)



feedback_data$conditions <-  factor(feedback_data$conditions,
                                    levels = c("unregulated_1000uM-ATP_feedback",
                                               "regulated_1000uM-ATP_feedback"),
                                    labels = c("Unregulated", "Regulated"))

n_df <- feedback_data |> group_by(conditions) |> summarize(n = n())



reg_data_for_memlet <- feedback_data  |> filter(conditions == "Regulated") |> dplyr::select(time_on_s, force_avg)

## write_csv(reg_data_for_memlet, file = here("data", "isometric-force-clamp", "reg-data-for-memlet.csv"))
# uncomment for a diagnostic plot per molecule
## ggplot()+
##     geom_point(data = filter(feedback_data, conditions == "Regulated"),
##   ## geom_point(data = feedback_data,
##                aes(force_avg,
##                    time_on_s,
##                    color = date_obs),
##                alpha = 0.5,
##                shape = 16,
##                size = 1.5,
##                show.legend = T) + facet_wrap(~date_obs)+ theme(legend.position = "none")+
##   scale_y_log10()

# define a function to fit bell equation with mle and bootstrap CI
fit_bell <- function(data, tmin = NULL){

  if(!is.null(tmin)){
    data <- data[which(data$time_on_s >= tmin),]
  }

  fit_bell_pdf_mle <- function(F, t, pars, tmin){
  # pass the pars through to the negative log likelihood function
  # optim will optimize these
  # the variables F, t will be inherited from parent function (no need to pass directly)
    nll_fun <- function(pars, tmin){
      k0 <- pars[1]
      d <- pars[2]
      if(is.null(tmin)){
    # PDF function from MEMLET https://github.com/michaelswoody/MEMLET/blob/master/Matlab%20Code/MEMLET/PDFs/bellsEqn.m
      -sum(log(((k0*exp(-(F*d)/4.1))*exp(-(k0*exp(-(F*d)/4.1))*t))))
    } else {
      -sum(
         log(
           ((k0*exp(-(F*d)/4.1))*exp(-(k0*exp(-(F*d)/4.1))*t)) / (exp(-(k0*exp(-(F*d)/4.1))*tmin))
         )
       )
     }
    }

    fit <- optim(pars, nll_fun, tmin = tmin)
    return(fit)
  } #close fit_bell_pdf_mle

  mod <- fit_bell_pdf_mle(F = data$force_avg, t = data$time_on_s, pars = c(50, 1), tmin = tmin)

  predict_bell <- function(mod, F){
    k0 <- 1/mod$par[1]
    d <- mod$par[2]
    dummy_F <- seq(min(F), max(F), by = 1/100)
    t <- k0*exp((dummy_F*d)/4.1)
    predict_df <- data.frame(F = dummy_F,
                             t = t)
  } #close predict_bell

  predict_df <- predict_bell(mod, data$force_avg)

#### BOOTSTRAP ####

  boostrap_bell_ci <- function(F, t, tmin){

    data <- data.frame(F, t)

    boot_bell_fit <- function(data, tmin){
      s <- sample(1:nrow(data), replace = TRUE)
      df <- data[s,]
      mod <- fit_bell_pdf_mle(F = df$F,
                              t = df$t,
                              pars = c(50, 1),
                              tmin = tmin)
    } #close boot_bell_fit

    boot <- replicate(1000, boot_bell_fit(data, tmin), simplify = FALSE)

    boot_df <- data.frame(k0 = sapply(boot, \(x) x$par[1]),
                          d = sapply(boot, \(x) x$par[2]))

    k0s <- sort(boot_df$k0)
    ds <- sort(boot_df$d)

    k0_lower <- k0s[25]
    k0_upper <- k0s[975]

    d_lower <- ds[25]
    d_upper <- ds[975]

    return(list(boot_df = boot_df,
                k0_95 = c(k0_lower, k0_upper),
                d_95 = c(d_lower, d_upper)))

  } #close bootstrap_bell_ci

  ci <- boostrap_bell_ci(data$force_avg, data$time_on_s, tmin)

  ci_low <- list(par = c(ci$k0_95[2], ci$d_95[1]))
  ci_high <- list(par = c(ci$k0_95[1], ci$d_95[2]))

  ## fit_ci_low <- predict_bell(ci_low, F = data$force_avg)
  ## fit_ci_high <- predict_bell(ci_high, F = data$force_avg)

  ## ci_df <- data.frame(F = fit_ci_low$F,
  ##                     tmin = fit_ci_low$t,
  ##                     tmax = fit_ci_high$t)


  k0_low_diff <- round(mod$par[[1]] - ci$k0_95[[1]], 0)
  k0_high_diff <- round(ci$k0_95[[2]] - mod$par[[1]], 0)
  d_low_diff <- format(round(mod$par[[2]] - ci$d_95[[1]], 2), nsmall = 2)
  d_high_diff <- format(round(ci$d_95[[2]] - mod$par[[2]], 2), nsmall = 2)

  label <- paste0("k<sub>0</sub> = ",
                     round(mod$par[1], 0),
                     " (+",
                     k0_high_diff,
                     "/-",
                     k0_low_diff,
                     ") s<sup>-1</sup>",
                     "<br>",
                     "d = ",
                     round(mod$par[2], 2),
                     " (+",
                     d_high_diff,
                     "/-",
                     d_low_diff,
                     ") nm")

  list(
    mod = mod,
    predict_df = predict_df,
    boot_df = ci$boot_df,
    boot_ci = list(k1_low = k0_low_diff,
                   k1_up = k0_high_diff,
                   d_low = d_low_diff,
                   d_up = d_high_diff),
    label = label
  )
}


# unravel to nested to data back to long format for ggplot
bell_df <-
  feedback_data |>
  dplyr::select(conditions, force_avg, time_on_s) |>
  group_by(conditions) |>
  nest() |>
  mutate(fit = map(data, ~fit_bell(.x)),
         predict = map(fit, `[[`, "predict_df"),
         label  = map_chr(fit, `[[`, "label")
         )

## for saving bootstrap parameteters for hypothesis testing in matlab
boot_bell_df <-
  bell_df |>
  mutate(boot = map(fit, `[[`, "boot_df"),
         label  = map_chr(fit, `[[`, "label")
         ) |>
  select(conditions, boot)

## walk2(boot_bell_df$conditions, boot_bell_df$boot, ~write_csv(x = .y,
##                                                              file = here("data",
##                                                                       "isometric-force-clamp",
##                                                                       paste0(tolower(.x), "_bootstrap-deadtime-10ms.csv"))))

## unravel prediction line from the nest
predict_df <-
  bell_df |>
  select(conditions, predict) |>
  unnest(cols = c(predict))

colorz <- unname(palette.colors()[c(1, 4)])

fig4_bottom <-
  ggplot()+
    geom_point(data = feedback_data,
               aes(x = force_avg,
                   y = time_on_s,
                   color = conditions),
               alpha = 0.4,
               shape = 16,
               size = 1)+
    geom_richtext(data = bell_df,
              aes(x = -Inf,
                  y = Inf,
                  label = label,
                  color = conditions),
              label.color = "white",
              fill = NA,
               hjust = 0,
              vjust = 1,
              size = textsize/.pt)+
    facet_wrap(~conditions)+
    scale_color_manual(values = colorz)+
    ylab("Time (s)")+
    xlab("Force (pN)")+
    scale_y_log10()+
    theme_cowplot(basesize)+
    theme(
    legend.position = "none",
    strip.background = element_rect("white"),
    strip.text = element_blank()
)


#### raw traces ####

read_greenberg <- function(file){

  raw_data <- fread(file, skip = 68)
  header_data <- fread(file, nrows = 67, header = FALSE)

  header_line_numbers <- c(15, 18, 20, 22, 24)
  header_data <- header_data[header_line_numbers,]
  options <- as.list(as.numeric(header_data$V2))
  names(options) <- c("hz", "pn_nm1", "pn_nm2", "nm_v1", "nm_v2")

  list(data = raw_data,
       options = as.data.frame(options))

}

start_unreg <- 8.75
stop_unreg <- 9.98

folder_path <- c(here("data/isometric-force-clamp/unregulated_1000uM-ATP_feedback/2022-08-26/obs-01"),
                 here("data/isometric-force-clamp/regulated_1000uM-ATP_feedback/2022-12-07/obs-03"))

plot_isometric_force_clamp <- function(folder_path, start, stop, title, color, offset = 0){

  file <- list.files(folder_path,
                     full.names = TRUE)
  trap <- read_greenberg(file[[1]])
  spasm <- read_excel(file[[2]], skip = 19)

  hz <- 20000
  start <- start*hz
  stop <- stop*hz

  spasm <- setDT(spasm)[`index before start (# pts)` >= start & `index before start (# pts)` <= stop]
  spasm[, `:=`(event_start = `index before start (# pts)` - start,
               event_stop = `index before end (# pts)` - start)]

  trap_data <- trap$data[start:stop,]
  trap_data$time_col <- 1:nrow(trap_data)/hz
  trap_data$Trap1X <- trap_data$Trap1X * trap$options$nm_v1 * trap$options$pn_nm1
  trap_data$Trap2X <- (trap_data$Trap2X * trap$options$nm_v2 * trap$options$pn_nm2 ) + offset

  gg <-
     ggplot(trap_data, aes(x = time_col))+
     geom_line(aes(y = Trap1X), color = color, linewidth = 0.2)+
     geom_line(aes(y = Trap2X), color = alpha(color, 0.6), size = 0.2)+
     xlab("Time (s)")+
     ylab("Force (pN)")+
     ggtitle(title)+
     theme_minimal(basesize)+
     theme(
       plot.title = element_text(face = "bold", hjust = 0.5))

  return(gg)

 }

unreg_n <- n_df$n[n_df$conditions == "Unregulated"]

(unreg_raw <-
  plot_isometric_force_clamp(folder_path[[1]], start_unreg, stop_unreg, "Unregulated", colorz[[1]] )+
  coord_cartesian(ylim = c(-3.5, NA))+
  annotate("segment", x = 0, xend = 0.25, y = -3, yend=-3)+
  annotate("segment", x = -0.05, xend = -0.05, y = 4, yend = 8)+
  annotate("text", x = 0.125, y = -3, label = "0.25 s", vjust = 1.3, size = textsize/.pt)+
  annotate("text", x = -0.05, y = 6, label = "4 pN", vjust = -0.4 , angle = 90, size =  textsize/.pt)+
  annotate("text", x = 1.25, y = 2.8, label = "M", hjust = 0 , vjust = 0 , angle = 0, size = textsize/.pt)+
  annotate("text", x = 1.25, y = -0.8, label = "T", hjust = 0 , vjust = 0 , angle = 0, size = textsize/.pt)+
  coord_cartesian(xlim = c(-0.13, 1.35), ylim = c(-3.7, NA))+
  theme_void(basesize)+
     theme(
       plot.title = element_text(face = "bold", hjust = 0.5))
)

unreg_raw <- unreg_raw+ggtitle(paste0("Unregulated (n = ", unreg_n, ")"))

reg_n <- n_df$n[n_df$conditions == "Regulated"]

(reg_raw <-
plot_isometric_force_clamp(folder_path[[2]], start = 97, stop = 98.2, "Regulated", colorz[[2]], offset = 1.5)+
  coord_cartesian(ylim = c(-3.5, NA))+
  annotate("segment", x = 0, xend = 0.25, y = -3, yend=-3)+
  annotate("segment", x = -0.05, xend = -0.05, y = 5, yend = 9)+
  annotate("text", x = 0.125, y = -3, label = "0.25 s", vjust = 1.3, size = textsize/.pt)+
  annotate("text", x = -0.05, y = 7, label = "4 pN", vjust = -0.4 , angle = 90, size =  textsize/.pt)+
  annotate("text", x = 1.22, y = 3.5, label = "M", hjust = 0 , vjust = 0 , angle = 0, size = textsize/.pt)+
  annotate("text", x = 1.22, y = -0.27, label = "T", hjust = 0 , vjust = 0 , angle = 0, size = textsize/.pt)+
  coord_cartesian(xlim = c(-0.13, 1.3), ylim = c(-3.7, NA))+
  theme_void(basesize)+
     theme(
       plot.title = element_text(face = "bold", hjust = 0.5, color = colorz[[2]]))
)

reg_raw <- reg_raw + ggtitle(paste0("Regulated (n = ", reg_n, ")"))

fig4_top <- plot_grid(unreg_raw, reg_raw, nrow = 1)

png(filename = "img/figure-4_feedback.png", width = 6.75, height = 4.5, units = "in", res = 300)
plot_grid(fig4_top, fig4_bottom, nrow = 2, rel_heights = c(0.4, 0.6), labels = letters)
dev.off()
svg(filename = "img/figure-4_feedback.svg", width = 6.75, height = 4.5)
plot_grid(fig4_top, fig4_bottom, nrow = 2, rel_heights = c(0.4, 0.6), labels = letters)
dev.off()
cairo_ps(filename = "img/figure-4_feedback.eps", width = 6.75, height = 4.5)
plot_grid(fig4_top, fig4_bottom, nrow = 2, rel_heights = c(0.4, 0.6), labels = letters)
dev.off()

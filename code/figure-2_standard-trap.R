library(data.table)
library(dplyr)
library(purrr)
library(readxl)
library(here)
library(cowplot)
library(ggtext)
library(ggpubr)

# define global text sizes for plots and tables
textsize <- 9
basesize <- 11
# get file locations for the SPASM output files
# these files contain all the data for the analyzed trap events
files <- list(
          Unregulated = here("data", "standard-trap", "unregulated_1uM-ATP_standard", "spasm-unregulated-standard.xlsx"),
          Regulated = here("data", "standard-trap" ,"regulated_1uM-ATP_standard", "spasm-regulated-standard.xlsx")
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

spasm_data$id <- factor(spasm_data$id, levels = c("Unregulated", "Regulated"))

# make new columns based on avg of 2 beads position / rename a few
spasm_data <- mutate(spasm_data,
                     "Displacements" = (`A total step (nm)` + `B total step (nm)`)/2,
                     "Attachment Times" = `duration (s)`,
                     "Substep 1" = (`A step 1 (nm)` + `B step 1 (nm)`) / 2,
                     "Substep 2" = (`A step 2 (nm)` + `B step 2 (nm)`) / 2)


spasm_sum <- spasm_data |> group_by(id) |> summarize(n = n())

n_unreg <- spasm_sum$n[spasm_sum$id == "Unregulated"]
n_reg <- spasm_sum$n[spasm_sum$id == "Regulated"]

#### Raw traces ###
colz <- c(unname(palette.colors())[c(1, 4)])

# helper function to read in data from the greenberg trap
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

# get file location
raw_unregulated_file <- here("data",
                             "standard-trap",
                             "unregulated_1uM-ATP_standard",
                             "2022-08-17",
                             "obs-01",
                             "220817 2.txt")

# read in raw data
raw_unreg_trace <- read_greenberg(raw_unregulated_file)

# obtains data point index to start and stop plotting of raw trace
start <- 101*raw_unreg_trace$options$hz
stop <- 104.5*raw_unreg_trace$options$hz

raw_unreg <- raw_unreg_trace$data[start:stop, .(Trap2X,
                                                dp = .I)]

# read it spasm data for event indentification bars
unreg_spasm <- read_xlsx(here("data",
                             "standard-trap",
                             "unregulated_1uM-ATP_standard",
                             "spasm-unregulated-standard.xlsx"),
                       sheet = 9,
                       skip = 19
                       ) |>
  filter(`index before start (# pts)` >= start & `index before end (# pts)` <= stop) |>
  mutate(starter = `index before start (# pts)` - start,
         stopper = `index before end (# pts)` - start)

raw_unreg[, y := Trap2X * -1 * raw_unreg_trace$options$nm_v2]


# make the plot and add scale bars with geom_segment
gg_raw_unreg <-
  ggplot(raw_unreg)+
  geom_line(aes(dp, y), linewidth = 0.15)+
   draw_line(x = c(-2000, 18000),
             y = c(rep(min(raw_unreg$y), 2)),
             linewidth = 0.7)+
   draw_line(x = rep(-2000, 2),
             y = c(mean(raw_unreg$y)+20,
                   mean(raw_unreg$y)-20),
             linewidth = 0.7)+
   geom_segment(data = unreg_spasm,
                aes(x = starter, xend = stopper,
                    y = -90, yend = -90),
                color = "black",
                alpha = 0.35,
                linewidth = 1)+
  annotate("text", x = -2500, y = mean(raw_unreg$y), label = "40 nm", angle = 90, vjust = -0, size = textsize/.pt)+
  annotate("text", x = 8000, y = min(raw_unreg$y)-1, label = "1 s", angle = 0, vjust = 1, size = textsize/.pt)+
  coord_cartesian(ylim = c(min(raw_unreg$y)-6, NA ))+
  ggtitle(paste0("Unregulated (n = ", n_unreg, ")"))+
  ylab("")+
  xlab("")+
  theme_void(basesize)+
  theme(
  plot.title = element_text(hjust = 0.50, face = "bold", size = basesize)
  )

# do the same as above, but for regulated
# probably would be better to define a function vs copy/paste...
raw_regulated_file <- here("data",
                           "standard-trap",
                           "regulated_1uM-ATP_standard",
                           "221208-f2-1um-atp-m3 15.txt")

raw_reg_trace <- read_greenberg(raw_regulated_file)

start2 <- 3.25*raw_reg_trace$options$hz
stop2 <- 6.9*raw_reg_trace$options$hz

raw_reg <- raw_reg_trace$data[start2:stop2, .(Trap2X,
                                                dp = .I)]

reg_spasm <- read_xlsx(here("data",
                            "standard-trap",
                            "regulated_1uM-ATP_standard",
                            "spasm-regulated-standard.xlsx"),
                       sheet = 1,
                       skip = 19
                       ) |>
  filter(`index before start (# pts)` >= start2 & `index before end (# pts)` <= stop2) |>
  mutate(start = `index before start (# pts)` - start2,
         stop = `index before end (# pts)` - start2)


raw_reg[, y := Trap2X * raw_reg_trace$options$nm_v2]

gg_raw_reg <-
  ggplot(raw_reg)+
  geom_line(aes(x = dp,
                y = (Trap2X * raw_reg_trace$options$nm_v2)),
            color = colz[2],
            linewidth = 0.15)+
   draw_line(x = rep(-2000, 2),
             y = c(mean(raw_reg$y)+20,
                   mean(raw_reg$y)-20),
             linewidth = 0.7)+
   draw_line(x = c(-2000, 18000),
             y = rep(min(raw_reg$y), 2),
             linewidth = 0.7)+
   geom_segment(data = reg_spasm,
                aes(x = start, xend = stop,
                    y = 185, yend = 185),
                color = "black",
                alpha = 0.35,
                linewidth = 1)+
  annotate("text",
           x = -2500,
           y = mean(raw_reg$y),
           label = "40 nm",
           angle = 90,
           vjust = -0,
           size = textsize/.pt)+
  annotate("text",
           x = 8000,
           y = min(raw_reg$y)-1,
           label = "1 s",
           angle = 0,
           vjust = 1,
           size = textsize/.pt)+
  coord_cartesian(ylim = c(min(raw_reg$y)-6, NA))+
  ggtitle(paste0("Regulated (n = ", n_reg, ")"))+
  ylab("")+
  xlab("")+
 theme_void(basesize)+
  theme(
    plot.title = element_text(hjust = 0.50, face = "bold", size = basesize, color = colz[2])
  )


# combine plots
gg_raw_plots <- plot_grid(gg_raw_unreg, gg_raw_reg, nrow = 1, labels = c("a", "b"))


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


# fit each conditions ensemble average
g1 <-
  spasm_plot_ensemble_average(here("data", "standard-trap", "ensemble-averages", "unregulated-combinedEnsembleAxesData.xlsx"),
                                 color = alpha(colz[1], 0.5),
                                 x_shift = 1.2,
                                 title = "Unregulated")

g2 <-
  spasm_plot_ensemble_average(here("data", "standard-trap", "ensemble-averages", "regulated-combinedEnsembleAxesData.xlsx"),
                            color = colz[2],
                            x_shift = 1.2,
                            title = "Regulated")


# combine enemble average plots
gg_ea <- plot_grid(g1$gg, g2$gg, nrow = 1, labels = c("c", "d"))

#### Step ####

# plot the displacements and substeps
# helper function to expedite the base plot
spasm_plot_ecdf <- function(spasm_data, var, colz, x_lab = "x_label"){

  ggplot(spasm_data)+
    stat_ecdf(aes(base::get(var),
                  color = id),
              linewidth = 1,
              pad = FALSE,
              show.legend = FALSE)+
    scale_color_manual(values = colz)+
    ylab("Cumulative Distribution")+
    xlab(x_lab)+
    ggtitle(var)+
    theme_cowplot(basesize)+
    theme(
      plot.title = element_text(hjust = 0.5, size = basesize)
    )
}

spasm_unreg <- filter(spasm_data, id == "Unregulated")
spasm_reg <- filter(spasm_data, id == "Regulated")


total_step_pvalue <- t.test(spasm_unreg$Displacements, spasm_reg$Displacements)$p.value


total_step_ecdf <- spasm_plot_ecdf(spasm_data, var = "Displacements", colz, x_lab = "nanometers")+
  ggtitle("Total Displacements")+
  ylab("Cumulative Probability")+
  xlab("Displacement (nm)")+
   annotate("text",
            label = paste0("P = ", round(total_step_pvalue, 2)),
            x = Inf,
            y = -Inf,
            hjust = 1,
            vjust = -0.15,
            size = textsize/.pt)+
 annotate("text",
          x = -Inf,
          y = 1,
          label = paste0("~bar(x)==",
                        round(mean(spasm_unreg$Displacements), 1),
                        "%+-%", round(sd(spasm_unreg$Displacements), 1),
                        "~nm"),
                        ## "<br>",
                        ## "<span style = 'color: #009E73;'>",
                        ## "x&#772; = ", round(mean(spasm_reg$Displacements), 1),
                        ## " &plusmn; ", round(sd(spasm_reg$Displacements), 1),
                        ## " nm"),
          parse = TRUE,
          hjust = 0,
          vjust = 1,
          size = textsize/.pt)+
 annotate("text",
          x = -Inf,
          y = 0.9,
          label = paste0("~bar(x)==",
                        round(mean(spasm_reg$Displacements), 1),
                        "%+-%", round(sd(spasm_reg$Displacements), 1),
                        "~nm"),
          color = "#009E73",
          parse = TRUE,
          hjust = 0,
          vjust = 1,
          size = textsize/.pt)+
   coord_cartesian(xlim = c(-32, NA))+
   theme(axis.title.x = element_text(color = "white"),
         axis.title.y = element_text(size = textsize+1))



# shapiro wilk tests
unreg_to_test_norm <- list(spasm_unreg$Displacements,
                           spasm_unreg$`Substep 1`,
                           spasm_unreg$`Substep 2`)

lapply(unreg_to_test_norm, shapiro.test)


reg_to_test_norm <- list(spasm_reg$Displacements,
                           spasm_reg$`Substep 1`,
                           spasm_reg$`Substep 2`)

lapply(reg_to_test_norm, shapiro.test)


one_step_pvalue <- t.test(spasm_unreg$`Substep 1`, spasm_reg$`Substep 1`)$p.value

one_step_ecdf <-
  spasm_plot_ecdf(spasm_data, var = "Substep 1", colz, x_lab = "nanometers")+
  ylab("Cumulative Probability")+
  xlab("Displacement (nm)")+
   annotate("text",
            label = paste0("P = ", round(one_step_pvalue, 2)),
            x = Inf,
            y = -Inf,
            hjust = 1,
            vjust = -0.15,
            size = textsize/.pt)+
 annotate("text",
          x = -Inf,
          y = 1,
          label = paste0("bar(x)==",
                        round(mean(spasm_unreg$`Substep 1`), 1),
                        "%+-%", round(sd(spasm_unreg$`Substep 1`), 1),
                        "~nm"),
          parse = TRUE,
          hjust = 0,
          vjust = 1,
          size = textsize/.pt)+
 annotate("text",
          x = -Inf,
          y = 0.9,
          label = paste0("bar(x)==",
                        round(mean(spasm_reg$`Substep 1`), 1),
                        "%+-%", round(sd(spasm_reg$`Substep 1`), 1),
                        "~nm"),
          color = "#009E73",
          parse = TRUE,
          hjust = 0,
          vjust = 1,
          size = textsize/.pt)+
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank())


two_step_pvalue <- t.test(spasm_unreg$`Substep 2`, spasm_reg$`Substep 2`)$p.value

two_step_ecdf <-
  spasm_plot_ecdf(spasm_data, var = "Substep 2", colz, x_lab = "nanometers")+
  ## ggtitle("Total Displacements")+
  ylab("Cumulative Probability")+
  xlab("Displacement (nm)")+
   annotate("text",
            label = paste0("P = ", round(two_step_pvalue, 2)),
            x = Inf,
            y = -Inf,
            hjust = 1,
            vjust = -0.15,
            size = textsize/.pt)+
 annotate("text",
          x = -Inf,
          y = 1,
          label = paste0("bar(x)==",
                        round(mean(spasm_unreg$`Substep 2`), 1),
                        "%+-%", round(sd(spasm_unreg$`Substep 2`), 1),
                        "~nm"),
          parse = TRUE,
          hjust = 0,
          vjust = 1,
          size = textsize/.pt)+
   annotate("text",
          x = -Inf,
          y = 0.9,
          label = paste0("bar(x)==",
                        round(mean(spasm_reg$`Substep 2`), 1),
                        "%+-%", round(sd(spasm_reg$`Substep 2`), 1),
                        "~nm"),
          color = "#009E73",
          parse = TRUE,
          hjust = 0,
          vjust = 1,
          size = textsize/.pt)+
   coord_cartesian(xlim = c(-25, NA))+
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(color = "white"))


# combine step plots, figure 2 bottom row
gg_steps <- plot_grid(total_step_ecdf,
                       one_step_ecdf,
                       two_step_ecdf,
                       nrow = 1,
                       rel_widths = c(1, 0.9, 0.9),
                       labels = c("e", "f", "g"))


#combine all plots
## plot_grid(gg_raw_plots, gg_ea, gg_steps, rel_heights = c(0.25, 0.5, 0.4), nrow = 3)

# save the plot in specified file
# combine all the plots into one to make figure 2
png(filename = "img/figure-2_standard-trap.png", width = 6.5, height = 6.5, units = "in", res = 300)
plot_grid(gg_raw_plots, gg_ea, gg_steps, rel_heights = c(0.25, 0.5, 0.4), nrow = 3)
dev.off()
svg(filename = "img/figure-2_standard-trap.svg", width = 6.5, height = 6.5)
plot_grid(gg_raw_plots, gg_ea, gg_steps, rel_heights = c(0.25, 0.5, 0.4), nrow = 3)
dev.off()
cairo_ps(filename = "img/figure-2_standard-trap.eps", width = 6.5, height = 6.5)
plot_grid(gg_raw_plots, gg_ea, gg_steps, rel_heights = c(0.25, 0.5, 0.4), nrow = 3)
dev.off()


# make a figure version for conference poster


## gg_steps2 <- plot_grid(total_step_ecdf,
##                        one_step_ecdf,
##                        two_step_ecdf,
##                        nrow = 1,
##                        rel_widths = c(1, 0.9, 0.9),
##                        labels = c("c", "d", "e"))
## #for poster
## png(file = "code/poster/img/step-size.png", width = 6.5, height = 4, units = "in", res = 500)
## plot_grid(gg_raw_plots,gg_steps2, rel_heights = c(0.25, 0.4), nrow = 2)
## dev.off()

## png(file = "code/poster/img/ea.png", width = 6.5, height = 3, units = "in", res = 500)
## plot_grid(g1$gg, g2$gg, nrow = 1, labels = c("b", "c"))
## dev.off()

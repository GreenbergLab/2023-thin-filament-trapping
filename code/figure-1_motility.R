# this is the only figure that is not 100% reprodicble in R alone.
# parts of the figure were already completed by co-authors

# this will generate the motility plot/stats
library(readxl)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggtext)
## library(ggpubr)


# helper function to read in motility track data
read_mot <- function(x){
  d <- fread(x)
  d$path <- x
  d
}

# find the file paths, read in, and combine into one dataframe
data_paths <- list.files("data/motility",
                  pattern = "MTrackJ-Tracks.csv",
                  recursive = T,
                  full.names = T)


dat <-  rbindlist(lapply(data_paths, read_mot))

# split the file path into new columns to get out conditions metadata
dat[, c("data", "mot", "actin", "pCa", "date", "video", "file") :=
        tstrsplit(path, "/", fixed = TRUE)]


# make the labels prettier by replacing dashes with spaces
dat$pCa_pretty <- gsub(x = dat$pCa, pattern = "-", replacement = " ")

# ordering conditions so they will plot in the desired order
dat$pCa_pretty <- factor(dat$pCa_pretty, levels = c("pCa 9", "pCa 6.25", "pCa 4"))
dat$actin <- factor(dat$actin, levels = c("regular-actin", "biotin-actin"))

# getting colors
colz <- c("#ff4d4dff", "#c4b100ff", unname(palette.colors())[4])

# making a label in HTML for use with ggtext pacakge so condition can be written in its color in the facet header
## dat[, pCa_pretty_html := paste0("<span style='color:", colz, "'>", pCa_pretty, "</span>")]

# summarizing data with average calculations
avg_by_date <- dat[, .(mean_v_nm_s = mean(`Mean v [nm/sec]`, na.rm = TRUE)), by  = .(actin, pCa, date)]

avg_by_video <- dat[, .(mean_v_nm_s = mean(`Mean v [nm/sec]`, na.rm = TRUE)), by  = .(actin, pCa, date, video)]

avg_by_pCa <- dat[, .(mean_v_nm_s = mean(`Mean v [nm/sec]`, na.rm = TRUE),
                      sd = sd(`Mean v [nm/sec]`, na.rm = TRUE),
                      n = .N), by  = .(actin, pCa_pretty)]

avg_by_pCa[, colz := ifelse(pCa_pretty == "pCa 9", colz[1],
                          ifelse(pCa_pretty == "pCa 6.25", colz[2],
                                 ifelse(pCa_pretty == "pCa 4", colz[3], "error")))]

avg_by_pCa[, pCa_pretty_html := paste0("<span style='color:", colz, "'>", pCa_pretty, "</span>")]

avg_by_pCa$pCa_pretty_html <- factor(avg_by_pCa$pCa_pretty_html, levels = c(
                                     "<span style='color:#ff4d4dff'>pCa 9</span>",
                                     "<span style='color:#c4b100ff'>pCa 6.25</span>",
                                     "<span style='color:#009E73'>pCa 4</span>"))

ggplot()+
  geom_errorbar(data = avg_by_pCa,
                 aes(x = actin,
                     y = mean_v_nm_s,
                     ymax = mean_v_nm_s+sd,
                     ymin = mean_v_nm_s-sd),
                width = 0.25,
           color = "black",
           show.legend = FALSE)+
  geom_col(data = avg_by_pCa,
           aes(x = actin,
               y = mean_v_nm_s,
               fill = pCa_pretty_html),
           color = "black",
           show.legend = FALSE)+
  facet_wrap(~pCa_pretty_html)+
  scale_x_discrete(labels = c("regular-actin"="Actin",
                              "biotin-actin"="Actin+Biotin"))+
  scale_y_continuous(expand = expansion(c(0, 0.1), c(0, 0.1)))+
  coord_cartesian(ylim = c(0, NA))+
  xlab("")+
  ylab("Speed (&#xb5;m/s)")+
  scale_color_manual(values = colz)+
  scale_fill_manual(values = colz)+
  scale_alpha_manual(values = c(0.5, 1))+
  theme_cowplot(18)+
  theme(
    axis.text.x = element_text(size = 10),
    ## axis.ticks.x = element_blank(),
    strip.text = element_markdown(),
    strip.background = element_blank(),
    axis.title.y = element_markdown(),
    legend.position = "none")

# save the plot
## ggsave("img/motility-pca.png", dpi = 500, bg = "white")

## Stats - uncomment sink lines to save to a file instead of printing to console
## sink("data/motility/motility-p-values.txt")
cat("This file was generated from the R script located in the project directory at ./code/figure-1_motility.R.
The contents of the file contain the p-values for the comparisons of motility velocities within a pCa.
There are two groups, actin vs biotin-actin.
The results follow: \n
pCa 4
\n")
wilcox.test(`Mean v [nm/sec]` ~ actin, data = dat[pCa == "pCa-4"])
cat("\n pCa 6.25 \n")
wilcox.test(`Mean v [nm/sec]` ~ actin, data = dat[pCa == "pCa-6.25"])
cat("\n pCa 9 \n")
wilcox.test(`Mean v [nm/sec]` ~ actin, data = dat[pCa == "pCa-9"])
## sink()


## fwrite(avg_by_pCa, file = "data/motility/motility-summary-values.csv")

#################################################
# pCa 6.25 trap trace
#################################################


# defines a helper function to help read in SPASM trap data
read_spasm_workbook <- function(spasm_file, id){
   map_df(excel_sheets(spasm_file),
           ~read_excel(spasm_file,
                       sheet = .x,
                       skip = 19)
        )|>
    mutate(id = id)

}

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
pca_625_file <- here("data",
                     "standard-trap",
                     "pCa-6.25",
                     "2023-03-21",
                     "obs-53",
                     "230321_f2_s3_m3 53.txt")

# read in raw data
pca_625_trace <- read_greenberg(pca_625_file)

# obtains data point index to start and stop plotting of raw trace
start <- 1*pca_625_trace$options$hz
stop <- 20*pca_625_trace$options$hz

raw_625 <- pca_625_trace$data[start:stop, .(Trap2X,
                                      dp = .I)]

# read it spasm data for event indentification bars
spasm_625 <- read_xlsx(here("data",
                             "standard-trap",
                             "pCa-6.25",
                             "2023-03-21",
                            "obs-53",
                            "spasm.xlsx"),
                       skip = 19
                       ) |>
  dplyr::filter(`index before start (# pts)` >= start & `index before end (# pts)` <= stop) |>
  dplyr::mutate(starter = `index before start (# pts)` - start,
         stopper = `index before end (# pts)` - start)

raw_625[, y := Trap2X * pca_625_trace$options$nm_v2]


# make the plot and add scale bars with geom_segment
(gg_625 <-
  ggplot(raw_625)+
  geom_line(aes(dp, y), linewidth = 0.45, color = "#c4b100")+
   draw_line(x = c(-2000, 98000),
             y = c(rep(min(raw_625$y), 2)),
             linewidth = 0.7)+
   draw_line(x = rep(-2000, 2),
             y = c(mean(raw_625$y)+25,
                   mean(raw_625$y)-25),
             linewidth = 0.7)+
   geom_segment(data = spasm_625,
                aes(x = starter, xend = stopper,
                    y = 105, yend = 105),
                color = "black",
                alpha = 0.35,
                linewidth = 1)+
  annotate("text", x = -2800, y = mean(raw_625$y), label = "40 nm", angle = 90, vjust = -0, size = 11/.pt)+
  annotate("text", x = 50000, y = min(raw_625$y)-1, label = "5 s", angle = 0, vjust = 1, size = 11/.pt)+
  coord_cartesian(ylim = c(min(raw_625$y)-6, NA ))+
  ggtitle("pCa 6.25")+
  ylab("")+
  xlab("")+
  theme_void(12)+
  theme(
  plot.title = element_text(hjust = 0.50, face = "bold", size = 12)
  )
)

# uncomment to save the data
## ggsave("img/pca-625-trace.png", dpi = 500)


######################################################################
# prior to the review - old #
######################################################################

# read in motility data and transform to long format got ggplot
## motility_file <- here("data/motility/WT & Biotin Regulated Motility 200nM Tm&Tn.xlsx")
## motility_data <- read_xlsx(motility_file, sheet = 4, skip = 1)
## motility_data <- melt(setDT(motility_data),
##                       measure.vars = c("regular-actin", "biotin-actin"),
##                       variable.name = "id",
##                      value.name = "velocity")

## motility_data$id <- factor(motility_data$id, levels = c("regular-actin", "biotin-actin"))

## motility_summary <- motility_data[, .(velocity_avg = mean(velocity, na.rm = TRUE),
##                                       velocity_sd = sd(velocity, na.rm = TRUE)),
##                                   by = id]


## # check normality
## shapiro.test(motility_data$velocity)

## # it fails - use non-parametrics
## t_test_results <- compare_means(velocity~id, data = motility_data, method = "wilcox.test")

## colorz <- c("black", "grey50")

## # Motility plot
## (
## gg_mot <-
## ggplot(data = motility_data)+
##   geom_boxplot(data = motility_data,
##                 aes(x = id,
##                     y = velocity,
##                      fill = id),
##                outlier.color = NA,
##                size = 1,
##                alpha = 0.5,
##                 width = 0.3)+
##   geom_jitter(data = motility_data,
##               aes(x = id,
##                   y = velocity,
##                   color = id),
##                width = 0.1,
##               size =  1.5,
##               alpha = 0.5,
##               shape = 16)+
##   stat_pvalue_manual(data = t_test_results, label = "p = {p.adj}", y.position=0.62, tip.length = 0, label.size = 7)+
##   ## scale_y_continuous(expand = expansion(c(0, 0.1)))+
##   scale_x_discrete(labels = c("regular-actin"="Actin", "biotin-actin"="Actin+Biotin"))+
##   scale_color_manual(values = colorz)+
##   scale_fill_manual(values = colorz)+
##   coord_cartesian(ylim = c(-0.05, 0.65))+
##   xlab("")+
##   ylab("Velocity (&#xb5;m/s)")+
##   ggtitle("Regulated Motility pCa 4")+
##   theme_cowplot(22)+
##   theme(
##     legend.position = "none",
##     axis.title.y = element_markdown()
## )
## )

## ## ggsave(here("code/poster/img", "motility-pca-4.png"), dpi = 500, bg = "white")

## # pca4 reg trace
## gg_pca4 <- readRDS(here("img", "pca4-trace.rds"))
## ## saveRDS(gg_pca4, here("img/pca4-trace.rds"))
## ## ggsave(here("img/pca4-trace.png"), dpi = 500)

## ## after generating for figures, I imported to inkscape to combine with
## ## an existing .ai file I inherited for this project containing a pCa 9 trace.

## # for poster

## pca9 <- data.frame(id = c("regular-actin-9", "biotin-actin-9"),
##                    velocity  = c(0, 0))


## pca9$id <- factor(pca9$id, levels = c("regular-actin-9", "biotin-actin-9"))

## all_mot <- rbind(pca9, motility_data)

## all_mot$id <- factor(all_mot$id, levels = c("regular-actin-9", "regular-actin", "biotin-actin-9", "biotin-actin"))
## ## ggplot(pca9)+
## ##   geom_point(aes(x = id,
## ##                  y = velocity),
## ##              shape = 4,
## ##              color = "red",
## ##              size = 10)+
## ##   coord_cartesian(ylim = c(-0.05, 0.65))+
## ##   scale_x_discrete(labels = c("regular-actin"="Actin", "biotin-actin"="Actin+Biotin"))+
## ##   xlab("")+
## ##   ylab("Velocity (&#xb5;m/s)")+
## ##   ggtitle("Regulated Motility pCa 9")+
## ##   theme_cowplot(22)+
## ##   theme(
## ##     legend.position = "none",
## ##     axis.title.y = element_markdown()
## ## )

## (
## gg_mot2 <-
## ggplot(data = all_mot)+
##   geom_boxplot(data = all_mot,
##                 aes(x = id,
##                     y = velocity,
##                      fill = id),
##                outlier.color = NA,
##                size = 1,
##                alpha = 0.5,
##                 width = 0.3)+
##   geom_jitter(data = motility_data,
##               aes(x = id,
##                   y = velocity,
##                   color = id),
##                width = 0.1,
##               size =  1.5,
##               alpha = 0.5,
##               shape = 16)+
##   ## geom_point(data = pca9,
##   ##            aes(x = id,
##   ##                y = (velocity+0.04)),
##   ##            shape = 4,
##   ##            color = "red",
##   ##            size = 10)+
##   stat_pvalue_manual(data = t_test_results, label = "p = {p.adj}", y.position=0.62, tip.length = 0, label.size = 7)+
##   ## scale_y_continuous(expand = expansion(c(0, 0.1)))+
##   scale_x_discrete(labels = c("regular-actin"="Actin (pCa 4)",
##                               "biotin-actin"="Actin+Biotin (pCa 4)",
##                               "regular-actin-9" = "Actin (pCa 9)",
##                               "biotin-actin-9" = "Actin+Biotin (pCa 9)"))+
##   scale_y_continuous(expand = expansion(mult = c(0, NA), add = c(0, NA)))+
##   scale_color_manual(values = rep(colorz, 2))+
##   scale_fill_manual(values = rep(colorz, 2))+
##   coord_cartesian(ylim = c(-0.001, 0.68))+
##   xlab("")+
##   ylab("Velocity (&#xb5;m/s)")+
##   ggtitle("Regulated Motility")+
##   theme_cowplot(30)+
##   theme(
##     legend.position = "none",
##     axis.title.y = element_markdown(),
##     plot.title = element_text(hjust = 0.5),
##     axis.text.x = element_text(size = 20)
## )
## )


## ggsave(here("code/poster/img", "motility-pca9.png"), dpi = 500, bg = "white")

library(data.table)
library(dygraphs)
library(ggplot2)
library(pracma)
library(zoo)

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


plot_greenberg_raw_data <- function(x, units){


  if(units == "volts"){
   b1 <- x$data$Trap1X
   b2 <- x$data$Trap2X
  } else if(units == "nm"){
   b1 <- x$data$Trap1X * x$options$nm_v1
   b2 <- x$data$Trap2X * x$options$nm_v2
  } else if(units == "pN"){
   b1 <- x$data$Trap1X * x$options$nm_v1 * x$options$pn_nm1
   b2 <- x$data$Trap2X * x$options$nm_v2 * x$options$pn_nm2
  }

  d <- data.table(x = 1:nrow(x$data)/x$options[["hz"]],
                  bead1 = b1,
                  bead2 = b2)

  dygraphs::dygraph(data = d) |> dyRangeSelector(fillColor = "")

}

files <- list.files("../brent-data/unregulated_1uM-ATP_standard",
                   pattern = ".txt",
                   full.names = TRUE,
                   recursive = TRUE)

trap_data <- read_greenberg(files[10])

## plot_greenberg_raw_data(trap_data, "nm")

bp <- seq(20000, nrow(trap_data$data), by = 100000)

t1 <- detrend((trap_data$data$Trap1X * trap_data$options$nm_v1), tt = "linear", bp = bp)
t2 <- detrend(trap_data$data$Trap2X * trap_data$options$nm_v2, tt = "linear", bp = bp)


xrange <- 40000:130000
t1_<-t1[xrange]
t2_<-t2[xrange]


covar <- ksmooth(1:length(t1_), rollapply(data.frame(t1_, t2_), width = 200, FUN = \(x) cov(x[,1], x[,2]), by.column = FALSE), bandwidth = 500)

ggplot()+
  geom_line(aes(x = 1:length(xrange)/20000,
                y = t1_),
            color = "#6c7373")+
  geom_line(aes(x = 1:length(xrange)/20000,
                y = t2_+90),
            color = "#a51417")+
  geom_line(aes(1:length(covar$y)/20000,
                y = covar$y-155),
            color = "#007360")+
  annotate("text",
           x = 0,
           y = 100,
           label = "Trap 1",
           hjust= 1.1)+
  annotate("text",
           x = 0,
           y = 0,
           label = "Trap 2",
           hjust = 1.1)+
  annotate("text",
           x = 0,
           y = -100,
           label = "Covariance",
           hjust = 1.1)+
  scale_x_continuous(expand = expansion(add = c(0.8, 0.1)))+
  xlab("seconds")+
  ylab("nanometers")+
  theme_minimal()

## ggsave("img/ggcovar.pdf")

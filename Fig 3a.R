library(readr)
library(deeptime)
library(ggplot2)
library(egg)
library(fANCOVA)
library(dplyr)
library(plyr)
library(reshape2)

# Set working directory (this will likely be the downloaded folder "Mo.U.Model")

setwd()

##==============================================
## NB - deeptime package required to reproduce these figures - if not yet installed need to use devtools to obtain from GitHub
## If not installed, comment out and run the lines of code directly below - install devtools if not installed already
##==============================================
# install.packages("devtools")
# devtools::install_github("willgearty/deeptime")

##==============================================
## Edit Deeptime timescale to fit 2012 geologic timescale and center labels wrt to x-axis
##==============================================

adapted_period_timescale <-  deeptime::periods
adapted_period_timescale[10,2] <- 443.8
adapted_period_timescale[11,3] <- 443.8
adapted_period_timescale[10,3] <- 440.8
adapted_period_timescale[11,2] <- 447.5

adapted_stages_timescale <-  deeptime::stages
adapted_stages_timescale[83,2] <- 443.8
adapted_stages_timescale[84,3] <- 443.8
adapted_stages_timescale[85,2] <- 447.5


##==============================================
## 3A - feux
##==============================================

round.factor <- 10

## Import Bartlett (with reconstructed age model)
load("Bartlett.full.RData")

## Omit any NA values
Bartlett <- filter(Bartlett, !is.na(U_iso))

## Fit LOESS regression to Bartlett
U.carb.fit <- loess.as(Bartlett$Age_est, Bartlett$U_iso,  criterion = c("aicc", "gcv"))

## Predict time-dependent geochemical data from Bartlett data
U.carb.model <- predict(U.carb.fit, seq(444.7,443.2, -0.02))

## Generate rounding factor
U.carb.round <- (range(Bartlett$U_iso)[2]-range(Bartlett$U_iso)[1])/round.factor

## Import Murzuq Basin data (with age model)
load("Murzuq.RData")

## Fit cross-validated LOESS regression to Murzuq U and Mo
Murzuq <- filter(Murzuq, !is.na(d98Mo_auth_Al))
Murzuq <- filter(Murzuq, !is.na(dU238_auth_a))

Mo.shale.fit <- loess.as(Murzuq$Age_est, Murzuq$d98Mo_auth_Al,  criterion = c("aicc"))
U.shale.fit <- loess.as(Murzuq$Age_est, Murzuq$dU238_auth_a,  criterion = c("aicc"))

## Predict time-dependent geochemical data from Murzuq data
Mo.shale.model <- predict(Mo.shale.fit, seq(443.8,440.8, -0.02))
U.shale.model <- predict(U.shale.fit, seq(443.8,440.8, -0.02))

## Generate rounding factors
U.shale.round <- (range(Murzuq$dU238_auth_a)[2]-range(Murzuq$dU238_auth_a)[1])/round.factor
Mo.shale.round <- (range(Murzuq$d98Mo_auth_Al)[2]-range(Murzuq$d98Mo_auth_Al)[1])/round.factor

## Import combined shale and carbonate model
load("Mo-U.mass.balance.RData")

## Round model
mass.balance$dU.carb <- plyr::round_any(mass.balance$dU.carb, U.carb.round, f=round)
mass.balance$dMo <- plyr::round_any(mass.balance$dMo, Mo.shale.round, f=round)
mass.balance$dU <- plyr::round_any(mass.balance$dU, U.shale.round, f=round)

## Initiate data frame for summary
f.eux.sum <- data.frame(timebin = double(), mean.f.eux.value = double(), 
                        median.f.eux.value = double(),
                        min.f.eux.value.05 = double(),
                        max.f.eux.value.95 = double(),
                        min.f.eux.value.10 = double(),
                        max.f.eux.value.90 = double(),
                        min.f.eux.value.20 = double(),
                        max.f.eux.value.80 = double(),
                        min.f.eux.value.30 = double(),
                        max.f.eux.value.70 = double(),
                        min.f.eux.value.40 = double(),
                        max.f.eux.value.60 = double())

## Round predictions to fit model rounding
## At 20kyr timesteps filter models based on LOESS models of geochemnical data for each percentile presented (plus some other statistics of interest)

# Bartlett data 
for(timebin in seq(444.7, 443.2, -0.02)){
  
  U.carb.value <- predict(U.carb.fit, timebin)
  U.carb.value.rounded <- plyr::round_any(U.carb.value, U.carb.round, f=round)
  f.eux.values <- filter(mass.balance, dU.carb == U.carb.value.rounded)
  mean.f.eux.value <- mean(f.eux.values$f.eux)
  median.f.eux.value <- median(f.eux.values$f.eux)

  min.f.eux.value.05 <- quantile(f.eux.values$f.eux, probs=0.05)
  max.f.eux.value.95 <- quantile(f.eux.values$f.eux, probs=0.95)
  min.f.eux.value.10 <- quantile(f.eux.values$f.eux, probs=0.1)
  max.f.eux.value.90 <- quantile(f.eux.values$f.eux, probs=0.9)
  min.f.eux.value.20 <- quantile(f.eux.values$f.eux, probs=0.2)
  max.f.eux.value.80 <- quantile(f.eux.values$f.eux, probs=0.8)
  min.f.eux.value.30 <- quantile(f.eux.values$f.eux, probs=0.3)
  max.f.eux.value.70 <- quantile(f.eux.values$f.eux, probs=0.7)
  min.f.eux.value.40 <- quantile(f.eux.values$f.eux, probs=0.4)
  max.f.eux.value.60 <- quantile(f.eux.values$f.eux, probs=0.6)
  
  f.eux.sum <- rbind(f.eux.sum, c(timebin, mean.f.eux.value, 
                                  median.f.eux.value,
                                  min.f.eux.value.05,
                                  max.f.eux.value.95,
                                  min.f.eux.value.10,
                                  max.f.eux.value.90,
                                  min.f.eux.value.20,
                                  max.f.eux.value.80,
                                  min.f.eux.value.30,
                                  max.f.eux.value.70,
                                  min.f.eux.value.40,
                                  max.f.eux.value.60))
  
}

# Murzuq Basin data 
for(timebin in seq(443.8, 440.8, -0.02)){
  U.shale.value <- predict(U.shale.fit, timebin)
  U.shale.value.rounded <- plyr::round_any(U.shale.value, U.shale.round, f=round)
  Mo.shale.value <- predict(Mo.shale.fit, timebin)
  Mo.shale.value.rounded <- plyr::round_any(Mo.shale.value, Mo.shale.round, f=round)
  f.eux.values <- filter(mass.balance, dU == U.shale.value.rounded 
                         & dMo == Mo.shale.value.rounded)
  mean.f.eux.value <- mean(f.eux.values$f.eux)
  median.f.eux.value <- median(f.eux.values$f.eux)
  min.f.eux.value.05 <- quantile(f.eux.values$f.eux, probs=0.05)
  max.f.eux.value.95 <- quantile(f.eux.values$f.eux, probs=0.95)
  min.f.eux.value.10 <- quantile(f.eux.values$f.eux, probs=0.1)
  max.f.eux.value.90 <- quantile(f.eux.values$f.eux, probs=0.9)
  min.f.eux.value.20 <- quantile(f.eux.values$f.eux, probs=0.2)
  max.f.eux.value.80 <- quantile(f.eux.values$f.eux, probs=0.8)
  min.f.eux.value.30 <- quantile(f.eux.values$f.eux, probs=0.3)
  max.f.eux.value.70 <- quantile(f.eux.values$f.eux, probs=0.7)
  min.f.eux.value.40 <- quantile(f.eux.values$f.eux, probs=0.4)
  max.f.eux.value.60 <- quantile(f.eux.values$f.eux, probs=0.6)
  
  f.eux.sum <- rbind(f.eux.sum, c(timebin, mean.f.eux.value, 
                                  median.f.eux.value,
                                  min.f.eux.value.05,
                                  max.f.eux.value.95,
                                  min.f.eux.value.10,
                                  max.f.eux.value.90,
                                  min.f.eux.value.20,
                                  max.f.eux.value.80,
                                  min.f.eux.value.30,
                                  max.f.eux.value.70,
                                  min.f.eux.value.40,
                                  max.f.eux.value.60))
}

names(f.eux.sum) <- c("timebin", "mean.f.eux.value", 
                      "median.f.eux.value",
                      "min.f.eux.value.05",
                      "max.f.eux.value.95",
                      "min.f.eux.value.10",
                      "max.f.eux.value.90",
                      "min.f.eux.value.20",
                      "max.f.eux.value.80",
                      "min.f.eux.value.30",
                      "max.f.eux.value.70",
                      "min.f.eux.value.40",
                      "max.f.eux.value.60")

f.eux.sum <- filter(f.eux.sum, !is.na(mean.f.eux.value))

# Fit cross-validated loess models to percentiles plotted (plus mean, median and confidence intervals for interest)
mean.f.eux.fit <- loess.as(f.eux.sum$timebin, f.eux.sum$mean.f.eux,  criterion = c("aicc"))
median.f.eux.fit <- loess.as(f.eux.sum$timebin, f.eux.sum$median.f.eux,  criterion = c("aicc"))
min.f.eux.fit.05 <- loess.as(f.eux.sum$timebin, f.eux.sum$min.f.eux.value.05,  criterion = c("aicc"))
max.f.eux.fit.95 <- loess.as(f.eux.sum$timebin, f.eux.sum$max.f.eux.value.95,  criterion = c("aicc"))
min.f.eux.fit.10 <- loess.as(f.eux.sum$timebin, f.eux.sum$min.f.eux.value.10,  criterion = c("aicc"))
max.f.eux.fit.90 <- loess.as(f.eux.sum$timebin, f.eux.sum$max.f.eux.value.90,  criterion = c("aicc"))
min.f.eux.fit.20 <- loess.as(f.eux.sum$timebin, f.eux.sum$min.f.eux.value.20,  criterion = c("aicc"))
max.f.eux.fit.80 <- loess.as(f.eux.sum$timebin, f.eux.sum$max.f.eux.value.80,  criterion = c("aicc"))
min.f.eux.fit.30 <- loess.as(f.eux.sum$timebin, f.eux.sum$min.f.eux.value.30,  criterion = c("aicc"))
max.f.eux.fit.70 <- loess.as(f.eux.sum$timebin, f.eux.sum$max.f.eux.value.70,  criterion = c("aicc"))
min.f.eux.fit.40 <- loess.as(f.eux.sum$timebin, f.eux.sum$min.f.eux.value.40,  criterion = c("aicc"))
max.f.eux.fit.60 <- loess.as(f.eux.sum$timebin, f.eux.sum$max.f.eux.value.60,  criterion = c("aicc"))

# Generate time dependent f.eux model based on loess models above
mean.f.eux.model <- predict(mean.f.eux.fit, seq(444.7,440.8, -0.02))
median.f.eux.model <- predict(median.f.eux.fit, seq(444.7,440.8, -0.02))
min.f.eux.model.05 <- predict(min.f.eux.fit.05, seq(444.7,440.8, -0.02))
max.f.eux.model.95 <- predict(max.f.eux.fit.95, seq(444.7,440.8, -0.02))
min.f.eux.model.10 <- predict(min.f.eux.fit.10, seq(444.7,440.8, -0.02))
max.f.eux.model.90 <- predict(max.f.eux.fit.90, seq(444.7,440.8, -0.02))
min.f.eux.model.20 <- predict(min.f.eux.fit.20, seq(444.7,440.8, -0.02))
max.f.eux.model.80 <- predict(max.f.eux.fit.80, seq(444.7,440.8, -0.02))
min.f.eux.model.30 <- predict(min.f.eux.fit.30, seq(444.7,440.8, -0.02))
max.f.eux.model.70 <- predict(max.f.eux.fit.70, seq(444.7,440.8, -0.02))
min.f.eux.model.40 <- predict(min.f.eux.fit.40, seq(444.7,440.8, -0.02))
max.f.eux.model.60 <- predict(max.f.eux.fit.60, seq(444.7,440.8, -0.02))

# summary data frame (for viewing only)
f.eux.sum.model <- as.data.frame(cbind(seq(444.7,440.8, -0.02), mean.f.eux.model,
                                       median.f.eux.model,
                                       min.f.eux.model.05,
                                       max.f.eux.model.95,
                                       min.f.eux.model.10,
                                       max.f.eux.model.90,
                                       min.f.eux.model.20,
                                       max.f.eux.model.80,
                                       min.f.eux.model.30,
                                       max.f.eux.model.70,
                                       min.f.eux.model.40,
                                       max.f.eux.model.60))

# summary of percentile bracket minima (for plotting as envelope minima)
f.eux.sum.model.mins <- as.data.frame(cbind(seq(444.7,440.8, -0.02),
                                            min.f.eux.model.05,
                                            min.f.eux.model.10,
                                            min.f.eux.model.20,
                                            min.f.eux.model.30,
                                            min.f.eux.model.40))

# summary of percentile bracket maxima (for plotting as envelope maxima)
f.eux.sum.model.maxs <- as.data.frame(cbind(seq(444.7,440.8, -0.02),
                                            max.f.eux.model.95,
                                            max.f.eux.model.90,
                                            max.f.eux.model.80,
                                            max.f.eux.model.70,
                                            max.f.eux.model.60))                                       

# rename timebins for plotting
names(f.eux.sum.model)[1] <- "timebin"
names(f.eux.sum.model.mins)[1] <- "timebin"
names(f.eux.sum.model.maxs)[1] <- "timebin"

# melt dataframe for plottimng
f.eux.sum.model.melted <- reshape2::melt(f.eux.sum.model, id.vars="timebin")
f.eux.sum.model.mins.melted <- reshape2::melt(f.eux.sum.model.mins, id.vars="timebin")
f.eux.sum.model.maxs.melted <- reshape2::melt(f.eux.sum.model.maxs, id.vars="timebin")

#annotate melted dataframe for generating fill scale when plotted
f.eux.sum.model.envelopes <- cbind(f.eux.sum.model.mins.melted,f.eux.sum.model.maxs.melted, 
                                   as.factor(c(rep("5", nrow(f.eux.sum.model.mins)), rep("10", nrow(f.eux.sum.model.mins)), rep("20", nrow(f.eux.sum.model.mins)), 
                                               rep("30", nrow(f.eux.sum.model.mins)), rep("40", nrow(f.eux.sum.model.mins)))))

names(f.eux.sum.model.envelopes)<- c("timebin", "Min.model", "Min", "timebin2", "Max.model", "Max", "low.bound")

# Relevel for plot fill scale
f.eux.sum.model.envelopes$low.bound <- relevel(f.eux.sum.model.envelopes$low.bound, 5)
f.eux.sum.model.envelopes$low.bound <- relevel(f.eux.sum.model.envelopes$low.bound, 10)
f.eux.sum.model.envelopes$low.bound <- relevel(f.eux.sum.model.envelopes$low.bound, 20)
f.eux.sum.model.envelopes$low.bound <- relevel(f.eux.sum.model.envelopes$low.bound, 30)
f.eux.sum.model.envelopes$low.bound <- relevel(f.eux.sum.model.envelopes$low.bound, 40)

f.eux.sum.model.envelopes <- filter(f.eux.sum.model.envelopes, !(is.na(Min) | is.na(Max)))

# plot time-dependent mass balance confidence envelopes
eux <- ggplot(f.eux.sum.model.envelopes, aes(x=timebin, ymax=Max*100, ymin=Min*100))+
  annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=0.11, ymax=0.3, fill="grey95", color="grey75")+
  annotate(geom="text", y=.15, x=441.52, label=expression("Modern f"[eux]), size=6)+
  geom_ribbon(aes(fill=low.bound),alpha=.8)+
  scale_y_log10(breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100),
                labels = c("0.1", rep("",8),"1.0", rep("",8),"10", rep("",8), "100"))+
  theme_minimal()+scale_x_reverse()+
  ylab(expression("f"[eux]*" (%)"))+
  xlab("Time (Ma)")+
  theme(plot.margin = margin(.1,1,.1,1,"cm"),panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.ticks = element_line(size=1, color="black"), 
        axis.line = element_line(lineend = 'square'), 
        axis.title = element_text(size=34),
        axis.text = element_text(size=26, color="black"),
        legend.title = element_text(size=24),
        legend.text = element_text( size=18),
        axis.ticks.length = unit(5, "points"),
        legend.justification=c(1,1), legend.position=c(.2,.58),
        legend.key = element_rect(color="grey20", size=1.5),
        legend.background = element_rect(fill=alpha('white', 1),color=alpha('white', 1)),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_fill_brewer(palette= "YlGnBu", direction=-1, "Percentile", labels=c(expression("5 to 95"^"th"), expression("10 to 90"^"th"), expression("20 to 80"^"th"),expression("30 to 70"^"th"), expression("40 to 60"^"th"))  )+
  annotate(geom="text", label="A", x=447.25, y=60, size=10)+
  coord_geo(xlim=c(447.5,440.8), ylim=c(.03,130),expand=FALSE, # Geologic timescale added for clarity
            pos = as.list(rep("bottom", 2)),
            abbrv=list(FALSE, FALSE),
            dat = list(adapted_stages_timescale, adapted_period_timescale),
            height = list(unit(2, "lines"), unit(2, "lines")),
            size=list(8,8),
            bord=list(c("left", "bottom", "right"),c("left", "bottom", "right")), lwd=as.list(c(1,1)))
  
  # plot temporal extent of geochemical data
top <- ggplot(f.eux.sum.model.envelopes)+
  annotate(geom="rect", ymin=2.5, ymax=4, xmax=444.7, xmin=443.2, fill="skyblue", colour="grey60")+
  annotate(geom="rect", ymin=0.5, ymax=2, xmax=443.8, xmin=440.8, fill="goldenrod2", colour="grey60")+ 
  annotate(geom="text", label=expression(delta^238*"U"[carb]), x=445.1, y=3.25, size=4.5)+
  annotate(geom="text", label=expression(delta^238*"U"[eux]*" + "*delta^98*"Mo"[eux]), x=444.78, y=1.25, size=4.5)+
  theme_minimal()+scale_x_reverse()+
  coord_cartesian(xlim=c(447.5,440.8),ylim= c(0,4.5), expand = c(0,0))+
  theme(plot.margin = margin(1,1,.1,1,"cm"),panel.border = element_rect(fill=NA,color=NA),
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# rerun if "polygon edge not found" error is returned
fig3a <- ggarrange2(top, eux,  ncol=1, heights=c(.17,1), margin = unit(0.1, "line"))

##==============================================
## Save as pdf
##============================================== 

ggsave("Figure 3a.pdf",fig3a, width = 11, height = 7)


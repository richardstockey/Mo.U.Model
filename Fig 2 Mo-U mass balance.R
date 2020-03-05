library(deSolve)
library(FME)
library(ggplot2)
library(MASS)
library(egg)
library(dplyr)
library(RColorBrewer)
library(fANCOVA)
library(readr)
library(ggridges)

# Set working directory (this will likely be the downloaded folder "Mo.U.Model")

setwd()

################################################
## Fully coupled Mo-U isotope mass balance model
## (see line 184 for start of model function)
################################################

# Model will loop through logarithmically-scaled scenarios of marine euxinia (in log10 sequence from 0.1% to 100%)
log_step_size <- 0.1
power_sequence <- c(seq(-3, 0, log_step_size))
set.seed(1993)

##==============================================
## Define fixed model variables
##==============================================

# Area of modern seafloor
A <- 3.6e+14  # m^2 (from Reinhard et al. 2013)

# Mass of seawater
M <- 1.41e+21 # mass of seawater in kg (from Lau et al 2017 calculation based on Hastings et al. 1996)

# Concentration of trace metals in modern seawater
conc.U <- 14e-9 # mol/kg (Morford and Emerson 1999) - nmol/kg value converted to mol/g
conc.Mo <- 105e-9 # mol/kg (Morford and Emerson 1999)

# Isotopic composition of trace metals in modern seawater
dU.sw0 <- -0.39 # Tissot and Dauphas (2015)
dMo.sw0 <- 2.34 # Kendall et al. 2017 and refs therein

# Mass of trace metals in modern seawater (number of moles)
NU0 <- conc.U * M # mol
NMo0 <- conc.Mo * M # mol

##==============================================
## Define Monte Carlo model variables (max and min)
##==============================================

# Broadly reducing seafloor area 
f.red.min <- 0 # f.red.max must be defined within loop so is hardcoded into parameter table starting line 250

# Riverine fluxes
FU.riv.min <- 2.75e7 #min estimate from Dunk et al. 2002 (42-14.5 Mmol/yr)
FU.riv.max <- 5.65e7 # max estimate from Dunk et al. 2002 (42+14.5 Mmol/yr)
FMo.riv.min <- 1.8e10/95.95 # min of scott, reinhard and chen (Scott et al. 2008, converted to mol/yr)
FMo.riv.max <- 30.0e7 # mol/yr  min of scott, reinhard and chen (Reinhard et al. 2013)

# Riverine isotopic compositions
dU.riv.min <- -0.40 # per mille - 5th percentile of ocean-draining rivers from Andersen et al. 2016 (Chem. Geo.) supplement
dU.riv.max <- -0.10 # per mille - 95th percentile of ocean-draining rivers from Andersen et al. 2016 (Chem. Geo.) supplement
dMo.riv.min <- 0.5 # per mille - lower estimate of Kendall et al. 2017 (pp.709)
dMo.riv.max <- 0.9 # per mille - upper estimate of Kendall et al. 2017

# Accumulation rates of trace metals in euxinic settings - mol per m^2 per year
bU.eux.min <- 5.4e-6 # min measurement from Dunk et al. 2002
bU.eux.max <-  46.2e-6 # max measurement from Dunk et al. 2002 - could also try 16.0 (Saanich Inlet)
bMo.eux.min <-  4800e-6/95.95 # min of scott, reinhard and chen (Reinhard et al. 2013, converted to mol/m2/yr)
bMo.eux.max <- 12000e-6/95.95 # max of scott, reinhard and chen (Scott et al. 2008, converted to mol/m2/yr)

# Accumulation rates of trace metals in broadly reducing settings
bU.red.min <- 0.92e-6 # min measurement from Dunk et al. 2002 (combined min diss + sed from California Margin)
bU.red.max <- 4.37e-6 # max measurement from Dunk et al. 2002 (combined max diss + sed from California Margin)
bMo.red.min <- 2500e-6/95.95 # min of scott, reinhard and chen (Scott et al. 2008, converted to mol/m2/yr)
bMo.red.max <- 2700e-6/95.95 # max of scott, reinhard and chen (Reinhard et al. 2013, converted to mol/m2/yr)

# Accumulation rates of trace metals in oxic settings
bU.ox.min <-  0.82e7/(A*0.8389) # min estimate from Dunk et al. 2002 (Table 5) with reinhard area
bU.ox.max <- 2.04e7/(A*0.8389) # max estimate from Dunk et al. 2002 (Table 5) with reinhard area
bMo.ox.min <- 20e-6/95.95 # min of scott, reinhard and chen (Scott et al. 2008, converted to mol/m2/yr)
bMo.ox.max <- 27.5e-6/95.95 # max of scott, reinhard and chen (Reinhard et al. 2013, converted to mol/m2/yr)

# Fractionation factors in euxinic settings
DU.eux.min <- 0.4 # min tested by Lau et al. 2017
DU.eux.max <- 0.8 # max tested by Lau et al. 2017
DMo.eux.min <- -0.8  # Nägler et al. 2008 max fractionation in euxinic settings
DMo.eux.max <- 0.0 # Quantitative reduction - directly records seawater

# Fractionation factors in broadly reducing settings
DU.red.min <- -0.23 # Negative fractionation is parameterized to mirror positive - see footnote of Table 1
DU.red.max <- 0.23 # Weyer et al. 2008 - max broadly reducing fractionation
DMo.red.min <- -2.8 # max. oxic fractionation used as boundary here 
DMo.red.max <- -0.8 # Nägler et al. 2008 max fractionation in euxinic settings

# Fractionation factors in oxic settings
DU.ox.min <- -0.043 # from Wei et al. 2018
DU.ox.max <- 0.0294118 # from Lau et al. 2017
DMo.ox.min <- -3.0 # from Kendall et al. 2017
DMo.ox.max <- -2.8 # from Dahl et al. 2017

# Local fractionation factors in euxinic settings - matches 'global' fractionation factors
loc.DU.eux.min <- DU.eux.min
loc.DU.eux.max <- DU.eux.max
loc.DMo.eux.min <- DMo.eux.min
loc.DMo.eux.max <- DMo.eux.max

# Local accumulation rates factors in euxinic settings - matches 'global' rates
loc.bU.eux.min <- bU.eux.min
loc.bU.eux.max <- bU.eux.max
loc.bMo.eux.min <- bMo.eux.min
loc.bMo.eux.max <- bMo.eux.max

# Local fractionation factors in carbonates
loc.DU.carb.min <- 0.2
loc.DU.carb.max <- 0.4

# Maximum area of oxic settings
f.ox.lim.min <- 0.8389 # Reinhard 2013
f.ox.lim.max <- 1 # Physical limit

##==============================================
## Define flux model from description in Reinhard et al. 2013
##==============================================

# Digitize the depth to area model of Menard & Smith 1966
depth_m <- c(0,0.2,1,2,3,4,5,6,7,8,9,10,11)*1000
Area_percent <- c(0,7.49,4.42,4.38,8.5,20.94,31.69,21.2,1.23,0.1,0.03,0.01,0.01) # a 0.01 added to last depth bin so that total = 100% (max ocean depth =10,994m)
Cum_Area_percent <- cumsum(Area_percent)

# Fit Loess model to Menard & Smith 1966 to allow us to call at regular depth intervals
fit.rev <- loess.as(Cum_Area_percent, depth_m)

# Flux equation from Middelburg et al. 1996 to local regression through Menard & Smith 1966
flux <- 1.58 - 0.16 * log(predict(fit.rev, c(seq(5,100,0.01), seq(4.99,0.01,-0.01)))) ## Following Reinhard et al. 2013, expand from 5% then only above when necessary

# Calculate cumulative flux per unit area incremement investigated
mean_flux_per_area <- cumsum(flux) / (seq(1,10000,1)) 

# Calculate vector of flux attenuation coefficients for feux scenarios, standardise to f.eux = 0.001, as the effect modeled 
# here is already incorporated into modern environmental measurements in oxygenated scenarios.
mean_flux_feux <- mean_flux_per_area[c(plyr::round_any(10^power_sequence*10000, 1, f=round))]/mean_flux_per_area[c(plyr::round_any(10^power_sequence[1]*10000, 1, f=round))]

# Initiate data frame for results
Mo.U.Sum <- data.frame(f.eux = double(), dMo=double(), dU=double(), dU.carb=double(), f.red=double(), f.ox=double(), f.tot=double())
Mo.U.Sum.Full <- numeric()
# Initiate loop - timed for utility (current runtime ~550 seconds on RGS MacBook Pro - comment out progress printing for speed!)
print(system.time(for(power in power_sequence){
  
  # increase area of euxinic seafloor in 30, log10-scaled steps between 0.1% and 100% euxinic seafloor
  f.eux <- 10^power 
  
  # carbon flux model for euxinic deposition, standardised to 1 at 0.1% euxinic seafloor
  alpha.eux <- mean_flux_feux[(power+3.1)*10]

  # Define initial parameters 
  # Note, these are not used in the current analyses (and can therefore be ignored), as we are performing a global sensitivity analyses by Monte Carlo. 
  # FME still requires a parameter list, however. 
  pars <- c(DU.eux  = runif(1, min=DU.eux.min, max=DU.eux.max),
            DMo.eux = runif(1, min=DMo.eux.min, max=DMo.eux.max),
            loc.DU.eux = runif(1, min=loc.DU.eux.min, max=loc.DU.eux.max),
            loc.DMo.eux = runif(1, min=loc.DMo.eux.min, max=loc.DMo.eux.max),
            DU.red = runif(1, min=DU.red.min, max=DU.red.max),
            DMo.red = runif(1, min=DMo.red.min, max=DMo.red.max),
            DU.ox = runif(1, min=DU.ox.min, max=DU.ox.max),
            DMo.ox = runif(1, min = DMo.ox.min, max = DMo.ox.max),
            f.red = runif(1, min=0, max=1-f.eux),
            FU.riv = runif(1, min = FU.riv.min, max = FU.riv.max),
            FMo.riv = runif(1, min = FMo.riv.min, max = FMo.riv.min),
            bU.eux = runif(1, min = bU.eux.min, max = bU.eux.max),
            bMo.eux = runif(1, min = bMo.eux.min, max = bMo.eux.max),
            bU.red = runif(1, min = bU.red.min, max = bU.red.max),
            bMo.red = runif(1, min = bMo.red.min, max = bMo.red.max),
            bU.ox = runif(1, min = bU.ox.min, max = bU.ox.max),
            bMo.ox = runif(1, min = bMo.ox.min, max = bMo.ox.max),
            loc.bU.eux = runif(1, min = loc.bU.eux.min, max = loc.bU.eux.max),
            loc.bMo.eux = runif(1, min = loc.bMo.eux.min, max = loc.bMo.eux.max),
            dU.riv = runif(1, min = dU.riv.min, max = dU.riv.max),
            dMo.riv = runif(1, min = dMo.riv.min, max = dMo.riv.max), 
            loc.DU.carb = runif(1, min = loc.DU.carb.min, max = loc.DU.carb.max),
            f.ox.lim = runif(1, min=f.ox.lim.min, max=f.ox.lim.max)
            
)

# Define mass balance function (to be used in sensRange analysis below)  
solve.mass.balance <- function(pars){ 
derivs <- function(t, y, pars){
  with(as.list(c(y,pars)), {
  NU <- y[1]
  NMo <- y[2]
  NU.dU <- y[3]
  NMo.dMo <- y[4]
  dU.sw <- NU.dU/NU
  dMo.sw <- NMo.dMo/NMo
  f.ox <- f.ox.lim - f.eux - f.red 

  if(f.ox < 0){f.ox <- 0} # allowing f.ox.lim to vary occassionally produces negative values, in these rare cases f.ox is set to zero
  
  red.and.eux.steps <- plyr::round_any((log10(f.eux+f.red)+3.1)*10, 1, f=round)
  alpha.red.and.eux <-  mean_flux_feux[red.and.eux.steps] # calculate overall mean flux coefficient to broadly reducing and euxinic settings
  
  # calculate mean flux coefficient to broadly reducing settings only
  alpha.red <- abs(alpha.eux * (power+3.1)*10 - alpha.red.and.eux * red.and.eux.steps)/( red.and.eux.steps - (power+3.1)*10 )
  
  if(is.na(alpha.red) == TRUE){ # In rare cases where f.red is sufficiently small for the above calculations to generate NA for alpha.red,
    # assume that alpha.red is best represented by next expansion step along from lower bound of prescribed euxinia (most representative)
    # If f.eux = 100% here, set alpha.red = 1, but doesnt matter what it is set to as there as f.red = 0 in this scenario.
    if(f.eux==1){alpha.red <- 1}else{
    red.and.eux.steps <- ((power+3.1)*10)+1
    alpha.red.and.eux <-  mean_flux_feux[red.and.eux.steps] 
    alpha.red <- abs(alpha.eux * (power+3.1)*10 - alpha.red.and.eux * red.and.eux.steps)/( red.and.eux.steps - (power+3.1)*10 )
    }
  }
  
  # define mass balance equations
  
  FU.ox <- bU.ox * NU/NU0 * A * f.ox
  FU.red <- bU.red * NU/NU0 * A * f.red * alpha.red
  FU.eux <- bU.eux * NU/NU0 * A * f.eux * alpha.eux
  
  FMo.ox <- bMo.ox * NMo/NMo0 * A * f.ox
  FMo.red <- bMo.red * NMo/NMo0 * A * f.red * alpha.red
  FMo.eux <- bMo.eux * NMo/NMo0 * A * f.eux * alpha.eux
  
  d.NU <- FU.riv - FU.ox - FU.red - FU.eux
  d.NMo <- FMo.riv - FMo.ox - FMo.red - FMo.eux
  
  d.NU.dU.sw <- dU.riv * FU.riv - (dU.sw + DU.ox) * FU.ox - (dU.sw + DU.red) * FU.red - (dU.sw + DU.eux) * FU.eux
  d.NMo.dMo.sw <- dMo.riv * FMo.riv - (dMo.sw + DMo.ox) * FMo.ox - (dMo.sw + DMo.red) * FMo.red - (dMo.sw + DMo.eux) * FMo.eux
  
  # see below for equations linking rock- and seawater-values
  return(list(c(U = d.NU, Mo = d.NMo, d.NU.dU.sw , d.NMo.dMo.sw), dU = dU.sw + loc.DU.eux , dMo = dMo.sw + loc.DMo.eux, 
              dU.carb = dU.sw + loc.DU.carb, f.eux = f.eux, f.ox = f.ox, f.red = f.red))
  })
}

# Define initial conditions based on modern
initial.conditions <- c(NU0, NMo0, NU0*dU.sw0, NMo0*dMo.sw0)

# Model is run dynamically for 10Myrs to ensure steady state - can run for 1Myr and get almost identical result, running for 1 Gyr gives fully identical result to 10Myrs
# Solving as steady state using rootSolve function 'steady' gives indistinguishable results, could change if prefered (ode gives more intuitive output)
time.step <- as.numeric(c(0,1e7))

return(ode(y = initial.conditions, t = time.step, func = derivs, parms=pars, d="vode", maxsteps=50000))
}

# Define minimum and maximum values of each parameter to be explored in global sensitivity analysis
# can print as table as desired
parRanges <- data.frame(min=c(DU.eux.min,
                              loc.DU.eux.min,
                              DMo.eux.min,
                              loc.DMo.eux.min, 
                              f.red.min,
                              FU.riv.min, 
                              FMo.riv.min,
                              bU.eux.min, 
                              bMo.eux.min, 
                              bU.red.min, 
                              bMo.red.min, 
                              bU.ox.min, 
                              bMo.ox.min, 
                              loc.bU.eux.min, 
                              loc.bMo.eux.min, 
                              dU.riv.min, 
                              dMo.riv.min,
                              DU.red.min,
                              DMo.red.min,
                              DU.ox.min,
                              DMo.ox.min,
                              loc.DU.carb.min,
                              f.ox.lim.min
                              
),
                        max=c(DU.eux.max,
                              loc.DU.eux.max,
                              DMo.eux.max,
                              loc.DMo.eux.max, 
                              {1-f.eux}, # Set absolute maximum bound on f.red - must be done inside loop therefore not set above
                              FU.riv.max, 
                              FMo.riv.max, 
                              bU.eux.max, 
                              bMo.eux.max, 
                              bU.red.max, 
                              bMo.red.max, 
                              bU.ox.max, 
                              bMo.ox.max,
                              loc.bU.eux.max, 
                              loc.bMo.eux.max, 
                              dU.riv.max,
                              dMo.riv.max, 
                              DU.red.max,
                              DMo.red.max,
                              DU.ox.max,
                              DMo.ox.max, 
                              loc.DU.carb.max,
                              f.ox.lim.max
))

rownames(parRanges) <- c("DU.eux", 
                         "loc.DU.eux", 
                         "DMo.eux", 
                         "loc.DMo.eux", 
                         "f.red", 
                         "FU.riv", 
                         "FMo.riv", 
                         "bU.eux", 
                         "bMo.eux", 
                         "bU.red", 
                         "bMo.red", 
                         "bU.ox", 
                         "bMo.ox", 
                         "loc.bU.eux", 
                         "loc.bMo.eux",
                         "dU.riv", 
                         "dMo.riv",
                         "DU.red",
                         "DMo.red",
                         "DU.ox",
                         "DMo.ox", 
                         "loc.DU.carb", 
                         "f.ox.lim"
)

  # perform global sensitivity analysis for this f.eux scenario
  Mo.U.Step <- sensRange(parms=pars, func=solve.mass.balance, 
                                           parRange = parRanges, dist = "unif", num=1000)
  
  # store key outputs of sensitivity analysis and add to summary data frame with f.eux values
  Step.sum <- cbind(rep(f.eux, length(Mo.U.Step$dMo1e.07)),  Mo.U.Step$dMo1e.07, Mo.U.Step$dU1e.07, Mo.U.Step$dU.carb1e.07, Mo.U.Step$f.red1e.07, Mo.U.Step$f.ox1e.07, Mo.U.Step$f.eux1e.07 + Mo.U.Step$f.red1e.07 + Mo.U.Step$f.ox1e.07)  
  Mo.U.Sum <- rbind(Mo.U.Sum, Step.sum) # produce summary results dataframe
  Mo.U.Sum.Full <- rbind(Mo.U.Sum.Full, Mo.U.Step) # produce full dataframe to investigate all outputs as desired
  
  # print loop progress (progress bar function is misbehaving...) - optional (remove for speed)
  print(paste((((power+3.1)*10)/length(power_sequence)*100),"% progress"))
  print(paste("feux = ", f.eux))

}))

# rename summary columns for plotting 
names(Mo.U.Sum) <- c("f.eux",  "dMo", "dU",  "dU.carb",  "f.red", "f.ox", "f.tot")

Mo.U.Sum <- filter(Mo.U.Sum, !is.na(dMo) & !is.na(dU))

dMo.map <- ggplot(Mo.U.Sum, aes(x=f.eux*100, y=dMo) ) +
  stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette= "YlGnBu", direction=-1, limits=c(0,1)) +
  theme_bw()+
  ylab(expression(delta^98*"Mo"[eux]*" (‰)"))+xlab(expression("f"[eux]*" (%)"))+
  scale_x_log10(breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100),
                labels = c("0.1", rep("",8),"1.0", rep("",8),"10", rep("",8), "100"), sec.axis = sec_axis(~.+0, labels = NULL))+
  scale_y_continuous(sec.axis = sec_axis(~.+0, labels = NULL),breaks=c(0,1,2,3), labels=c("0", "1.0", "2.0", " 3.0"))+
  coord_cartesian(xlim= c(.1,100), ylim=c(min(Mo.U.Sum$dMo),max(Mo.U.Sum$dMo)), expand = c(0,0))+
  theme( panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.line = element_line(size=2, lineend = 'square'), 
        axis.ticks = element_line(size=1, colour="black"), 
        axis.title = element_text(size=34),
        axis.text =element_text( size=26, colour="black"),
        legend.text = element_text( size=16),
        legend.title = element_text( size=16),
        legend.justification=c(1,1), legend.position=c(.98,.98),
        plot.margin = margin(5,5,5,5),
        legend.background = element_rect(fill=alpha('white', 0.6)),
        axis.title.x = element_blank(),
        axis.text.x =  element_blank(),
        axis.ticks.x.top = element_blank(), 
        axis.ticks.y.right = element_blank(),
        axis.ticks.length = unit(5, "points"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  guides(fill=guide_legend(title="Density"))

dU.map <- ggplot(Mo.U.Sum, aes(x=f.eux*100, y=dU) ) +
  stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette= "YlGnBu", direction=-1, limits=c(0,1)) +
  theme_bw()+
  ylab(expression(delta^238*"U"[eux]*" (‰)"))+xlab(expression("f"[eux]*" (%)"))+
  scale_x_log10(breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100),
                labels = c("0.1", rep("",8),"1.0", rep("",8),"10", rep("",8), "100"), sec.axis = sec_axis(~.+0, labels = NULL))+
  scale_y_continuous(sec.axis = sec_axis(~.+0, labels = NULL))+
  coord_cartesian(xlim= c(.1,100), ylim=c(min(Mo.U.Sum$dU),max(Mo.U.Sum$dU)), expand = c(0,0))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.line = element_line(size=2, lineend = 'square'), 
        axis.ticks = element_line(size=1, colour="black"), 
        axis.title = element_text(size=34),
        axis.text =element_text( size=26, colour="black"),
        legend.text = element_text( size=16),
        legend.title = element_text( size=16),
        legend.justification=c(1,1), legend.position=c(.98,.98),
        legend.background = element_rect(fill=alpha('white', 0.6)),
        plot.margin = margin(5,5,5,5),
        axis.ticks.x.top = element_blank(), 
        axis.ticks.y.right = element_blank(),
        axis.ticks.length = unit(5, "points"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
  guides(fill=guide_legend(title="Density"))

##==============================================
## Isotope model plots and Rhuddanian prediction
##==============================================

# Rhuddanian.Data.RData contains two dataframes: Murzuq (euxinic Murzuq Basin files), Bartlett_Data_Rhuddanian (Rhuddanian  Anticosti carbonates)
load("Rhuddanian.Data.RData")

dMo_lower <- min(Murzuq$d98Mo_auth_Al, na.rm=T) # d98Mo detrital correction using [Al]
dMo_upper <- max(Murzuq$d98Mo_auth_Al, na.rm=T) # d98Mo detrital correction using [Al]

dU_lower <- min(Murzuq$dU238_auth_a, na.rm=T) # d238.02891U detrital correction using [Th]
dU_upper <- max(Murzuq$dU238_auth_a, na.rm=T) # d238.02891U detrital correction using [Th]

dU.carb_lower <- min(Bartlett_Data_Rhuddanian$U_iso, na.rm=T)
dU.carb_upper <- max(Bartlett_Data_Rhuddanian$U_iso, na.rm=T)

rhudd.all <- filter(Mo.U.Sum, dU <= dU_upper # Summarise mass balance results compatible with U and Mo euxinic shale values and U carbonate values
                    & dU >= dU_lower
                    & dMo <= dMo_upper
                    & dMo >= dMo_lower
                    & dU.carb <= dU.carb_upper
                    & dU.carb >= dU.carb_lower)

rhudd.dMo <- filter(Mo.U.Sum, dMo <= dMo_upper # Summarise mass balance results compatible with Mo euxinic shale values 
                    & dMo >= dMo_lower)

rhudd.dU <- filter(Mo.U.Sum, dU <= dU_upper # Summarise mass balance results compatible with U euxinic shale values
                   & dU >= dU_lower)

rhudd.dU.carb <- filter(Mo.U.Sum, dU.carb <= dU.carb_upper # Summarise mass balance results compatible with U carbonate values
                    & dU.carb >= dU.carb_lower)

rhudd.iso.sum <- rbind(rhudd.all, rhudd.dMo, rhudd.dU, rhudd.dU.carb) # Combine into one dataframe for plotting

# Generate and assign names to results in rhudd.iso.sum
names <- c(rep("all", nrow(rhudd.all)), rep("dMo only", nrow(rhudd.dMo)), rep("dU only", nrow(rhudd.dU)), rep("dU carb only", nrow(rhudd.dU.carb)))
rhudd.iso.sum <- cbind(rhudd.iso.sum, names) 

# Generate density plots of each proxy and all proxies combined
dens.ridge <- ggplot(rhudd.iso.sum, aes(x=f.eux*100))+
  annotate(geom="rect", ymin=-Inf, ymax=Inf, xmin=0.11, xmax=0.3, fill="grey90", color="grey60")+
  geom_density_ridges( alpha=.7, size=1,aes(y=names, fill=names), scale=.95) +
  theme_bw()+
  ylab("freq")+
  xlab(expression("f"[eux]*" (%)"))+
  scale_x_log10(breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100),
                labels = c("0.1", rep("",8),"1.0", rep("",8),"10", rep("",8), "100"))+coord_cartesian(xlim= c(.1,100), clip="off", expand =c(0,0))+
  theme(panel.border = element_rect(fill=NA,color="black", size=NA,linetype="solid"),
        axis.ticks = element_line(size=1, colour="black"), 
        axis.title = element_text(size=34),
        axis.text = element_text( size=26, colour="black"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text( size=18),
        legend.justification=c(1,1), legend.position='none',
        axis.title.x = element_blank(),
        axis.text.x =  element_blank(),
        axis.line.x = element_line(size=2,colour = "black"),
        plot.margin = margin(20,5,5,10),
        axis.ticks.y =  element_blank(),
        axis.ticks.length = unit(5, "points"),
        legend.background = element_rect(fill=alpha('white', 0.6)),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_fill_manual(values=c(rgb(166,220,181,maxColorValue=255),"goldenrod",  "skyblue", "goldenrod" ))

## Generate data distributions for right side of plots (and white box for next to density distributions)

dMo.right <- ggplot() +
  theme_bw()+
  coord_cartesian(xlim= c(.1,100), ylim=c(min(Mo.U.Sum$dMo),max(Mo.U.Sum$dMo)), expand = c(0,0))+
  theme(panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position='none',
        element_blank(),
        plot.margin = margin(5,5,5,0),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=dMo_lower, ymax=dMo_upper,fill=rgb(166,220,181,maxColorValue=255), color=rgb(166,220,181,maxColorValue=255))

dU.right <- ggplot() +
  theme_bw()+
coord_cartesian(xlim= c(.1,100), ylim=c(min(Mo.U.Sum$dU),max(Mo.U.Sum$dU)), expand = c(0,0))+
  theme(panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position='none',
        element_blank(),
        plot.margin = margin(5,5,5,0),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=dU_lower, ymax=dU_upper, fill=rgb(166,220,181,maxColorValue=255), color=rgb(166,220,181,maxColorValue=255))

dens.right <- ggplot() +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position='none',
        element_blank(),
        plot.margin = margin(20,5,5,0),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())

fig2.dens.ridge <- ggarrange(dens.ridge, dens.right,dens.right, dMo.map, dMo.right,dens.right, dU.map,dU.right,dens.right, ncol=3, nrow=3,heights=c(.6,1,1), widths=c(1,.05,0))

##==============================================
## Save as pdf
##==============================================

pdf("Figure 2.pdf", width = 6.5, height = 13,  encoding = "MacRoman")
fig2.dens.ridge
dev.off()

# Generate statistics to describe distributions for comparison in text

mean(rhudd.all$f.eux)
median(rhudd.all$f.eux)
quantile(rhudd.all$f.eux, 0.05)
quantile(rhudd.all$f.eux, 0.95)

nrow(rhudd.all)
nrow(rhudd.dU)
nrow(rhudd.dU.carb)
nrow(rhudd.dMo)

mean(rhudd.dU.carb$f.eux)
median(rhudd.dU.carb$f.eux)
quantile(rhudd.dU.carb$f.eux, 0.05)
quantile(rhudd.dU.carb$f.eux, 0.95)

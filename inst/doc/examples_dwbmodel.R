## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=TRUE, results='hide'------------------------------------------------
library(DWBmodelUN)
library(raster)
# Verify that the coordinates of the databases match
Coord_comparison(P_sogamoso, PET_sogamoso)
# Load geographic info of GRU and parameters per cell
data(GRU, param)
# Construction of parameter maps from values by GRU
GRU.maps <- buildGRUmaps(GRU, param)
alpha1_v <- GRU.maps$alpha1
alpha2_v <- GRU.maps$alpha2
smax_v <- GRU.maps$smax
d_v <- GRU.maps$d

# Establish the initial modeling conditions
init <- init_state(GRU.maps$smaxR)
g_v <- init$In_ground
s_v <- init$In_storage
rm(init)

# Load general characteristics of modeling
setup_data <- readSetup(Read = TRUE)
Dates <- seq(as.Date( gsub('[^0-9.]','',colnames(P_sogamoso)[3]), format = "%Y.%m.%d"), 
             as.Date(gsub('[^0-9.]','',tail(colnames(P_sogamoso),1)) , format = "%Y.%m.%d"), by = "month")
Start.sim <- which(Dates == setup_data[8,1]); End.sim <- which(Dates == setup_data[10,1])
# Sim.Period: the 1st two columns of the P and PET are the coordinates of the cells
Sim.Period <- c(Start.sim:End.sim)+2  

# Run DWB model
DWB.sogamoso <- DWBCalculator(P_sogamoso[ ,Sim.Period], 
                              PET_sogamoso[ ,Sim.Period],
                              g_v, s_v, alpha1_v, alpha2_v, smax_v, d_v)

## ----echo=TRUE, results='hide'------------------------------------------------
library(DWBmodelUN)
library(raster)
# Load P and PET databases
data(P_sogamoso, PET_sogamoso)

# Verify that the coordinates of the databases match
Coord_comparison(P_sogamoso, PET_sogamoso)

# Load geographic info of GRU and basins where calibration will be performed
data(GRU,basins)
cellBasins <- cellBasins(GRU, basins)

# Establish the initial modeling conditions
GRU.maps <- buildGRUmaps(GRU, param)
init <- init_state(GRU.maps$smaxR)
g_v <- init$In_ground
s_v <- init$In_storage
rm(init)

# Load general characteristics of modeling
setup_data <- readSetup(Read = TRUE)
Dates <- seq(as.Date( gsub('[^0-9.]','',colnames(P_sogamoso)[3]), format = "%Y.%m.%d"), 
             as.Date(gsub('[^0-9.]','',tail(colnames(P_sogamoso),1)) , format = "%Y.%m.%d"), by = "month")

# For this calibration exercise, the last date of simulation is 
# the same as the final date of calibration
Start.sim <- which(Dates == setup_data[8,1])
End.sim <- which(Dates == setup_data[11,1])
# the first two columns of the P and PET are the coordinates of the cells
Sim.Period <- c(Start.sim:End.sim)+2 
Start.cal <- which(Dates == setup_data[9,1])
End.cal <- which(Dates == as.Date("2004-12-01"))
# the first two columns of the P and PET are the coordinates of the cells
Cal.Period <- c(Start.cal:End.cal)+2  

#Load observed runoff
data(EscSogObs)

# Function that runs the DWB model
NSE_Sogamoso_DWB <- function(parameters, P, PET, g_v,s_v, Sim.Period, EscObs, Cal.Period){
  
  parameters <- as.vector(parameters)
  # Transform the parameters to the format that the model needs
  param <- matrix(parameters, nrow = raster::cellStats(GRU,stat="max"))  
  
  # Construction of parameter maps from values by GRU
  GRU.maps <- buildGRUmaps(GRU, param)
  alpha1_v <- GRU.maps$alpha1
  alpha2_v <- GRU.maps$alpha2
  smax_v <- GRU.maps$smax
  d_v <- GRU.maps$d
  DWB.sogamoso <- DWBCalculator(P_sogamoso[ ,Sim.Period], PET_sogamoso[ ,Sim.Period],
                                g_v,s_v, alpha1_v, alpha2_v, smax_v,d_v, calibration = TRUE)
  Esc.Sogamoso <- varBasins(DWB.sogamoso$q_total, cellBasins$cellBasins)
  
  # model evaluation; in case of possible NA results in the simulation, 
  # add a conditional assingment to a very high value
  sim <- Esc.Sogamoso$varAverage[Cal.Period - 2, ]
  obs <- EscSogObs[Cal.Period - 2, ]
  
  if (sum(!is.na(sim)) == prod(dim(sim))){
    numer <- apply((sim - obs)^2, 2, sum, na.rm = TRUE)
    demom <- apply((obs - apply(obs, 2, mean, na.rm = TRUE))^2, 2, sum, na.rm = TRUE)
    nse.cof <- 1 - numer / demom
  } else {
    nse.cof <- NA
  }
  
  Perf <- (-1)*nse.cof
  if(!is.na(mean(Perf))){ 
    Mean.Perf <- mean(Perf)
  } else {Mean.Perf <- 1e100}
  return(Mean.Perf)
}

# coupling with the DDS algorithm
xBounds.df <- data.frame(lower = rep(0, times = 40), upper = rep(c(1, 2000), times = c(30, 10)))
result <- dds(xBounds.df = xBounds.df, numIter=2, OBJFUN=NSE_Sogamoso_DWB,
              P = P_sogamoso, PET = PET_sogamoso, g_v = g_v, s_v = s_v, Sim.Period = Sim.Period, 
              EscObs = EscSogObs, Cal.Period = Cal.Period)

## ----echo=TRUE, fig.align= "center", fig.height= 4, fig.width= 8--------------
library(DWBmodelUN)
library(dygraphs)
data(P_sogamoso)
P.est <- ts(c(t(P_sogamoso[1, -2:-1])), star = c(2012, 1), frequency = 12)
var <- list("Precipitation" = P.est)
 graphDWB(var, tp = 1, main = "Precipitation Lat:7.0 Lon:-72.94")

## ----echo=TRUE, fig.align= "center", fig.height=5, fig.width=7, dpi = 600-----
library(DWBmodelUN)
library(dygraphs)
data(P_sogamoso, simDWB.sogamoso, EscSogObs)
P.est <- ts(c(t(P_sogamoso[1, -2:-1])), star = c(2012, 1), frequency = 12)
runoff.sim <- ts(simDWB.sogamoso[c(131:192) ,1], star = c(2012, 1), frequency = 12)
runoff.obs <- ts(EscSogObs[c(131:192) ,1] , star = c(2012, 1), frequency = 12)
var <- list("Precipitation" = P.est,"Runoff.sim" = runoff.sim, "Runoff.obs" = runoff.obs)
graphDWB(var, tp = 3, main = "DWB results at Sogamoso Basin")

## ----echo=TRUE, fig.align= "center", fig.height=5, fig.width=7, dpi = 600-----
library(DWBmodelUN)
library(dygraphs)
data(P_sogamoso, PET_sogamoso, simDWB.sogamoso)
P <- ts(c(t(P_sogamoso[1, -2:-1])), star = c(2012, 1), frequency = 12)
PET <- ts(c(t(PET_sogamoso[1, -2:-1])), star = c(2012, 1), frequency = 12)
runoff.sim <- ts(simDWB.sogamoso[c(131:192), 1], star = c(2012, 1), frequency = 12)
var <- list("P" = P,"PET" = PET, "Runoff.sim" = runoff.sim)
graphDWB(var, tp = 4, main = "General Comparison Sogamoso Basin")


##  ========================================================================  ##
##  Ronny Vallejos, Clemente Ferrer, Andrea Vásquez, Rodrigo Díaz             ##
##  and Teresita Marín.                                                       ##
##                                                                            ##
##  Copyright (C) 2024                                                        ##
##  ------------------------------------------------------------------------  ##
##  This program is free software; you can redistribute it and/or modify      ##
##  it under the terms of the GNU General Public License as published by      ##
##  the Free Software Foundation; either version 2 of the License, or         ##
##  (at your option) any later version.                                       ##
##                                                                            ##
##  This program is distributed in the hope that it will be useful,           ##
##  but WITHOUT ANY WARRANTY; without even the implied warranty of            ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             ##
##  GNU General Public License for more details.                              ##
##                                                                            ##
##  You should have received a copy of the GNU General Public License         ##
##  along with this program; if not, a copy is available at                   ##
##  http://www.r-project.org/Licenses/                                        ##
##  ========================================================================  ##

library(readxl)
library(spdep)
library(terra)
library(sp)
library(tidyterra)
library(Matrix)
library(RColorBrewer)

## Read data and add covariates for multivariate analysis ##
## Here, IDP = PI, IDC = MDI, PPC = USC, RCP = RPS ##

valp_2018 <- read_xlsx("prom_com_valp_2018_1.xlsx")
valp_2020 <- read_xlsx("prom_com_valp_2020_1.xlsx")
valp_2022 <- read_xlsx("prom_com_valp_2022_1.xlsx")

valp_2018$IDP <- as.numeric(as.character(valp_2018$IDP))
valp_2018$IDC <- as.numeric(as.character(valp_2018$IDC))
valp_2018$PPC <- as.numeric(as.character(valp_2018$PPC))
valp_2018$RCP <- as.numeric(as.character(valp_2018$RCP))

data_valparaiso  <- data.frame(MAT021 = c(valp_2018$MAT021, valp_2020$MAT021, valp_2022$MAT021),
                               FIS100 = c(valp_2018$FIS100, valp_2020$FIS100, valp_2022$FIS100),
                               IWI131 = c(valp_2018$IWI131, valp_2020$IWI131, valp_2022$IWI131),
                               QUI010 = c(valp_2018$QUI010, valp_2020$QUI010, valp_2022$QUI010),
                               PI = c(valp_2018$IDP, valp_2018$IDP, valp_2018$IDP),
                               MDI = c(valp_2018$IDC, valp_2018$IDC, valp_2018$IDC),
                               USC = c(valp_2018$PPC, valp_2018$PPC, valp_2018$PPC),
                               RPS = c(valp_2018$RCP, valp_2018$RCP, valp_2018$RCP))

data_valparaiso$idx <- 1:length(data_valparaiso$MAT021)

santiago_2018 <- read_xlsx("prom_com_rm_2018_1.xlsx")
santiago_2020 <- read_xlsx("prom_com_rm_2020_1.xlsx")
santiago_2022 <- read_xlsx("prom_com_rm_2022_1.xlsx")

santiago_2018$IDP <- as.numeric(as.character(santiago_2018$IDP))
santiago_2018$IDC <- as.numeric(as.character(santiago_2018$IDC))
santiago_2018$PPC <- as.numeric(as.character(santiago_2018$PPC))
santiago_2018$RCP <- as.numeric(as.character(santiago_2018$RCP))

data_santiago  <- data.frame(MAT021 = c(santiago_2018$MAT021, santiago_2020$MAT021, santiago_2022$MAT021),
                               FIS100 = c(santiago_2018$FIS100, santiago_2020$FIS100, santiago_2022$FIS100),
                               IWI131 = c(santiago_2018$IWI131, santiago_2020$IWI131, santiago_2022$IWI131),
                               QUI010 = c(santiago_2018$QUI010, santiago_2020$QUI010, santiago_2022$QUI010),
                               PI = c(santiago_2018$IDP, santiago_2018$IDP, santiago_2018$IDP),
                               MDI = c(santiago_2018$IDC, santiago_2018$IDC, santiago_2018$IDC),
                               USC = c(santiago_2018$PPC, santiago_2018$PPC, santiago_2018$PPC),
                               RPS = c(santiago_2018$RCP, santiago_2018$RCP, santiago_2018$RCP))

data_santiago$idx <- 1:length(data_santiago$MAT021)

## Load maps of Valparaiso Region and Santiago Metropolitan Region ##

rds_santiago <- readRDS("gadm41_CHL_3_pk.rds")
rds_santiago <- st_as_sf(rds_santiago) %>% filter(NAME_1 == "Santiago Metropolitan")
adj_santiago <- poly2nb(rds_santiago)
W_santiago <- as(nb2mat(adj_santiago, style = "B"), "Matrix")

rds_valparaiso <- readRDS("gadm41_CHL_3_pk.rds")
rds_valparaiso <- st_as_sf(rds_valparaiso) %>% filter(NAME_1 == "Valparaíso")
rds_valparaiso <- rds_valparaiso[-c(1), ] #Remove "Isla de Pascua"
rds_valparaiso <- subset(rds_valparaiso, NAME_3 != "Juan Fernández") # Remove "Juan Fernandez"
adj_valparaiso <- poly2nb(rds_valparaiso)
W_valparaiso <- as(nb2mat(adj_valparaiso, style = "B"), "Matrix")

# Number of areas Santiago
n <- nrow(W_santiago)

IDP <- matrix(NA, ncol = 3, nrow = 3 * n)
IDP[1:n, 1] <- santiago_2018$IDP
IDP[n + 1:n, 2] <- santiago_2018$IDP
IDP[2*n + 1:n, 3] <- santiago_2018$IDP

IDC <- matrix(NA, ncol = 3, nrow = 3 * n)
IDC[1:n, 1] <- santiago_2018$IDC
IDC[n + 1:n, 2] <- santiago_2018$IDC
IDC[2*n + 1:n, 3] <- santiago_2018$IDC

PPC <- matrix(NA, ncol = 3, nrow = 3 * n)
PPC[1:n, 1] <- santiago_2018$PPC
PPC[n + 1:n, 2] <- santiago_2018$PPC
PPC[2*n + 1:n, 3] <- santiago_2018$PPC

RCP <- matrix(NA, ncol = 3, nrow = 3 * n)
RCP[1:n, 1] <- santiago_2018$RCP
RCP[n + 1:n, 2] <- santiago_2018$RCP
RCP[2*n + 1:n, 3] <- santiago_2018$RCP

list_santiago <- as.list(data_santiago)

list_santiago$IDP <- IDP
list_santiago$IDC <- IDC
list_santiago$PPC <- PPC
list_santiago$RCP <- RCP


# Number of areas Valparaiso
n <- nrow(W_valparaiso)

IDP <- matrix(NA, ncol = 3, nrow = 3 * n)
IDP[1:n, 1] <- valp_2018$IDP
IDP[n + 1:n, 2] <- valp_2018$IDP
IDP[2*n + 1:n, 3] <- valp_2018$IDP

IDC <- matrix(NA, ncol = 3, nrow = 3 * n)
IDC[1:n, 1] <- valp_2018$IDC
IDC[n + 1:n, 2] <- valp_2018$IDC
IDC[2*n + 1:n, 3] <- valp_2018$IDC

PPC <- matrix(NA, ncol = 3, nrow = 3 * n)
PPC[1:n, 1] <- valp_2018$PPC
PPC[n + 1:n, 2] <- valp_2018$PPC
PPC[2*n + 1:n, 3] <- valp_2018$PPC

RCP <- matrix(NA, ncol = 3, nrow = 3 * n)
RCP[1:n, 1] <- valp_2018$RCP
RCP[n + 1:n, 2] <- valp_2018$RCP
RCP[2*n + 1:n, 3] <- valp_2018$RCP

list_valparaiso <- as.list(data_valparaiso)

list_valparaiso$IDP <- IDP
list_valparaiso$IDC <- IDC
list_valparaiso$PPC <- PPC
list_valparaiso$RCP <- RCP

###############  Models in Santiago Metropolitan Region  ############### 
## Please, in order to reproduce analysis change manually the subject ##

set.seed(1)

library("INLAMSM")
library("INLA")

# Number of years
k <- 3
# Range of autocorrelation parameter
alpha.min <- 0
alpha.max <- 1

A <- kronecker(Diagonal(k, 1), Matrix(1, ncol = nrow(W_santiago), nrow = 1))
e = rep(0, k)

# IMCAR model
model_imcar <- inla.IMCAR.model(k = k, W = W_santiago)
IMCAR <- inla(IWI131 ~ 0 + IDP + IDC + log(PPC)+ RCP + f(idx, model = model_imcar, extraconstr = list(A = as.matrix(A), e = e)),
              family="gaussian", data = list_santiago, control.compute = list(dic = TRUE, waic = TRUE))

summary(IMCAR)

# MCAR model
model_mcar <- inla.MCAR.model(k = k, W = W_santiago, alpha.min = alpha.min, alpha.max = alpha.max)
MCAR <-  inla(IWI131 ~ 0 + IDP + IDC + log(PPC) + RCP + f(idx, model = model_mcar),
              family="gaussian", data = list_santiago, control.compute = list(dic = TRUE, waic = TRUE))

summary(MCAR)

# M-Model
model_m <- inla.Mmodel.model(k = k, W = W_santiago, alpha.min = alpha.min, alpha.max = alpha.max)
MModel <-  inla(IWI131 ~ 0 + IDP + IDC + log(PPC) + RCP + f(idx, model = model_m),
                family= "gaussian", data = list_santiago, control.compute = list(dic = TRUE, waic = TRUE))

summary(MModel)

## Plot predictions ##

rds_santiago$IWI131_2018 <- santiago_2018$IWI131
rds_santiago$IWI131_2020 <- santiago_2020$IWI131
rds_santiago$IWI131_2022 <- santiago_2022$IWI131

rds_santiago$IWI131_2018_IMCAR <- IMCAR$summary.fitted.values[1:52, "mean"]
rds_santiago$IWI131_2020_IMCAR <- IMCAR$summary.fitted.values[52 + 1:52, "mean"]
rds_santiago$IWI131_2022_IMCAR <- IMCAR$summary.fitted.values[104 + 1:52, "mean"]

rds_santiago$IWI131_2018_MCAR <- IMCAR$summary.fitted.values[1:52, "mean"]
rds_santiago$IWI131_2020_MCAR <- IMCAR$summary.fitted.values[52 + 1:52, "mean"]
rds_santiago$IWI131_2022_MCAR <- IMCAR$summary.fitted.values[104 + 1:52, "mean"]

rds_santiago$IWI131_2018_MModel <- MModel$summary.fitted.values[1:52, "mean"]
rds_santiago$IWI131_2020_MModel <- MModel$summary.fitted.values[52 + 1:52, "mean"]
rds_santiago$IWI131_2022_MModel <- MModel$summary.fitted.values[104 + 1:52, "mean"]

bordercolor <- "gray"
bluecols <- brewer.pal(9, 'Blues')  
newcol <- colorRampPalette(bluecols)
ncols <- 1000
bluecols2 <- newcol(ncols)
  
min_max_breaks <- seq(0, 100, length.out = 1000)
  
# Transform to spatial object
rds_santiago <- as_Spatial(rds_santiago)
  
names_attr <- c("IWI131 2018", "IWI131 IMCAR 2018", "IWI131 MCAR 2018", "IWI131 M-Model 2018", 
                "IWI131 2020", "IWI131 IMCAR 2020", "IWI131 MCAR 2020", "IWI131 M-Model 2020",
                "IWI131 2022", "IWI131 IMCAR 2022", "IWI131 MCAR 2022", "IWI131 M-Model 2022")

spplot(rds_santiago, c("IWI131_2018", "IWI131_2018_IMCAR", "IWI131_2018_MCAR",
                        "IWI131_2018_MModel","IWI131_2020", "IWI131_2020_IMCAR",
                        "IWI131_2020_MCAR", "IWI131_2020_MModel", "IWI131_2022",
                        "IWI131_2022_IMCAR", "IWI131_2022_MCAR", "IWI131_2022_MModel"),
        col = bordercolor, lwd = 0.1, col.regions = bluecols2, at = min_max_breaks,
        names.attr = names_attr, sp.layout = list(list("sp.polygons", rds_santiago, first = TRUE, fill = "orange")))

###############  Models in Valparaiso Region  ############### 
## Please, in order to reproduce analysis change manually the subject ##

set.seed(1)

library("INLAMSM")
library("INLA")

# Number of years
k <- 3
# Range of autocorrelation parameter
alpha.min <- 0
alpha.max <- 1

A <- kronecker(Diagonal(k, 1), Matrix(1, ncol = nrow(W_valparaiso), nrow = 1))
e = rep(0, k)

# IMCAR model
model_imcar <- inla.IMCAR.model(k = k, W = W_valparaiso)
IMCAR <- inla(IWI131 ~ 0 + IDP + IDC + log(PPC)+ RCP + f(idx, model = model_imcar, extraconstr = list(A = as.matrix(A), e = e)),
              family="gaussian", data = list_valparaiso, control.compute = list(dic = TRUE, waic = TRUE))

summary(IMCAR)

# MCAR model
model_mcar <- inla.MCAR.model(k = k, W = W_valparaiso, alpha.min = alpha.min, alpha.max = alpha.max)
MCAR <-  inla(IWI131 ~ 0 + IDP + IDC + log(PPC) + RCP + f(idx, model = model_mcar),
              family="gaussian", data = list_valparaiso,  verbose=TRUE, control.compute = list(dic = TRUE, waic = TRUE))

summary(MCAR)

# M-Model
model_m <- inla.Mmodel.model(k = k, W = W_valparaiso, alpha.min = alpha.min, alpha.max = alpha.max)
MModel <-  inla(IWI131 ~ 0 + IDP + IDC + log(PPC) + RCP + f(idx, model = model_m),
                family= "gaussian", data = list_valparaiso, control.compute = list(dic = TRUE, waic = TRUE))

summary(MModel)

## Plot predictions ##

rds_valparaiso$IWI131_2018 <- valp_2018$IWI131
rds_valparaiso$IWI131_2020 <- valp_2020$IWI131
rds_valparaiso$IWI131_2022 <- valp_2022$IWI131
  
rds_valparaiso$IWI131_2018_IMCAR <- IMCAR$summary.fitted.values[1:36, "mean"]
rds_valparaiso$IWI131_2020_IMCAR <- IMCAR$summary.fitted.values[36 + 1:36, "mean"]
rds_valparaiso$IWI131_2022_IMCAR <- IMCAR$summary.fitted.values[72 + 1:36, "mean"]

rds_valparaiso$IWI131_2018_MCAR <- IMCAR$summary.fitted.values[1:36, "mean"]
rds_valparaiso$IWI131_2020_MCAR <- IMCAR$summary.fitted.values[36 + 1:36, "mean"]
rds_valparaiso$IWI131_2022_MCAR <- IMCAR$summary.fitted.values[72 + 1:36, "mean"]

rds_valparaiso$IWI131_2018_MModel <- MModel$summary.fitted.values[1:36, "mean"]
rds_valparaiso$IWI131_2020_MModel <- MModel$summary.fitted.values[36 + 1:36, "mean"]
rds_valparaiso$IWI131_2022_MModel <- MModel$summary.fitted.values[72 + 1:36, "mean"]

bordercolor <- "gray"
bluecols <- brewer.pal(9, 'Blues')  
newcol <- colorRampPalette(bluecols)
ncols <- 1000
bluecols2 <- newcol(ncols)#apply the function to get 1000 colours
  
min_max_breaks <- seq(0, 100, length.out = 1000)
  
#Transform to spatial object
rds_valparaiso <- as_Spatial(rds_valparaiso)
  
names_attr <- c("IWI131 2018", "IWI131 IMCAR 2018", "IWI131 MCAR 2018", "IWI131 M-Model 2018", 
                  "IWI131 2020", "IWI131 IMCAR 2020", "IWI131 MCAR 2020", "IWI131 M-Model 2020",
                  "IWI131 2022", "IWI131 IMCAR 2022", "IWI131 MCAR 2022", "IWI131 M-Model 2022")
  
spplot(rds_valparaiso, c("IWI131_2018", "IWI131_2018_IMCAR", "IWI131_2018_MCAR",
                        "IWI131_2018_MModel","IWI131_2020", "IWI131_2020_IMCAR",
                        "IWI131_2020_MCAR", "IWI131_2020_MModel", "IWI131_2022",
                        "IWI131_2022_IMCAR", "IWI131_2022_MCAR", "IWI131_2022_MModel"),
        col = bordercolor, lwd = 0.1, col.regions = bluecols2, at = min_max_breaks,
        names.attr = names_attr, sp.layout = list(list("sp.polygons", rds_valparaiso, first = TRUE, fill = "orange")))
  

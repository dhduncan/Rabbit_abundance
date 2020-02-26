###########################################################################
## 01 - Microclimate model for WNP    ##########
###########################################################################

## Background: 

# Set up version of the microclimate model for and ...

## Required output: Want to feed into model of rabbit abundance, soil moisture/potential plant growth (??), 

# Interested in variation between transects + variation through time. 


#################################################################################
# Set-up                                 ########################################
#################################################################################
install.packages('microclima')
devtools::install_github('mrke/NicheMapR')

library(NicheMapR)
library(lubridate)

#################################################################################
# Parameters                             ########################################
#################################################################################

# Simulation settings
loc<-c(143.463718, -37.431736)
ystart<-2016
yfinish<-2017

# Env parameters for simulation  
Usrhyt<-0.1
minshade<-0
maxshade<-30
cap=1 # add organic layer of soil at top? See 2014 microclim paper for discussion. 

# soil mositure 
runmoist<-1 # run soil moisture model? If have yes best to add burn in time before period of interest. 
LAI = 0.1 #leaf area index, used to partition traspiration/evaporation from PET
#microclima.LAI = 0 # leaf area index, used by package microclima for radiation calcs - need to look into these
#microclima.LOR = 1 # leaf orientation for package microclima radiation calcs - need to look into these


#Source of env data
opendap = 0 # opendap = 1, query met grids via opendap - not sure if this works??
dailywind = 0 # 0 = no, 1 = yes, will get error if don't have access to data and try to run with it on
microclima = 1 # Use microclima and elevatr package to adjust solar radiation for terrain? 1=yes, 0=no
soildata = 1 # Extract emissivities from gridded data? 1=yes, 0=no
soilgrids = 1 # query soilgrids.org database for soil hydraulic properties?
terrain = 1 # Use 250m resolution terrain data? 1=yes, 0=no


#################################################################################
# Run model                              ########################################
#################################################################################
setwd("~/Students/Dave_WP")

micro <- micro_aust(loc = loc, ystart=ystart, yfinish=yfinish, Usrhyt = Usrhyt, minshade=minshade, maxshade=maxshade,cap=cap, LAI=LAI,
                    runmoist=runmoist,opendap = opendap, dailywind = dailywind, microclima = microclima, soildata=soildata, soilgrids = soilgrids, 
                    terrain=terrain, write_input=1) 

metout<-as.data.frame(micro$metout)
shadmet<-as.data.frame(micro$shadmet)
soil<-as.data.frame(micro$soil)
shadsoil<-as.data.frame(micro$shadsoil)  
soilmoist<-as.data.frame(micro$soilmoist) 

## 
metout$dates<-seq(dmy_hm(paste0("01/01/",ystart,"00:00")), dmy_hm(paste0("31/12/",yfinish,"23:00")), by="hour") 
shadmet$dates<-soil$dates<-shadsoil$dates<-metout$dates



#################################################################################
# Old code to run model - query databases  ######################################
#################################################################################

ystart <- 2016 # start year
yfinish <- 2017 # end year
nyears<-yfinish-ystart+1 # integer, number of years for which to run the microclimate model

stD<-paste0(ystart,"-01-01 00:00:00") # start date
endD<-paste0(yfinish,"-12-31 23:00:00") # end date
dates<-seq(ymd_hms(stD),ymd_hms(endD),by='hours') 

#### STEP 1: Run microclimate model ####
## Run the micro_aust function for each site.

longlat<-c(142.04, -35.38) # Wyperfield
Usrhyt <- 0.1 # m
minshade <- 0 # min simulated shade #
maxshade <- 30 # max simulated shade #
LAI <- 0.1 # increased from default 0.1, reduce water in soil
REFL <- 0.20 #Not used if soildata=1
soildata <- 1 
soiltype <- 6 # Not used if soildata=1. Use http://www.asris.csiro.au/mapping/viewer.html if you want to confirm soil type 
ERR = 1 
RUF <- 0.004 # Default
SLE <- 0.98 #Not used if soildata=1
Thcond <- 2.5 # soil minerals thermal conductivity (W/mC) defult = 2.5
SpecHeat <- 870 
Density = 2560
DEP <- c(0., 1.,  2.5, 7.5,  10,  20.,  40.,  60.,  100.,  200.) # Soil nodes (cm) - set to match logger depths


## Get depth-specific soil properties 
#set source to wherever you extract depth-specific soil properties from
source("C:/Users/nbriscoe/OneDrive - The University of Melbourne/Documents/Students/Will/Microclimate_modelling/soil_hydro.R")
#get this file?
#option to extract these properties? - Nat to look into
prevdir<-getwd()
setwd('Y:')
cmd<-paste("R --no-save --args ",loc[1]," ",loc[2]," < extract.R",sep='')
system(cmd)
soilpro<-read.csv('data.csv')
setwd(prevdir)
soilpro[,1]<-c(2.5,7.5,22.5,45,80,150)
colnames(soilpro)[1] <- 'depth'

#pedotransfer functions
soil.hydro<-soil.hydro(soilpro = as.data.frame(soilpro), DEP = DEP)
PE<-soil.hydro$PE
BB<-soil.hydro$BB
BD<-soil.hydro$BD
KS<-soil.hydro$KS
BulkDensity <- BD[seq(1,19,2)]*1000 #soil bulk density, kg/m3

# search through observed textures and find the nearest match to Campell and Norman's Table 9.1
stypes<-NULL
for(m in 1:nrow(soilpro)){
  ssq<-(CampNormTbl9_1[,2]-soilpro[m,4]/100)^2 + (CampNormTbl9_1[,3]-soilpro[m,3]/100)^2
  stypes[m]<-which(ssq==min(ssq))
}

# produce a table of the qualitative soil profile
soils<-as.character(CampNormTbl9_1[,1])
profile<-as.data.frame(cbind(soilpro[,1],soils[stypes]), stringsAsFactors=FALSE)
profile[,1]<-as.numeric(profile[,1])
colnames(profile)<-c("depth","soiltype")
print(profile)

### run microclimate model
## cap is organic surface layer depth and maxpool is water pooling depth
micro_02<-micro_aust(loc = longlat, LAI = LAI, ystart = ystart, yfinish = yfinish, runmoist=1, maxshade=maxshade,
                    minshade=minshade, evenrain=0,Usrhyt = Usrhyt, DEP = DEP, spatial = "W:/", ERR = 1,dailywind=FALSE,
                    opendap=0, cap=cap, soildata=soildata,Thcond=Thcond,Density=Density, SpecHeat=SpecHeat, maxpool=10, PE=PE,
                    BB=BB, BD=BD, KS=KS, BulkDensity=BulkDensity, RUF=RUF, write_input=1) # Run model + save parameters


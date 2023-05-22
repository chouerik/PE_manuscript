# Max Entropy approach to modelling species distributions through Maxent.
# The lines below was used to infer suitability surfaces served as input for estimating lineage ranges in Choueri et al (in prep).
# to install maxnet package:
#      require(devtools)
#      install_github("mrmaxent/maxnet")

# to install ENMEval version prior to 2.0
#      require(devtools)
#      install_version("ENMeval", version = "0.3.1", repos = "http://cran.us.r-project.org")

# to install other necessary packages
# install.packages(c('raster', 'sp', 'rgeos', 'rgdal', 'dismo', 'fields', 'rJava', 'knitr', 'tidyverse', 'dplyr', 'landscapetools', 'adehabitatHR', 'rgbif', 'rgeos', 'sf', 'CoordinateCleaner', 'rnaturalearth', 'rnaturalearthdata','stringr', 'spThin', 'tidyr', 'Matrix'))

rm(list=ls())
library(maxnet)      
library(raster)
library(sp)
library(rgeos)
library(rgdal)
library(dismo)
library(fields)
library(rJava)
library(ENMeval) # version prior to 2.0
library(knitr)
library(tidyverse)
library (dplyr)
library(landscapetools)
library(adehabitatHR)
library(rgbif)
library(rgeos)
library(sf)
library(CoordinateCleaner)
library(rnaturalearth)
library(rnaturalearthdata)
library(stringr)
library(spThin)
library(tidyr)
library(Matrix)

# This script performs the following steps for each occurrence file in the 'occurrence_data' folder:
# Defines the analysis extent as 5-degree buffers around all occurrence points.
# Filters occurrence points to avoid overlap and redundancy.
# Generates background points within the defined extent, with the quantity determined by the number of spatial cells in the extent.
# Selects the cross-validation method depending on the number of available occurrence points for the species.
# Runs Maxent considering five feature classes, 10 regularization multipliers, and the adopted cross-validation method.
# Selects the best models based on deltaAICc and trained AUC.
# Infers raw and cloglog predictions weighted by the trained AUC of the best models.

#Setting directories
  path <- '' #set '02_species_distribution_models' folder path here!
  setwd(path) #example data
  dir.create(paste0(path,'output_sdms'))

  occ_folder <- paste(path,'sdm_example_data/occurrence_data', sep="")  #folder containing species occurrences points
  var_folder <- paste (path,'sdm_example_data/bioclimatic_variables', sep ="")  #folder containing bioclimatic variables
  amazonia_sl_folder <- paste(path,'sdm_example_data/amz_sensulatissimo', sep="") #folder containing a shapefile of Amazonia sensu latissimo
  output_folder <- paste(path,'output_sdms', sep="") #output folder

#Creating a list of species that will have modeled distributions
  spp_pts = list.files(paste(occ_folder, sep="/"))
  species_list = word(gsub("_"," ",spp_pts), 1,2, sep=" ")

#Loading topographic and bioclimatic variables
  setwd(var_folder)
  var <- dir(pattern = "tif$")
  var_stack <- raster::stack(var)
  projection(var_stack) <- "+proj=longlat + datum=WGS84"
  #plot(var_stack, raster::labels(var_stack)) # plot of variables

#Loading occurrence points
  setwd(occ_folder)
  spp_pts = list.files(occ_folder)
  species_list = word(gsub(".csv","",spp_pts), 1,2, sep=" ")

# Looping to model suitability surfaces for each species
  for (i in 1:length(species_list)){

# STEP 1: selecting target species and extracting occurence points
  target_species = species_list[i]
  print(paste0("Start of species distribution model: ","sp #", i," - ", species_list[i]))
  genus <- word(target_species, 1, sep=" ")
  species <- word(target_species, 2, sep=" ")
  full_presence_data_clean <- read.csv(paste(occ_folder,spp_pts[i], sep = "/")) # occurrence points
  occs.xy <- full_presence_data_clean[c('longitude', 'latitude')]
  sp::coordinates(occs.xy) <- ~ longitude + latitude
  occs.xy_UTM <- occs.xy
  projection(occs.xy_UTM) <- '+proj=utm +zone=20 +datum=WGS84'
  
  # Setting the extent as a buffer around occurrence points
  bgExt <- rgeos::gBuffer(occs.xy_UTM, width = 5) # Considering a buffer of 5° 
  
  # Using the extent to crop the environmental variables
  envsBgCrop <- raster::crop(var_stack, bgExt)
  envsBgMsk <- raster::mask(envsBgCrop, bgExt)
  envsBgMsk_scale <- scale(envsBgMsk) # forcing layers to the same scale
  
  # Ploting occurrence points over cropped environmental layers
  fun <- function() {
  plot(bgExt, add=TRUE, lty=2, lwd=0.5)
  points(occs.xy@coords, bg = "red", pch = 21)
  }
  
  setwd(output_folder)
  pdf(paste0('Variables_cropped_by_buffer_around_', genus,'_', species, '_occurrence_points.pdf'))
  plot(envsBgMsk, main= labels(envsBgMsk), addfun=fun) 
  dev.off()

# STEP 2: Selecting occurrence points 10 km away from neighboring points and converting it to Spatial Points Data Frame
  thinned <- thin( loc.data = full_presence_data_clean,
                lat.col = "latitude", long.col = "longitude",
                spec.col = "species",
                thin.par = 10, reps = 1,
                locs.thinned.list.return = TRUE,
                write.files = FALSE,
                write.log.file = FALSE)

  print(paste0(nrow(full_presence_data_clean)-nrow(as.data.frame(thinned)), ' occurrence points was/were thinned in ', genus, ' ', species, ' dataset' ))
  nr <- nrow(as.data.frame(thinned))
  full_presence_data_clean2 <- cbind.data.frame(as.data.frame(thinned),rep(target_species,nr)) 
  colnames(full_presence_data_clean2)=c("longitude","latitude","species"  )

  bck.na<-envsBgMsk[[1]]
  bck.na[]<-NA
  r<-rasterize(coordinates(full_presence_data_clean2[,c('longitude', 'latitude')]),bck.na,fun='count')
  occ_data_clean_ok.pa<-rasterToPoints(r,fun=function(x){x>0}, spatial=T) 
  
  # Ploting thinned occurrence points over cropped environmental layers
  fun2 <- function() {
    plot(bgExt, add=TRUE, lty=2, lwd=0.5)
    points(occ_data_clean_ok.pa@coords, bg = "blue", pch = 21)
  }
  
  pdf(paste0('Thinned_occurrence_points_', genus,'_', species, '.pdf'))
  plot(envsBgMsk, main= labels(envsBgMsk), addfun=fun2) 
  dev.off()

# STEP 3: Generating background points. Its amount will depend on the size of the extent defined for the species.
  if (ncell(envsBgMsk)>200000){
    bg <- randomPoints(envsBgMsk,20000)
  } else {
    bg <- randomPoints(envsBgMsk,ncell(envsBgMsk)*0.1)
  }

# STEP 4: Defining the partition method for estimating model cross-validation based on the number of occurrence points. 
  method <- ifelse(nrow(occ_data_clean_ok.pa@data)<21, "jackknife" ,"block") #setting the partition method for estimating model cross-validation based on the number of occurrence points.
  mmm <- ifelse(nrow(occ_data_clean_ok.pa@data)<25, 10 ,5) # setting the number of k-folds, if "block method" is selected
  
# STEP 5: Running Maxent 
  mx3 <- ENMevaluate(env=envsBgMsk_scale, occ=coordinates(occ_data_clean_ok.pa), bg.coords=bg,
                 RMvalues=seq(0.5,5,0.5),
                 fc=c("L","LQ","LQP","T","LQT"),
                 method=method, kfolds=mmm, parallel = TRUE, numCores=3) # set the number of cores to be used

  # Saving results in a table and ordering it by best AUCs values.
  ENMresults <- mx3@results[order(mx3@results$delta.AICc),]
  write.table(ENMresults, paste(genus,species,'ENMevaluate.csv', sep="_"), sep="\t", row.names=FALSE)
  
# STEP 6: Selecting models with deltaAICc < 2 and trained AUC > 0.5.
  nn <- ifelse(sum(mx3@results$delta.AICc < 2 &  mx3@results$train.AUC>0.5,na.rm = TRUE)==0,10,2) 
  ENMv_dAICc <- subset(ENMresults,delta.AICc <=nn)
  ENMv_dAICc <- subset(ENMv_dAICc,train.AUC > 0.5)

# STEP 7: Calculating raw predictions weighted by trained AUC.
  mx_raw_pred_AIC2 <- mx3@predictions[[which (mx3@results$delta.AICc<=nn  & mx3@results$train.AUC>0.5)]]
  weights_AUC <- mx3@results[which (mx3@results$delta.AICc<=nn & mx3@results$train.AUC>0.5),]$train.AUC

  if (length(weights_AUC)>1){
    raw_pred <- weighted.mean(x=mx_raw_pred_AIC2, w=weights_AUC)
  } else {
    raw_pred <- mx_raw_pred_AIC2
  }

  projection(raw_pred) <- "+proj=longlat + datum=WGS84"
  
  #Calculating cloglog prediction
  mx_pred_AIC2 <- mx3@models[which (mx3@results$delta.AICc<=nn  & mx3@results$train.AUC>0.5)]
  weights_AUC <- mx3@results[which (mx3@results$delta.AICc<=nn & mx3@results$train.AUC>0.5),]$train.AUC

  norm_fun <- function(mx_pred_AIC2) {
  p <- maxnet.predictRaster(mx_pred_AIC2, envsBgMsk_scale, type="cloglog",clamp = TRUE)
    return(p)
  }
  present_norm_list <- lapply(mx_pred_AIC2,norm_fun)
  present_norm <- stack(present_norm_list)

  #Calculating the weighted average
  cloglog_pred <- weighted.mean(x=present_norm, w=weights_AUC)
  projection(cloglog_pred) <- "+proj=longlat + datum=WGS84"
  
  #Saving raw prediction raster
  outputname <- paste(paste(genus, species, paste(ENMv_dAICc$settings,collapse='&'), "raw", sep="_"), '.tif', sep="")
  writeRaster(raw_pred,outputname, bylayer = TRUE, options = c("COMPRESS=DEFLATE"), format="GTiff", overwrite=TRUE)
  
  #Saving Cloglog prediction raster
  outputname <- paste(paste(genus, species, paste(ENMv_dAICc$settings,collapse='&'), "cloglog", sep="_"), '.tif', sep="")
  writeRaster(cloglog_pred,outputname, bylayer = TRUE, options = c("COMPRESS=DEFLATE"), format="GTiff", overwrite=TRUE)

  #Preparing PDF files
  setwd(amazonia_sl_folder)
  world <- rnaturalearth::ne_countries(scale = "medium", type = "map_units", returnclass = 'sf')
  amz <-shapefile('amazonia_sensulatissimo.shp')
  amz_sl <- gUnaryUnion(amz)
  proj4string(amz_sl) = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  sf_amz_sl<- st_as_sf(amz_sl)
  tib_amz_sl <- as_tibble (sf_amz_sl)
  sf_bgExt_sl <- st_as_sf(bgExt)
  tib_bgExt_sl <- as_tibble(sf_bgExt_sl)

  #Saving raw prediction PDF
  setwd(output_folder)
  raw_pred_df <- as.data.frame(raw_pred, xy = TRUE )%>%drop_na()

  raw_pred_plot <- ggplot () +
                geom_sf(data = world, fill='white', color='gray50', lwd=0.3) +
                ggtitle (paste(genus, species, 'raw prediction',"\n", paste(ENMv_dAICc$settings,collapse=' ; '),  sep=" "))+
                geom_raster(aes(x=x, y=y, fill=raw_pred_df[,3]), data=raw_pred_df) + scale_fill_gradientn(colours=c("#FFFFFFFF","#00F3FF","#000B73")) +
                geom_sf(data = tib_amz_sl, aes(geometry=geometry), color = 'darkgreen',lwd=0.5, fill= NA) +
                geom_point (data = mx3@occ.pts, 
                            mapping = aes(x = LON, y = LAT), pch=20, lwd=2, col='black') +
                coord_sf ( xlim= c(-81, -36), ylim = c(-28, 11)) +
                theme_bw() +
                theme(plot.title = element_text(hjust = 0.5, size= 15), axis.title.y=element_blank(),axis.title.x=element_blank(), legend.title = element_blank())

  pdf(paste(genus, species, paste(ENMv_dAICc$settings,collapse='&'), 'raw.pdf', sep="_"))
  plot(raw_pred_plot)
  dev.off()

  #Saving cloglog prediction PDF
  cloglog_pred_df <- as.data.frame(cloglog_pred, xy = TRUE )%>%drop_na()
  cloglog_pred_plot <-ggplot () +
                geom_sf(data = world, fill='white', color='gray50', lwd=0.3) +
                ggtitle (paste(genus, species, 'cloglog prediction',"\n", paste(ENMv_dAICc$settings,collapse=' ; '),  sep=" "))+
                geom_raster(aes(x=x, y=y, fill=cloglog_pred_df[,3]), data=cloglog_pred_df) + scale_fill_gradientn(colours=c("#FFFFFFFF","#00F3FF","#000B73")) +
                geom_sf(data = tib_amz_sl, aes(geometry=geometry), color = 'darkgreen',lwd=0.5, fill= NA) +
                geom_point (data = mx3@occ.pts, 
                            mapping = aes(x = LON, y = LAT), pch=20, lwd=2, col='black') +
                coord_sf ( xlim= c(-81, -36), ylim = c(-28, 11)) +
                theme_bw() +
                theme(plot.title = element_text(hjust = 0.5, size= 15), axis.title.y=element_blank(),axis.title.x=element_blank(), legend.title = element_blank())

  pdf(paste(genus, species, paste(ENMv_dAICc$settings,collapse='&'), 'cloglog.pdf', sep="_"))
  plot(cloglog_pred_plot)
  dev.off()
  
  #printing best models
  col_int<- c('settings', 'AICc', 'delta.AICc', 'w.AIC','train.AUC', 'avg.test.AUC')
  ENMv_dAICc[,col_int]
  
  print(paste0("Distribution modeling of ", genus, ' ', species, ' complete!'))
}
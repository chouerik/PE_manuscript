# This script uses maximum sensitivity + sensibility method to convert suitability-based LDMs into binary presence/absence surfaces.
# The LDMs resulting from this process serve as spatial input for the phylogeogrpahic betadiversity analyzes performed in Choueri et al (in prep).

# to install the necessary packages
# install.packages(c('raster', 'sdm', 'rgdal', 'rnaturalearth', 'rgeos', 'stringr', 'dismo'))

rm(list = ls())

library(raster)
library(sdm)
library(rgdal)
library(rnaturalearth)
library(rgeos)
library(stringr)
library(dismo)

# setting directories
path <- '' #set '06_phylobetadiversity' folder path here

suitability_ldm_folder <- paste0(path, 'phylobetadiversity_example_data/suitability_LDMs') #folder containing the original LDMs
occ_points_folder <- paste0(path, 'phylobetadiversity_example_data/occ_points_folder') #folder containing the occurrence points of the samples
output_folder <- paste0(path, 'binary_LDMs')
dir.create(paste0(path, 'binary_LDMs'))
amz_sl_folder <- paste0(path,'phylobetadiversity_example_data/amz_sensulatissimo/') #folder containing amazonia sensulatissimo extent

#loading some files and creating some objects
setwd(occ_points_folder)
occ_table <- read.csv('occ_points_lizards_example.csv')

setwd(suitability_ldm_folder)
ldm_list <- list.files(patt='.asc')
lineage_names <- gsub('.asc', '', ldm_list)
lin_nam_graph <-  gsub('_', ' ', lineage_names)
lin_nam_graph_2 <-gsub('Lineage', '- Lineage', lin_nam_graph)

# loading complemetary information to generate pdf maps
#Amazon basin limits
setwd(amz_sl_folder)
amz <-shapefile('amazonia_sensulatissimo.shp')
amz_sl <- gUnaryUnion(amz)
proj4string(amz_sl) = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

#Geopolitical boundaries
world <- rnaturalearth::ne_countries(scale = "medium", type = "map_units", returnclass = 'sf')

#Rivers
#rivers_folder <- '/home/erik/analises_betadiversidade/LDMs_betadiversidade/00_lineages_occ/complementar_information_to_graphs/rivers_naturalearth/'
#setwd(rivers_folder)
#rivers <- shapefile('ne_10m_rivers_lake_centerlines.shp')
rivers <- ne_download(scale=10, type='rivers_lake_centerlines', category = 'physical')
rivers_amz <- crop(rivers, amz_sl)
proj4string(rivers_amz) = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# creating some empty objects to keep the values that will be used in summary statistics table
AUC <- c()
TSS <- c()
Kappa <- c()
max_se_sp <- c()
sum.stats <- c()
zza <- c()

# Starting looping to convert each of the lineages
for (i in 1:length(ldm_list)){
#i=19 
setwd(suitability_ldm_folder)
occs.xy.z <- occ_table[occ_table$lineage == lineage_names[i],] #getting lineage samples 
sp1 <- occs.xy.z[,c("longitude","latitude")] #kept the coordinates columns
sp1$Occurrence <- 1 # add Occurrence column
coordinates(sp1) <- c('longitude','latitude')
proj4string(sp1) = "+proj=longlat +datum=WGS84 +no_defs"
 
ldm_raster <- raster(ldm_list[i])
proj4string(ldm_raster) = "+proj=longlat +datum=WGS84 +no_defs"

# creating random sampling equal to the used on SDMs
if (ncell(ldm_raster)>200000){
  zza[i] <- 20000       # 10% of 200,000
}  else {
    zza[i] <- ncell(ldm_raster)*0.1    # 10% of ldm ncell
  }

psudo <- sampleRandom(ldm_raster,na.rm=TRUE,zza[i],xy=TRUE,sp=TRUE)    

names(psudo) = c("long","lat","var")
psudo@data$Occurrence <- 0
psudo@data <- psudo@data[,"Occurrence",drop=F] # we only keep the column Occurrence
proj4string(psudo) = "+proj=longlat +datum=WGS84 +no_defs"
species <- rbind(sp1,psudo)

obs <- species$Occurrence

pred <- raster::extract(ldm_raster,species) 
ev <- evaluates(obs,pred)
ev1 <- ev@threshold_based$threshold
ev2 <- ev1[2] #  max(specificity+sensitivity). 

# Getting some statistcs to report
AUC[i] <- ev@statistics$AUC
TSS[i] <- ev@threshold_based$TSS[2] 
Kappa[i] <- ev@threshold_based$Kappa[2]
max_se_sp[i] <- ev2
sum.stats[[i]] <- cbind.data.frame(lineage_names[i],AUC[i],TSS[i],Kappa[i],max_se_sp[i])

# using ifelse to convert predicted probabilities into presence-absence (threshold=ev2)
pr.pa <- raster(ldm_raster) # creating an empty raster
pr.pa[] <- ifelse(ldm_raster[] >= ev2,1,0)
pr.pa[] <- ifelse(ldm_raster[] >= ev2,1,NA)

# Bellow, I will only kepp the predictions that are close to the points (for example, when you want to remove overpredictions when you want to create a species distribution map)
#Detect continous clamps in the rasters
clumped <- clump(pr.pa, directions=8)
sp_buf <- raster::buffer(sp1,width=100000) ## here we use a buffer remove suitable areas that are not continuous with the species records

inter <- raster::extract(clumped, sp_buf,na.rm=TRUE) 
inters <- na.exclude(inter[[1]])

raster1 <- match(clumped,inters)

setwd(output_folder)
pdf(paste(lineage_names[i], 'thresholded_map.pdf', sep="_"))
pr.pa2 <- raster(raster1) # creating an empty raster
pr.pa2[] <- ifelse(raster1[] >= 1,1,0)
raster::plot(pr.pa2, legend=FALSE, col = 'green')
maps::map(rivers_amz, add=TRUE, col = 'lightblue')
points(sp1, pch=21)
maps::map(amz, add=TRUE, fill = FALSE)
maps::map('world',add=TRUE, lwd=0.3)
title(paste(lin_nam_graph_2[i], "Thresholded_map"))
dev.off()

final_species_polygons = rasterToPolygons(pr.pa2,dissolve=TRUE)    

#creating outputs
outputname = paste(paste(lineage_names[i], "_threshold_ras" , '.tif', sep=""))
writeRaster(pr.pa2,outputname, bylayer = TRUE, options = c("COMPRESS=DEFLATE"), format="GTiff", overwrite=TRUE)
writeOGR(final_species_polygons,paste0(lineage_names[i], "_threshold_pol", ".shp"), driver="ESRI Shapefile",layer=final_species_polygons@data$layer, overwrite_layer = TRUE)

print(paste0("Conversion complete for ", lineage_names[i]))
}

# creating summary statistics table
final.data.frame <- do.call('rbind', sum.stats)
colnames(final.data.frame) <- c('lineage','AUC','TSS','Kappa','max. spec. sens.')

setwd(output_folder)
write.csv(final.data.frame,"00_LDMs_threshold_evaluation.csv",row.names = FALSE)
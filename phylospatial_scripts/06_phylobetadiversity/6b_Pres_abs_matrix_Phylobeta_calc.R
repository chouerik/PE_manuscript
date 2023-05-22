# This script produces presence/absence matrices and calculates the phylogeographic total dissimilarity, turnover and nestedness
# and was applied in Choueri et al (in prep). To do so, we run the lines below for five different spatial resolutions: 1, 1.25, 
# 1.5, 1.75, and 2 degrees. These values were configured on line 46.
# The data frames generated here were submitted to phylobetadiversity regionalization analysis.

#to install speciesgeocodeR
#require("devtools")
#devtools::install_github("azizka/speciesgeocodeR")

#to install biomapME
#devtools::install_github("GregGuerin/biomapME")
#

# to install other necessary packages
# install.packages(c('ape','plyr', 'RColorBrewer','betapart','maps','raster', 'sp', 'maptools', 'spdep', 'prabclus', 'geosphere', 'betapart', 'rgdal','viridis'), dependencies=TRUE)

rm(list = ls())

library(raster)
library(sp)
library(maptools)
library(spdep)
library(prabclus)
library(geosphere)
library(speciesgeocodeR)
library(biomapME)
library(rgdal)
library(viridis)
library(maps)
library(betapart)
library(RColorBrewer)
library(plyr)
library(ape)

# setting directories
path <- '' #set '06_phylobetadiversity' folder path here
input_folder <- paste0(path, 'binary_LDMs')
output_folder <- paste0(path, 'phylobeta_components_results')
dir.create(output_folder) # creating output folder
amz_sl_folder   <- paste0(path, 'phylobetadiversity_example_data/amz_sensulatissimo') # extent folder
tree_folder <- paste0(path, 'phylobetadiversity_example_data/lineage_tree')

# 1. Creating a grid for Amazônia with the desired resolution
    ex <- extent(c(-80,-44,-21, 9)) # Defining the extent   
    # Producing an empty raster
    r <- raster(ex, nrow=30, ncol=34, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") #check the xmax, xmin, ymax, and ymin of 'ex' object to calculate the approximate number of columns needed. I considered xmax-xmin to ncol and ymax-ymin to nrow
    res(r) <- 1 # Set the spatial resolution here!
    r[] <- 1:ncell(r)

    # Loading Amazonia sensulatissimo mask
    setwd(amz_sl_folder)
    amazonia_shapefile=readShapePoly("amazonia_sensulatissimo.shp")

    plot(r)
    plot(amazonia_shapefile, add = T)

    # masking the grid for Amazonia extent
    r2<-mask(r,buffer(amazonia_shapefile,0.60))
    plot(r2)
    plot(amazonia_shapefile, add = T)
    
    #Polygonizing the grid
    grid_amazonia <- rasterToPolygons(r2)
    plot(grid_amazonia)
    setwd(output_folder)
    writePolyShape(grid_amazonia, "grid_amazonia_buf06.shp")

# 2. Extracting geographic coordinates of occurrences based on the binary LDMs
    setwd(input_folder)
    # Joining the lineages distribution polygons
    shp_files <- list.files(pattern='.shp')
    lin_list <- list()
        for (i in 1: length(shp_files)){
            lin_list[[i]] <- shapefile(shp_files[i])}

    lin_ranges = do.call(rbind, lin_list) 
    gsub("_threshold_pol.shp", "", shp_files)
    lin_ranges@data$layer = gsub("_threshold_pol.shp", "", shp_files)

    crossing<- raster::intersect(lin_ranges,grid_amazonia) 
    colnames(crossing@data) = c("species","layer")

    crossing_dat<-crossing@data # extracting data from polygons
    crossing_dat<- na.exclude(crossing_dat)
    plot(amazonia_shapefile)

    # Getting coordinates from the grid cells centroids
    grid_coord<- centroid(grid_amazonia)
    grid_coord<- as.data.frame(grid_coord)
    grid_coord$layer<-grid_amazonia@data$layer
    colnames(grid_coord) = c("x","y","layer")
    
    # linking lineages occurences to centroid coordinates and creating a data frame
    crossing_dat$long <- grid_coord$x[match(crossing_dat$layer,grid_coord$layer)]     
    crossing_dat$lat <- grid_coord$y[match(crossing_dat$layer,grid_coord$layer)]

    lin_coords<- as.data.frame(crossing_dat[,c(1,3,4)])
    colnames(lin_coords)= c("species", "decimallongitude", "decimallatitude")

    # converting point occurrences to polygons and preparing presence/absence data
    SpGeoCod_object = SpGeoCod(lin_coords, grid_amazonia, areanames = "layer")
    pres_abs1 = SpGeoCod_object$spec_table
    pres_abs2 <-pres_abs1[,order(as.numeric(colnames(pres_abs1)),decreasing = FALSE ) ]
    pres_abs3 = pres_abs2[,!colnames(pres_abs2) %in% "not_classified"]

    # converting for 1/0.
    pres_abs4 = pres_abs3
    pres_abs4<-1*(pres_abs4>0)
    
    #Defining the minimum amount of lineages in a grid cell as 3. Cells containing less than 3 lineages will not to be considered in the analysis
    pres_abs5 = pres_abs4[,which(colSums(pres_abs4)>3)]

# 3. Calculating phylogeographic beta-diversity and its components
    setwd(tree_folder)
    # loading phylogeographic tree
    tree_file <- read.nexus('Lineage_tree_example.tree')
    dichotomic_tree <- multi2di(tree_file, random=FALSE) # resolving basal multichotomies by inserting zero-length branches
    
    # calculating...
    all_spp.phylobetapair<-phylo.beta.pair(t(pres_abs5), dichotomic_tree, index.family="sorensen") 

    # creating data.frames
    all_spp.phylobetapair.beta.total <- as.data.frame(as.matrix(all_spp.phylobetapair$phylo.beta.sor))        
    all_spp.phylobetapair.beta.nestdness <- as.data.frame(as.matrix(all_spp.phylobetapair$phylo.beta.sne))    
    all_spp.phylobetapair.beta.turnover <- as.data.frame(as.matrix(all_spp.phylobetapair$phylo.beta.sim)) 

    # saving dataframes
    setwd(output_folder)
    write.csv(all_spp.phylobetapair.beta.total, 'phylobetadiv_total_dissimilarity_matrix.csv')
    write.csv(all_spp.phylobetapair.beta.nestdness, 'phylobetadiv_nestedness_component_matrix.csv')
    write.csv(all_spp.phylobetapair.beta.turnover, 'phylobetadiv_turnover_component_matrix.csv')
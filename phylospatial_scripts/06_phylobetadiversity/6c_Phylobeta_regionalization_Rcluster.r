# This script infers regionalizations for the phylobetadiversity components promoted in Choueri et al, in prep.
# It has been executed for each of the spatial resolutions defined in the previous script.

# to install the necessary packages
# install.packages(c('recluster','dplyr', 'phytools','geiger','dendextend','viridis', 'raster', 'stringr', 'rgdal'), dependencies=TRUE)

rm(list = ls())

library(recluster)
library(dplyr)
library(phytools)
library(geiger) 
library(dendextend)
library(viridis)
library(raster)
library(stringr)
library(rgdal)

# setting directories
path <- '' #set '06_phylobetadiversity' folder path here

input_folder <- paste0(path, 'phylobeta_components_results')
output_folder <- paste0(path,'phylobeta_regionalizations_results')
dir.create(output_folder)
amz_sl_folder   <- paste0(path, 'phylobetadiversity_example_data/amz_sensulatissimo')
grid_name <- 'grid_amazonia_buf06.shp' # grid created in the previous script

# setting some objects
  diversity_index <- c('total_dissimilarity','nestedness_component','turnover_component')
  tested_method  <- c('average')

# beggining a loop for estimate regionalizations for each phylobetadiversity index computed  
for (z in 1:length(diversity_index)){

# Loading dissimilarity matrix
  setwd(input_folder)
  dissim_matrix <- read.csv(paste0('phylobetadiv_',diversity_index[z],'_matrix.csv'), row.names = 1, header = TRUE,check.names=FALSE)                 # LINHA ORIGINAL:   pres_abs6_p # polygons # 

#Run R-cluster
  turn_cl<-recluster.region(t(dissim_matrix),tr=50,rettree=TRUE,mincl=2,maxcl=50,method=tested_method)
  tcsol_df <- as.data.frame(turn_cl$solutions)
  
  # Selecting K based on the conditions:
  # 1- if the explained dissimilarity is between 90% and 95%, the cluster with best silhouette were selected;
  # 2- if clusters exhibits explained dissimilarity lower than 0.9, the k with max explained dissimilarity were selected;
  # 3- if clusters exhibits explained dissimilarity higher than 0.95, the k with the explained dissimilarity nearest to 0.95 were selected.

  if (max(tcsol_df[,3]) < 0.9){
  filt_k_exdiss <- tcsol_df %>% dplyr::filter(as.character(tcsol_df[,3]) == as.character(max(tcsol_df[,3])))
  } else if (any(between(tcsol_df[,3], 0.90, 0.95)) == TRUE){ 
    filt_k_exdiss <- tcsol_df[between(tcsol_df[,3], 0.90, 0.95),]
  } else { filt_k_exdiss <- tcsol_df[which(abs(tcsol_df[,3]-0.95)==min(abs(tcsol_df[,3]-0.95))),]}

  filt_k_silhou <- filt_k_exdiss %>% dplyr::filter(as.character(filt_k_exdiss[,2]) == as.character(max(filt_k_exdiss[,2])))
  best_k <- as.numeric(filt_k_silhou[,1])

  #if we have a tie among the silhouette values of the best_k proposals, the lower best_k will be selected
  if (length(best_k) > 1){
  filt_k_silhou <- filt_k_silhou %>% dplyr::filter((filt_k_silhou[,1]) ==  min(filt_k_silhou$k))
  best_k <- as.numeric(filt_k_silhou[,1])} 

  #saving best_k infos
  setwd(output_folder)
  write.csv (filt_k_silhou[filt_k_silhou[,1]==best_k,], paste0('best_k_info_',diversity_index[z],'_',tested_method, '.csv'), row.names = FALSE )

  # Creating explained dissimilarity graphics
  pdf(paste0('ExplainedDiss_', diversity_index[z],'_', tested_method,'.pdf'))
  if (5 < best_k & best_k < 45){
  plot(turn_cl$solutions[as.character(best_k-5):as.character(best_k+5),1],turn_cl$solutions[as.character(best_k-5):as.character(best_k+5),3], xlab="Number of clusters", ylab="Explained dissimilarity", lab=c(8,8,7))
  } else if (best_k >=45){
  plot(turn_cl$solutions[40:49,1],turn_cl$solutions[40:49,3], xlab="Number of clusters", ylab="Explained dissimilarity", lab=c(8,8,7))
  } else {
  plot(turn_cl$solutions[1:10,1],turn_cl$solutions[1:10,3], xlab="Number of clusters", ylab="Explained dissimilarity", lab=c(8,8,7))}

  abline(h=.9,col="blue",lwd=0.5, lty=5)
  abline(h=.95,col="red",lwd=0.5, lty=5)
  points(x=best_k, y=turn_cl$solutions[best_k-1,3], pch=3, col='black')
  dev.off()

# savin r-cluster result
  saveRDS(turn_cl,paste0("turn_cl_",tested_method,".rds"))
  write.csv(turn_cl$solutions, paste0('turn_cl_solutions_',diversity_index[z],'_', tested_method,'.csv'), row.names=FALSE)

# Preparing a dendrogram inferred through Ward method
  tree <- cluster::agnes(dissim_matrix, method = "ward")
  dend <- as.dendrogram(tree)                                
  
  cluster_membership<- as.factor(cutree(tree, k = best_k)) # Getting cluster membership classification
  if((identical (names(cluster_membership), colnames(dissim_matrix))) == 'FALSE') stop("Cluster membership and dissimilarity matrix has different names")
  cluster_membership2=cluster_membership[order.dendrogram(dend)] # Change cluster membership label order according to the 'dend' label order

  labels(dend) <- names(cluster_membership2) # Set the names of dendrogram labels to be equal to cluster_membership2

  # Check if cluster_membership is in the same order as dend labels
  if((identical (names(cluster_membership2), labels(dend))) == 'FALSE') stop("cluster_membership2 is not in the same order than dendrogram labels")
  dend2=branches_attr_by_clusters(dend,cluster_membership2)
  # Check the order again...
  if((identical (names(cluster_membership2), labels(dend2))) == 'FALSE') stop("cluster_membership2 is not in the same order than dendrogram-2 labels")

  tree = ladderize(as.phylo(dend2)) # preparing the tree
  if((identical (tree$tip.label, names(cluster_membership))) == 'FALSE') stop("tree$tip.label and names(cluster_membership) do not match!!!")
  
  tip.label<- names(cluster_membership[!duplicated(cluster_membership)]) # setting tip labels
  clade.label<- paste0("PR ", unique(cluster_membership)) #Phylogeographic Region (PR)

  N = sqrt(as.numeric(table(cluster_membership)))   # aesthetic alteration

  sub_tree = ladderize(keep.tip(tree,tip.label))

  tip.label = sub_tree$tip.label # Reset tip labels of the subtree

  # set crown node depth to 1/2 the maximum depth
  depth<-sapply(tip.label,function(x,y) 0.5*y$edge.length[which(y$edge[,2]==which(y$tip.label== x))],y=sub_tree)
  trans<-data.frame(tip.label,clade.label,N,depth)
  rownames(trans)<-NULL
  rm(tip.label,clade.label,N,depth)
  tt<-phylo.toBackbone(sub_tree,trans)

  # Plot colors in the dendrogram
  sub_tree$tip.label = trans$clade.label

  is_tip <- sub_tree$edge[,2] <= length(sub_tree$tip.label)
  ordered_tips <- sub_tree$edge[is_tip, 2]
  order_in_tree = sub_tree$tip.label[ordered_tips]
  
  # Ordering the trans dataframe to add vector of colors
  trans2 = trans
  trans3 = trans2[order(match(trans2$clade.label,order_in_tree)),]

  # Now add the colors of the palette (ex, viridis below)
  trans3$col = turbo(best_k, begin = 1, end = 0.01) # opções de paleta do pacote viridis: viridis, magma, inferno, plasma, cividis, rocket, mako, turbo. Dentro de cada uma é possível definir o começo e o fim do range de cores

  # Back to the original order
  trans4 = trans3[order(match(trans3$clade.label, trans$clade.label)),]

  # Saving the dendrogram
  pdf(paste0('Ward_rcluster_', diversity_index[z],'_', tested_method,'.pdf'))
  par(fig = c(0.3, 0.7, 0, 1)) # change each number per time to see the best adjust
  plot(tt, cex=.8, col=trans4$col, fixed.height = FALSE, print.clade.size = FALSE, lwd = 1)
  dev.off()

# Mapping diversity index
  # Below we'll use the same colors for the grid cells
  setwd(input_folder)
  grid_amazonia <- shapefile (grid_name)

  classified_areas = grid_amazonia[grid_amazonia$layer %in% names(turn_cl$grouping[,best_k-1]),]
  classified_areas$class = turn_cl$grouping[,best_k-1]
  classified_areas$col =  plyr::mapvalues(classified_areas$class, from=unique(turn_cl$grouping[,best_k-1]), to= trans4$col) 

  diversity_index_name1 <- word(diversity_index[z], 1,1,"_")
  diversity_index_name2 <- word(diversity_index[z], 2,2,"_")

  #ploting regionalization results
  setwd(amz_sl_folder)
  amazonia_shapefile=readShapePoly("amazonia_sensulatissimo.shp")
  
  setwd(output_folder)
  pdf(paste0('Mapping_', diversity_index[z],'_',tested_method,'.pdf'))
  plot(amazonia_shapefile, border='transparent', col=rgb(red=0, green=0, blue=0, alpha = 0.07),
       xlim=c(-80,-44),
       ylim=c(-21,9),
       axes=TRUE)
  plot(classified_areas, col = classified_areas$col, border='transparent', add=TRUE)
  plot(grid_amazonia,border="white",add=T)
  maps::map('world',add=TRUE)
  dev.off()

  # Saving the grid
  writeOGR(classified_areas,paste0("grid_rcluster_", diversity_index[z],'_',tested_method,".shp"),driver="ESRI Shapefile",layer =classified_areas@data$class )
}
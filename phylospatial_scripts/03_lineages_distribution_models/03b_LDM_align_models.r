# This script aligns lineages distribution models to the extent that will be used in endemism and diversity analysis.
# The original version of the script was developed by Rosauer et al (2015) and is available at https://github.com/DanRosauer/phylospatial/tree/master/LineageModels/align_models.r
# The lines below were slightly changed for method implementation in Choueri et al., in prep.

rm(list=ls())
library(raster)
library(stringr)

#define directories
base_dir <- '/home/erik/Desktop/Post_doc/manuscritos/endemismo_filogeografico/Scripts_to_publish/03_lineages_distribution_models/'   # set '03_lineages_distribution_models' folder path here

input.dir    = paste0(base_dir, '/output_ldm/misaligned')  # location of existing lineage distribution models
setwd(paste0(base_dir, 'output_ldm'))
dir.create('aligned')
output.dir   = paste0(base_dir, 'output_ldm/aligned')  # location to save aligned lineage distribution models
template_ext = paste0(base_dir, 'ldm_example_data/extent/extent.asc')   # an .asc grid with the extent to which all models will be cropped and aligned.

new_only <- TRUE  # if true, skip grids which are already in the output directory
file.pattern    <- '*.asc$'

input_files = list.files(path=input.dir, pattern=file.pattern, full.names=FALSE, recursive=FALSE, ignore.case=TRUE, include.dirs=FALSE)
output_files= list.files(path=output.dir, pattern=file.pattern, recursive=FALSE, ignore.case=TRUE, include.dirs=FALSE)

template.ras = raster(template_ext)
new_extent =      extent(template.ras)

setwd(input.dir)

for (tfile in input_files) {

  filepath=paste(input.dir,tfile,sep='/')
  outname = paste0(output.dir,'/', word(tfile, 2, sep="LDM_"))

  if ((!tfile %in% output_files) | (! new_only)) {
    grid.ras = raster(tfile)
    grid_ext.ras = extend(grid.ras,new_extent,value=0) # extend to the union of current grid and new extent
    grid_ext.ras = crop(grid_ext.ras,new_extent) # crop back to new extent
    writeRaster(grid_ext.ras,outname,overwrite=TRUE, NAflag=-9999)
    cat("\nExtended asc written for",tfile)

      raster_names <- outname

  } else {
    cat("\nSkipped",tfile)
  }
}

setwd(output.dir)
# now make a raster stack
raster_names<-list.files(pattern='Lineage')
lin.stack <- stack(raster_names)
writeRaster(lin.stack, "lin_models.stack")

maxval  <- stackApply(lin.stack,rep(1,nlayers(lin.stack)),fun=max,filename="max_val.asc")
sumval     <- stackApply(lin.stack,rep(1,nlayers(lin.stack)),fun=sum,filename="sum.asc")
maxprop <- maxval / sumval
writeRaster(maxprop,"max_prop.asc")
maxlin  <- which.max(lin.stack)
writeRaster(maxlin,"max_lin.asc")
stack_names <- data.frame(cbind(1:nlayers(lin.stack),names(lin.stack)))
names(stack_names) <- c("layer_num","layer_name")
write.csv(stack_names,"layer_names")
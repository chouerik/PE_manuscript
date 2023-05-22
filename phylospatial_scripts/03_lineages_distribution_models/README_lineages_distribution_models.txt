This folder contains the scripts used to apply the approach developed by Rosauer et al (2015) in Choueri et al (in prep.). Briefly, the species suitability surface (inferred in Maxent) is partitioned among lineages based on their occurrence points and cost-distances functions. 

For this, we use the cloglog prediction rasters obtained for each species as input together with the corresponding spreadsheet of occurrence of the specific lineages in the script '03_LDM_suitability_partitioning.r'.

To execute the subsequent diversity and endemism analysis, the lineage distribution models needs to be aligned for the extent encompassing the study area. We use the script '03b_LDM_align_models.r' to align all the lineages models to an extent covering the Amazon.


Both the original R scripts were designed by Rosauer et al (2015) and are avaliable at https://github.com/DanRosauer/phylospatial/tree/master/LineageModels/multi-lineage range method.r

The scripts contained here have been subtly altered for the present study. For more information on the scope of our study, methods, results, and references, see Choueri et al, in prep. All scripts were architected on Linux/Mac systems.

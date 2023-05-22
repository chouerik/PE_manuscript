# bGMYC multiple threshold approach applied to intraspecific lineage delimitation.
# The method below was used to propose hypothesis of intraspecific lineages based on evolutionary uniqueness in Choueri et al (in prep).

rm(list=ls())
library(splits) #install.packages("splits", repos="http://R-Forge.R-project.org")
library(bGMYC) # install from the file bGMYC_1.0.2.tar.gz, available in http://nreid.github.io/software/
library(ape) # needs to be a version prior to ape 4. Works fine in ape 3.5. 
             # require(remotes)
             # install_version("ape", version = "3.5", repos = "http://cran.us.r-project.org")


# This script runs bGMYC method to propose intraspecific lineages following the steps below:
# 1. Selecting conspecific terminals from a phylogenetic trees inferred for the genus as a hole;
# 2. Adjustment of priors through a delimitation analyzes using a single tree;
# 3. Computation of delimitation analysis considering multiple trees. Graphics of parameters settings, convergence of runs, and heatmaps of posterior probabilities are also created here;
# 4. Generates a list of proposed lineages for pre-defined multiple thresholds. 
# Obs: The data used here as example is the ultrametric tree composed of samples of the genus Alopoglossus inferred by Choueri et al., in prep).

#setting directories
path <- '' #set '01_lineages_delimitation' folder path here
dir.create(paste0(path,'output_bgmyc'))
setwd(paste0(path, 'lin_del_example_data')) #example data

#STEP 1: Preparing trees that will be input for the intraspecific delimitation analyzes. 

  # The list of terminals listed below correspond to samples classified as Alopoglossus angulatus by Ribeiro-Júnior et al., 2020, 2021), focal species of this example.
  SPlevel_terminals <- c('Alop_BOHE010_ND4',
                         'Alop_PG314_ND4',
                         'Alop_MPEG29468_ND4',
                         'Alop_MPEG27609_ND4',
                         'Alop_MPEG21853_ND4',
                         'Alop_MPEG28239_ND4',
                         'Alop_MPEG27542_ND4',
                         'Alop_MPEG27543_ND4',
                         'Alop_MPEG29762_ND4',
                         'Alop_MPEG28502_ND4',
                         'Alop_MPEG29602_ND4',
                         'Alop_MPEG29603_ND4')
   
   # Single-tree: based on Maximum Clade Credibity tree.
   GENUS_tree <- read.nexus(file="Alopoglossus_BEAST.tree")
   tips_to_drop <- setdiff(GENUS_tree$tip.label, SPlevel_terminals) 
   SPlevel_tree <- drop.tip(GENUS_tree, tips_to_drop)                                    

   # Multiple-trees: sampling of the raw output from BEAST.
   GENUS_trees <- read.nexus(file="Alopoglossus_BEAST.trees") # it is necessary to unzip/extract the file.
   GENUS_trees_samples<-sample(GENUS_trees[1000:10000],100) 
   SPlevel_trees<-lapply(GENUS_trees_samples,drop.tip,tip=tips_to_drop)
   class(SPlevel_trees)<-"multiPhylo"

#STEP 2: Adjusting priors through single-tree analysis. 
   # For details about the arguments of the bgmyc.singlephy function, see http://nreid.github.io/assets/bGMYC_instructions_14.03.12.txt

   results.single <- bgmyc.singlephy(SPlevel_tree, mcmc=50000, burnin=5000, thinning=5,py2=0.75,pc2=1.25, t1=2, t2=100, start=c(0.5,0.5,10),scale=c(5,25,0.5)) # Check if acceptance rates are between .2 and .8. If they are not, regulate values of the "scale" argument according to the instructions of the file above. 
   
   bgmycratessingle <- checkrates(results.single) # Output all parameter values from the run in a table. Distribution of ratios of the Coalescence to Yule rates sampled in the analysis must be greater than 0, with no negative values.
   
   #saving convergence results   
   setwd(paste0(path,'output'))
   pdf('bgmycrates_single.pdf')
   plot(bgmycratessingle)
   dev.off()

# STEP 3. Applying the same settings from single-tree analysis to multi-tree analysis.
   results.multi <- bgmyc.multiphylo(SPlevel_trees, mcmc=50000, burnin=5000, thinning=5,py2=0.75,pc2=1.25, t1=2, t2=100, start=c(0.5,0.5,10),scale=c(5,25,0.5))

   results.spec <- bgmyc.spec(results.multi) # generates a table listing the sampled lineages over the runs and their subsequent probabilities. This step can take a few minutes. 
   write.csv(results.spec, "results.specs.multi.csv") 

   results.probmat <- spec.probmat(results.multi) # calculation of posterior probabilities that samples are from the same lineage. This step can take a few hours. 

   bgmycrates <- checkrates(results.multi) # This step can take a few hours. 
   write.table(bgmycrates, "bgmycrates.multi.csv")

   # saving convergence results
   pdf('bgmycrates_multi.pdf')     
   plot(bgmycrates)
   dev.off()

   # creating heatmap file
   pdf('posterior_probability_heatmap.pdf')
   plot(results.probmat, SPlevel_tree)
   dev.off()

   save.image(file = 'lineages_delimitation.RData') #saving R data

# STEP 4. Creating and exporting a list of proposed lineages for pre-defined multiple thresholds.

   #load('lineages_delimitation.RData')
   #library('bGMYC')

   pts <- seq(0.05,0.95,0.05) #Select the posterior probability thresholds to be considered. 

   # Generating the list of delimitation proposals for each threshold
   for(z in pts) {
   cat(z,file="Alopoglossus_angulatus_bgmyc_points_delimit.txt",append=T, fill=T, labels="point =")
   capture.output(bgmyc.point(results.probmat, z), file="Alopoglossus_angulatus_bgmyc_points_delimit.txt", append=T)
   }
   
   
# For more information on the scope of our study, methods, results, and references, see Choueri et al., in prep.
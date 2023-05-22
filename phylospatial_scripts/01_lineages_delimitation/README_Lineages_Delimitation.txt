This folder constains the scripts used to delimit lineages based on the hypothesis of intraspecific evolutionary uniqueness achieved by three distinct methods (Choueri et al, in prep). Here the scripts used to generate the bGMYC and mPTP results will be presented. We also provide data from Alopoglossus angulatus as an executable example. All the scripts were architected on Linux/Mac systems.

1. BGMYC:
The bGMYC script computed in R environment is available in ‘01_Lineages_delimitation_bGMYC.R’. bGMYC requires R version 3 and ‘ape’ package version 3.5. Details about packages and functions can be found in the R script comments. 
In addition to the conventional bGMYC results, the R script generates a list with the proposed lineages delimited for each threshold of posterior probability.

2. MPTP:

To install mPTP we follow the instructions available at https://github.com/Pas-Kapli/mptp. The implementation of the analysis considered the suggestions proposed by the authors in https://github.com/Pas-Kapli/mptp/wiki/.
We convert the MrBayes genus-level trees in rooted Newick trees to posteriorly crops theyr branches accordingly the delimitation target species. In this way, intraspecific delimitations were proposed for each species independently. 
One of the prerequisites for delimiting using mPTP is the definition of the Minimum Branch Lenght (minbr), which prevents deviations in the analysis considering branch lengths and distances between DNA sequences. For its calculation we submit a fasta file containing the specific sequences and the Newick tree. The following command line was used, exemplified for A. angulatus:


   mkdir output_mptp
  
   mptp --tree_file lin_del_example_data/Aangulatus_mPTP.newick --minbr_auto lin_del_example_data/Aangulatus_mPTP.fasta --output_file output_mptp/Alopoglossus_angulatus_mbrl.out
   
Subsequently, the value of minbr was considered in the execution of the analysis through the following command line:

   mptp --tree_file lin_del_example_data/Aangulatus_mPTP.newick --output_file output_mptp/Alopoglossus_angulatus_mptp --mcmc 1000000 --multi --minbr 0.0024918110 --mcmc_sample 100 --mcmc_runs 4 --mcmc_log


We considered the Average Standard Deviation of Delimitation Support Values ​​(ASDDSV) printed on the screen to verify the congruence between the results of the independent MCMC runs. We consider that the value indicated in "Number of delimited species" as the most probable number of lineages of our approach. In addition, we performed a visual analysis of the tree of the ".combined.svg" file, in which red branches connect terminals exhibiting evolutionary independence.

For more information on the scope of our study, methods, results, and references, see Choueri et al, in prep.

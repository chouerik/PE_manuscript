This folder contains the scripts used for the analyses of phylobetadiversityregionalization applied in Choueri et al., in prep.
As an example, we provide the lineages distribution models and its respective phylogenetic branches obtained for some of the lizard species in the present study.

The script '6a_LDM_converter.R' applies the maximum sensitivity + specificity method to convert suitability-based lineage distribution models to presence/absence surfaces. To do so, it uses the LDMs previously generated together with the occurrence points of its respective lineages. The results generated here are used in the following script to build presence/absence matrices and calculate phylobetadiversity indices.

The script '6b_Pres_abs_matrix_Phylobeta_calc.R' uses the binary LDMs generated in the previous step to build a presence/absence matrix that is used together with the lineage tree to calculate the total dissimilarity, turnover and nestedness of the phylobetadiversity at different spatial resolutions. The results of these indices are applied in the following script for the regionalization of these components.

The script '6c_Phylobeta_regionalization_Rcluster.R' proposes spatial regionalizations of the phylobetadiversity components generated in the previous script, as well as a dendrogram revealing the hierarchical relationship between the proposed regions. The selection of the best number of regionalizations is done automatically through conditions imposed on the silhouette values ​​and explained dissimilarity.


For more information on the scope of our study, methods, results, and references, see Choueri et al, in prep.

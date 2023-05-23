# Spatial correlations between phylogeographic surfaces and their potential determining agents using Spatial 
# Autoregressive Models and Generalized Additive Mixed Models. Choueri et al., in prep.

# This script performs the analyzes following these steps:
# 1. Extracting values from phylogeographic surfaces and their predictors through a standardized resolution 
# grid;
# 2. Tests correlations of models composed of one or multiple predictor variables through SAR
# 3. Tests correlations of models composed of one or multiple predictor variables through GAMM
# 4. Summarizes the results of both methods in a spreadsheet.
# The procedures described between steps 2 and 4 must be performed independently for each phylogeographic 
# surface. To do so, change the values of the 'main_var' and 'select_formula' objects. Details about this 
# procedure are found starting on line 158.

#to install the necessary packages
#install.packages(c('raster', 'rgdal', 'spdep', 'spatialreg', 'mgcv', 'wiqid', 'pbapply', 'robustHD', 'dplyr', 'rcompanion', 'geosphere', 'stringr', 'MuMIn'))

rm(list = ls())

library (raster)
library (rgdal)
library(spdep) 
library(spatialreg)
library(mgcv)
library(wiqid)
library(pbapply)
library(robustHD)
library(dplyr)
library(rcompanion)
library(geosphere)
library(stringr)
library(MuMIn)

# STEP 1: Database preparation:
#   1.a Defining directories
    base.dir <- '' # #set '05_spatial_correlation_analysis' folder path here
    
    diversity_layers.dir <- paste0(base.dir, 'spatial_corr_example_data/lizards_phylogeographic_surfaces/')
    predictor_layers.dir <- paste0(base.dir, 'spatial_corr_example_data/predictor_layers/')              
    output.dir   = paste0(base.dir, 'output_spatial_corr_analysis')              
    setwd(base.dir)
    dir.create(paste0(base.dir, 'output_spatial_corr_analysis'))
    output.sar <- paste0(output.dir, '/out_sar')
    dir.create(paste0(output.dir, '/out_sar'))
    output.gamm <- paste0(output.dir, '/out_gamm')
    dir.create(paste0(output.dir, '/out_gamm'))
    
#   1.b Defining and standardizing a grid to extract the values of the variables by pixel.
    setwd(diversity_layers.dir)
    diversity_files <- list.files (patt='final.asc')
    template_ext = raster(diversity_files[1])  # using PD file to define the extent
    
    res_atual <- raster::res(template_ext)[1]
    res_adjust <- 2 # Original resolution: 0.2
    agg_fact <- res_adjust/res_atual
    
    template_ext.aggr = raster::aggregate(template_ext,agg_fact)    # grid resolution standardization
    proj4string(template_ext.aggr) = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
    template_grid <- rasterToPolygons(template_ext.aggr)            #Standardized grid

#   1.c Loading surfaces
        # Phylogeographic surfaces
         diversity_stack <- stack(diversity_files)
         proj4string(diversity_stack) = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

        # Predictors variables
         setwd(predictor_layers.dir)
         pred_list_files <- list.files(patter='*.asc')
         pred_rasters <- lapply(pred_list_files,raster)
         
         for (i in 1:length(pred_list_files)) {
           proj4string(pred_rasters[[i]]) = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
         }
            
#   1.d Setting resolutions
        # Phylogeographic surfaces
         diversity_stack_new_res = raster::resample(diversity_stack, template_ext.aggr)
         proj4string(diversity_stack_new_res) = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
         
        # Predictor surfaces
         for (i in 1:length(pred_list_files)){
           if (names(pred_rasters[[i]]) == 'hydrobasins_southam_masked_final'){ 
             pred_rasters[[i]] <- raster:: resample(pred_rasters[[i]], template_ext.aggr, method='ngb') # method='ngb' assumes nearest neighbor values, being more suitable for categorical variables
             
           }  else {
             pred_rasters[[i]] <- raster:: resample(pred_rasters[[i]], template_ext.aggr, method='bilinear')
           }
         }
         
         pred_rasters_new_res <- stack(pred_rasters)
         proj4string(pred_rasters_new_res) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

#   1.e Creating data frames based of surfaces values
        for (i in 1:length(diversity_files)){
          diversity_df = raster::extract(diversity_stack_new_res[[i]], template_grid, df=TRUE)
          diversity_col = diversity_df[,2]
          assign(paste(names(diversity_stack_new_res[[i]]),'df',sep='.'),diversity_col)
        }
        
        for (i in 1:length(pred_list_files)){
          preditor_df = raster::extract(pred_rasters_new_res[[i]], template_grid, df=TRUE)
          preditor_col = preditor_df[,2]
          assign(paste(names(pred_rasters_new_res[[i]]),'df',sep='.'),preditor_col)
        }
        
#   1.f Merging data frames
         preditor_list = ls(pattern="final.df")
         dados_ambient = as.data.frame(do.call(cbind,mget(preditor_list)))
         
         is.na.data.frame <- function(x)
         do.call(cbind, lapply(x, is.na))
         dados_ambient[is.na(dados_ambient)] <- 0     
         
#   1.g Standardizing data frame
         # standardizing values
         dados_ambient_std_X <- robustHD::standardize(dados_ambient, centerFun = mean, scaleFun = sd)
         dados_ambient_std_X$hydrobasins_southam_masked_final.df <- dados_ambient$hydrobasins_southam_masked_final.df 
         dados_ambient_std_X$S23stability_masked_final.df <- dados_ambient$S23stability_masked_final.df 
         dados_ambient_std_desordered <- dados_ambient_std_X %>% select(-contains("X"))

         # standardizing column names
         colnames(dados_ambient_std_desordered)[which(names(dados_ambient_std_desordered) == 'lizards_SR_cropped_final.df')]  <- "SR"
         colnames(dados_ambient_std_desordered)[which(names(dados_ambient_std_desordered) == 'lizards_PD_cropped_final.df')]  <- "PD"
         colnames(dados_ambient_std_desordered)[which(names(dados_ambient_std_desordered) == 'lizards_WE_cropped_final.df')]  <- "WE"
         colnames(dados_ambient_std_desordered)[which(names(dados_ambient_std_desordered) == 'lizards_PE_cropped_final.df')]  <- "PE"
         colnames(dados_ambient_std_desordered)[which(names(dados_ambient_std_desordered) == 'lizards_PE_WE_resid_cropped_final.df')]  <- "PE_WE"
         colnames(dados_ambient_std_desordered)[which(names(dados_ambient_std_desordered) == 'lizards_PD_SR_resid_cropped_final.df')]  <- "PD_SR"
         colnames(dados_ambient_std_desordered)[which(names(dados_ambient_std_desordered) == 'dhi_fpar8qa_f_masked_final.df')]  <- "fpar"
         colnames(dados_ambient_std_desordered)[which(names(dados_ambient_std_desordered) == 'dhi_gppqa_f_croped_masked_final.df')]  <- "gpp"
         colnames(dados_ambient_std_desordered)[which(names(dados_ambient_std_desordered) == 'KrigeSoil_fernRcropNEW_masked_final.df')]  <- "fern"
         colnames(dados_ambient_std_desordered)[which(names(dados_ambient_std_desordered) == 'hydrobasins_southam_masked_final.df')]  <- "hydBas"
         colnames(dados_ambient_std_desordered)[which(names(dados_ambient_std_desordered) == 'S23stability_masked_final.df')]  <- "climStb"
         
         # standardizing column position
         dados_ambient_std <- dados_ambient_std_desordered[, c(9, 10, 5, 7, 6, 8, 1, 2, 3, 4, 11)] # numbers refer to the original position of the columns.
        
         # converting NA to 0
         is.na.data.frame <- function(x)
         do.call(cbind, lapply(x, is.na))
         dados_ambient_std[is.na(dados_ambient_std)] <- 0    
         
         # adding latitude and longitude of spatial cells' centroids (required for GAMM)
         lat_long <- centroid(template_grid)
         colnames(lat_long) <- c("longitude","latitude")
         dados_ambient_std <- cbind(dados_ambient_std, lat_long)

#   1.h Creating a vector with the combinations of tested predictor variables
         neigh_df=dados_ambient_std
         neigh_df$hydBas <- as.factor(neigh_df$hydBas)   # classifying categorical values as factors.
         
         x <- colnames(neigh_df[,7:11])  # columns corresponding to phylogeographic surfaces and geographic coordinates should be disregarded.
         env_var_comb = stringi::stri_list2matrix(
           do.call(c, lapply(seq_along(x), combn, x = x, simplify = FALSE)),
           byrow = TRUE
         )

                            #############################################  
                            ### Calculations of spatial correlations  ###
                            #############################################
         
# Here the possible correlations will be checked independently for each phylogeographic surface.
# For that, you must change the value of the 'main_var' object for each run, also changing the value of 'select_formula'(line 172).
# The options are: SR, PD, WE, PE, PD_SR, PE_WE.
# Also, it is necessary to define the amount of cells neighboring the focal cell that will be considered in the search for spatial
# correlations. To that, change the k number in of the object 'my_nb_list' below.

         main_var= "SR"
         my_nb_list <- knn2nb(knearneigh(coordinates(na.omit(template_grid)), k = 8))          

         
# STEP 2: Spatial Autoregressive Models:
#   2.a Setting formulas and running SAR models
        select_formula= as.formula(SR~.) 
        y= 1:nrow(env_var_comb)
  
        all_sar_models = function(y){
          COL.lag.eig <- errorsarlm(select_formula, data=neigh_df[,c(main_var,na.exclude( env_var_comb[y,]))], nb2listw(na.omit(my_nb_list,zero.policy=TRUE , style="W")), method="eigen", quiet=FALSE)
          return(COL.lag.eig)
        }
        
        sar_summaries <- pblapply(y,all_sar_models,cl=3) # set the number of clusters here.
        
        setwd(output.sar)        
        saveRDS(sar_summaries, paste0(main_var,'_SAR_summaries.rds'), ascii = TRUE)

#   2.b Getting AICs from models
        model_list = sar_summaries
        get_AICc = function(model_list){ 
          AICCCC <- wiqid::AICc(model_list)
          return(AICCCC)
        }

        AICc_SAR_env = pblapply(model_list,get_AICc,cl=3) # set the number of clusters here.
        
        AICc_SAR_env2= cbind(as.data.frame(env_var_comb),round(unlist(AICc_SAR_env),2)) # join AIC values with the list of variables
        AICc_SAR_env_top = AICc_SAR_env2[order(AICc_SAR_env2$`round(unlist(AICc_SAR_env), 2)`,decreasing = FALSE),] # Order dataframe to select the best model 
        AICc_SAR_env_top$var <- trimws(gsub("\\+ NA","",paste(AICc_SAR_env_top[,1],AICc_SAR_env_top[,2],AICc_SAR_env_top[,3],AICc_SAR_env_top[,4],AICc_SAR_env_top[,5], sep =" + ")))
   
        #saving AICs file.
        write.csv (AICc_SAR_env_top[,7:6], paste(main_var, 'SAR_correlations_AICc.csv', sep="_"))
        
#   2.c Generating correlograms for the best model
        best_sar_model <- rownames(AICc_SAR_env_top[1,])
        m_nub_sar = sar_summaries[[as.integer(best_sar_model)]]
        sar_cor_best_model = sp.correlogram(my_nb_list,m_nub_sar$residuals,method="I",style="W",order=10,zero.policy=TRUE)
        save(sar_cor_best_model, file = paste(main_var, "SAR_correlogram_MoranI", sep="_"))
        #load(paste(main_var, "SAR_correlogram_MoranI", sep="_")) # to load the object created in the line above

        
# 3. Generalized Additive Mixed Models:
# GAMM applies smoothing functions to continuous variables. As the hydBas variable is categorical, we must 
# create formulas that do not apply the smoothing function when hydBas is composing the tested model.

#   3.a Setting formulas and running GAMM models
    env_var_comb_df <- as.data.frame(env_var_comb)
    env_var_comb_hydBas_NA = env_var_comb
    env_var_comb_hydBas_NA[env_var_comb_hydBas_NA %in% 'hydBas'] <- NA

    formula_list <- list()
    for (i in 1:nrow(env_var_comb)){
      if( "hydBas" %in% env_var_comb_df[i,] )  {
        formula_list[i] = paste0(main_var,"~",paste(paste0("s(", na.exclude( env_var_comb_hydBas_NA[i,]),", k=10",")",collapse="+")), '+hydBas', sep='')}
      else{
        formula_list[i] = paste0(main_var,"~",paste0("s(", na.exclude( env_var_comb[i,]),", k=10",")",collapse="+"))}
    }
    
    formula_list[[3]] <- paste0(main_var, '~s(as.integer(hydBas), k=10)') # When hydBas is the only variable to compose the model, it must be converted to integral values.

    all_gamm_models = function(y){
      gam_formula= as.formula(formula_list[[y]])
      tryCatch({  
        gam_results <-  mgcv::gamm(gam_formula, data=neigh_df, method = 'REML', correlation = corSpatial(form = ~ longitude + latitude, type = 'gaussian'))# linha GAMM
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      return(gam_results)
    }

    gamm_summaries <- pblapply(y,all_gamm_models,cl=3)
    
    setwd(output.gamm)
    saveRDS(gamm_summaries, paste0(main_var, '_GAMM_summaries.rds'), ascii = TRUE)

#   3.b Getting AICs from models
    model_list <- list()
    for (i in 1: length(gamm_summaries)){
      model_list[[i]] = gamm_summaries[[i]]
    }
    
    get_AICc = function(model_list){ 
      AICCC <- MuMIn::AICc(model_list)
      return(AICCC)
    }
    
    AICc_gam_env = pblapply(model_list,get_AICc,cl=2) # set the number of clusters here!
    
    # join AIC values with the list of variables
    AICc_gam_env2= cbind(as.data.frame(env_var_comb),round(unlist(AICc_gam_env),2))
    # Order dataframe to select the best model 
    AICc_gam_env_top = AICc_gam_env2[order(AICc_gam_env2$`round(unlist(AICc_gam_env), 2)`,decreasing = FALSE),]
    AICc_gam_env_top$var <- trimws(gsub("\\+ NA","",paste(AICc_gam_env_top[,1],AICc_gam_env_top[,2],AICc_gam_env_top[,3],AICc_gam_env_top[,4],AICc_gam_env_top[,5], sep =" + ")))
    
    # saving AICs
    write.csv (AICc_gam_env_top[,7:6], paste(main_var, 'GAMM_correlations_AIC.csv', sep="_"))        

#   3.c Generating correlogram for the best model
    best_gamm_model <- rownames(AICc_gam_env_top[1,])
    m_nub_gamm = gamm_summaries[[as.integer(best_gamm_model)]] 
    gam_cor_best_model = sp.correlogram(my_nb_list,as.numeric(resid(m_nub_gamm$lme, type = "normalized")),method="I",style="W",order=10,zero.policy=TRUE) 
    saveRDS(gam_cor_best_model, file = paste(main_var, "GAMM_correlogram_MoranI.rds", sep="_"))
    plot(gam_cor_best_model)
    
    col.W <- nb2listw(my_nb_list, style="W")
    moran.test(as.numeric(resid(m_nub_gamm$lme, type = "normalized")),col.W)

    
# 4. Summarizing spatial correlation results
    sar_gamm_tab_function= function(m_nub_sar,m_nub_gamm,main_var,neigh_df,best_sar_model, best_gamm_model){
      
    #SAR results
    sar_val_names= names(m_nub_sar$coefficients)[-1]
    err_coef = round(m_nub_sar$coefficients[-1],3) # Coeficients
    err_sd_dev = round(m_nub_sar$rest.se[-1],3) # Coef error
    names(err_sd_dev)= NULL # sd from coeficients
    AICc_sar_model = round(wiqid::AICc(m_nub_sar),1) # AIC
    summary_sarlm = summary(m_nub_sar, Nagelkerke=TRUE) #pseudo-R²
    pseudo_sar_r_squared <- summary_sarlm$NK
    pseudo_sar_r_squared <- round(pseudo_sar_r_squared, digits = 3)
    p_value_sar = as.numeric(summary(m_nub_sar)$Coef[-1,4])
    p_value_sar_res= ifelse(p_value_sar >= 0.05, p_value_sar_res <- paste0(""), ifelse(p_value_sar <= 0.05 & p_value_sar >= 0.01, p_value_sar_res <- paste0("*"), ifelse(p_value_sar < 0.009 & p_value_sar >= 0.001, p_value_sar_res <- paste0("**"),  p_value_sar_res <- paste0("***"))))

    #GAMM results
    gamm_val_names= row.names(summary(m_nub_gamm$gam)$s.table)
    f_values = round(summary(m_nub_gamm$gam)$s.table[,3], 3 ) # Coeficients
    names(f_values)= NULL 
    AICc_gamm_model = as.numeric(AICc_gam_env_top[1,6]) # AIC
    gamm_r_squared= round(summary(m_nub_gamm$gam)$r.sq, 3) # R-squared
    p_value_gamm = round(summary(m_nub_gamm$gam)$s.table[,4], 3)  # P-values
    p_value_gamm_res= ifelse(p_value_gamm >= 0.05, p_value_gamm <- paste0(""), ifelse(p_value_gamm <= 0.05 & p_value_gamm >= 0.01, p_value_gamm <- paste0("*"), ifelse(p_value_gamm < 0.009 & p_value_gamm >= 0.001, p_value_gamm <- paste0("**"),  p_value_gamm <- paste0("***"))))  ## p_value_sar refers to the asterisks
      
    #Table organization
    estimate_column_sar = paste0(sar_val_names,", ", err_coef," (",err_sd_dev,")",p_value_sar_res)
    estimate_column_gam= paste0(gamm_val_names,", ", f_values ,p_value_gamm_res)

    col_sar1= c("SARerr",estimate_column_sar[1],pseudo_sar_r_squared,AICc_sar_model)
    col_sar2= cbind(rep("",length(estimate_column_sar)-1),estimate_column_sar[-1],rep("",length(estimate_column_sar)-1),rep("",length(estimate_column_sar)-1))
    col_gam1=c("GAMM",estimate_column_gam[1],gamm_r_squared,AICc_gamm_model)
    col_gam2=cbind(rep("",length(estimate_column_gam)-1),estimate_column_gam[-1],rep("",length(estimate_column_gam)-1),rep("",length(estimate_column_gam)-1))
      
    final_table = rbind(col_sar1,col_sar2,col_gam1,col_gam2)
    colnames(final_table)= c("Model","SAR Coeficients and GAM F-values","R-squared","AICc")
    return(final_table)
    }
    
    #Assembling and saving table
    setwd(output.dir)
    final_table_1dr = sar_gamm_tab_function(m_nub_sar,m_nub_gamm,main_var,neigh_df,best_sar_model, best_gamm_model)
    write.csv(final_table_1dr, paste0(main_var, "_corr_summary_output.csv"))
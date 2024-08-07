frescalo = function(
    occ_data,
    neigh_wts,
    R_star = 0.2703,
    Phi = 0.74,
    missing_data = 2,
    in_parallel = TRUE,
    n_cores = parallel::detectCores(logical = FALSE)-1,
    chunkSize = 250,
    filter_wts = FALSE
){
  # Set object to determine if function created a cluster
    cl_create = FALSE

  # If in_parallel TRUE then check if parallel backend already registered with %dopar%
  if(in_parallel & !foreach::getDoParRegistered()){
    # If not then create a cluster
      cl = parallel::makeCluster(n_cores)
    # Change value of cl_create to determine cluster has been created by this function
      cl_create = TRUE
    # Register with dopar
      doParallel::registerDoParallel(cl)
    # Show message to indicate a cluster has been created and registered
      message("No parallel backend registered with 'dopar'. Creating & registering a cluster of ",n_cores," cores.",sep="")
    # Set on.exit so that cluster will be closed when function exists
      on.exit({
        parallel::stopCluster(cl = cl)
        env <- foreach:::.foreachGlobals
        rm(list=ls(name=env), pos=env)
        rm(env)
      })
  }

  # Determine locations & species in occupancy dataset
    occ_locs = unique(occ_data$location)
    occ_spp = as.character(unique(occ_data$species))  # create list of unique species

  # For each region record presence/absence of each species
    loc_grp = as.factor(rep(c(1:ceiling(length(occ_locs)/chunkSize)),each=chunkSize))
    sSplit = split(s, loc_grp[match(occ_data$location, occ_locs)])  # Split species data up into hectads
    if(in_parallel){
      speciesList <- foreach(i = 1:length(sSplit), .inorder=T, .combine='c', .export = "speciesListFun") %dopar% {
        devtools::load_all()
        speciesListFun(spList = sSplit[[i]], species = occ_spp) # which species are in each location?
      }
    } else {
      speciesList <- foreach(i = 1:length(sSplit), .inorder=T, .combine='c') %do% {
        speciesListFun(spList = sSplit[[i]], species = occ_spp) # which species are in each location?
      }
    }
    # Cleanup
    rm(sSplit, loc_grp)

  # For each focal regional calculate the sampling intensity multiplier alpha_i
    # Determine focal locations from neight_wts
    if(filter_wts){
      # If filter_wts TRUE then subset neigh_wts to keep only those focal locations in occ data
      neigh_wts = neigh_wts[which(neigh_wts$location1 %in% occ_locs),]
    }
    foc_locs = as.character(unique(neigh_wts$location1))
    foc_grp = as.factor(rep(c(1:ceiling(length(foc_locs)/chunkSize)),each=chunkSize))
    dSplit = split(neigh_wts, foc_grp[match(neigh_wts$location1, foc_locs)])

  # Run fres_scaling algorithm across chunks of neighbourhood data
    if(in_parallel){
      out_freq <- foreach(i=1:length(dSplit), .inorder=T, .combine='cfun', .multicombine=TRUE) %dopar% {
        devtools::load_all()
        fresc_scaling(neigh_wts = dSplit[[i]], spp_pa = speciesList, all_loc = occ_locs, all_spp = occ_spp, Phi = Phi, R_star=R_star, missing_data = missing_data)
      }
    } else {
      out_freq <- foreach(i=1:length(dSplit), .inorder=T, .combine='cfun', .multicombine=TRUE) %do% {
        fresc_scaling(neigh_wts = dSplit[[i]], spp_pa = speciesList, all_loc = occ_locs, all_spp = occ_spp, Phi = Phi, R_star=R_star, missing_data = missing_data)
      }
    }

    # Cleanup
    rm(dSplit, foc_grp)

  # Do the Frescalo trend analysis if there are more than 1 year bins (use same location groups as sSplit)
    sSplit2 = split(occ_data, as.factor(occ_data$time))  # Split species data up into year bins
    if(in_parallel){
      trend_out <- foreach(i=1:length(sSplit2), .inorder=T, .combine='cfunTrend', .multicombine=TRUE) %dopar% {
        devtools::load_all()
        fresc_trend(s_data = sSplit2[[i]], out_freq$freq, all_spp = occ_spp)
      }
    } else {
      trend_out <- foreach(i=1:length(sSplit2), .inorder=T, .combine='cfunTrend', .multicombine=TRUE) %do% {
        fresc_trend(s_data = sSplit2[[i]], out_freq$freq, all_spp = occ_spp)
      }
    }

    # Cleanup
    rm(sSplit2)

  # Create output object
    out_obj = c(out_freq,trend_out)
  # Return output
    return(out_obj)
}

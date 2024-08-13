#' Apply frescalo model to occupancy data
#'
#' @inheritParams fresc_scaling
#' @param occ_data occupancy dataset. Supplied as a data.frame with the
#'   following named columns; `time`, `location`, `species`. The dataset should
#'   contain a row corresponding to the recording or presence of the specified
#'   species at the specified location during the specified time period.
#' @param in_parallel Logical determining whether the function is to be run in
#'   parallel. Default `in_parallel = TRUE` will run the frescalo model
#'   functions in parallel where possible.
#' @param n_cores The number of cores/worker processes to spawn when running in
#'   parallel
#' @param chunkSize The size of chunks/batches passed to the workers when
#'   running in parallel.
#' @param filter_wts Logical determining whether the neighbourhoods weights
#'   dataset should be filtered to only keep focal locations that are in the
#'   occupancy dataset or if all locations in `neigh_wts` should be retained.
#'
#' @return A `list` with four elements each of containing a `data.frame`. This
#'   list effectively collates the outputs from `fresc_scaling()` and
#'   `fresc_trend()` functions. The four elements are; `loc`, `freq`,`trend` and
#'   `site_time`. for more details see the documentation for `fresc_scaling()`
#'   and `fresc_trends()`
#'
#' @export
#'
#' @description This function is a wrapper function that processes the data and
#'   calls other functions from the package to fit the two stages of the
#'   frescalo model to the occupancy data supplied. The two stages of the model
#'   are fitted using the functions `fresc_scaling()` and `fresc_trend()` from
#'   the package. The first `fresc_scaling()` determines the rescaling parameter
#'   `alpha` required to adjust the neighbourhood frequency for each
#'   neighbourhood to match the target frequency `Phi`. This alpha value
#'   effectively quantifies the recording effort for that location and is used
#'   to produce adjusted frequencies.
#'
#'   The second `fresc_trend()` uses the outputs from `fresc_scaling()` to
#'   estimate time factors for each species to provide an estimate of how the
#'   adjusted frequencies have changed between time periods.
#'
#' @seealso [fresc_scaling()] for the rescaling part of the frescalo model or
#'   [fresc_trend()] for the estimation of the time factors
#'
#' @examples
#'
#' ## Use test dataset (s) and weights data (d) included with the package
#' out_fres = frescalo(s,d,in_parallel = FALSE,filter_wts = TRUE)
#'
#' ## View outputs
#' # Location
#' head(out_fres[["loc"]])
#'
#' # Frequency
#' head(out_fres[["freq"]])
#'
#' # Trend
#' head(out_fres[["trend"]])
#'
#' # Site Time
#' head(out_fres[["site_time"]])
#'
#'
#'
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

  # If in_parallel TRUE then check if parallel backend already registered with %dopar%
  if(in_parallel & !foreach::getDoParRegistered()){
    # If not then create a cluster
      cl = parallel::makeCluster(n_cores)
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
    sSplit = split(occ_data, loc_grp[match(occ_data$location, occ_locs)])  # Split species data up into hectads
    if(in_parallel){
      speciesList <- foreach(i = 1:length(sSplit), .inorder=T, .combine='c', .export = "speciesListFun") %dopar% {
        #devtools::load_all() # lines need for local testing remove for installed package
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
      freq_out = foreach(i=1:length(dSplit), .inorder=T, .combine='cfun', .multicombine=TRUE) %dopar% {
        #devtools::load_all() # lines need for local testing remove for installed package
        fresc_scaling(neigh_wts = dSplit[[i]], spp_pa = speciesList, all_loc = occ_locs, all_spp = occ_spp, Phi = Phi, R_star=R_star, missing_data = missing_data)
      }
    } else {
      freq_out = foreach(i=1:length(dSplit), .inorder=T, .combine='cfun', .multicombine=TRUE) %do% {
        fresc_scaling(neigh_wts = dSplit[[i]], spp_pa = speciesList, all_loc = occ_locs, all_spp = occ_spp, Phi = Phi, R_star=R_star, missing_data = missing_data)
      }
    }

    # Cleanup
    rm(dSplit, foc_grp)

  # Do the Frescalo trend analysis if there are more than 1 year bins (use same location groups as sSplit)
    sSplit2 = split(occ_data, as.factor(occ_data$time))  # Split species data up into year bins
    if(in_parallel){
      trend_out = foreach(i=1:length(sSplit2), .inorder=T, .combine='cfunTrend', .multicombine=TRUE) %dopar% {
        #devtools::load_all() # lines need for local testing remove for installed package
        fresc_trend(s_data = sSplit2[[i]], freq_out$freq, all_spp = occ_spp)
      }
    } else {
      trend_out = foreach(i=1:length(sSplit2), .inorder=T, .combine='cfunTrend', .multicombine=TRUE) %do% {
        fresc_trend(s_data = sSplit2[[i]], freq_out$freq, all_spp = occ_spp)
      }
    }

    # Cleanup
    rm(sSplit2)

  # Create output object
    out_obj = c(freq_out,trend_out)
  # Return output
    return(out_obj)
}

#' Calculate distances between two sets of sf spatial objects
#'
#' @param sf_foc The focal `sf` spatial object for which you want sets of distances returned
#' @param sf_comp The `sf` object containing the features you want to determine the distances between
#'
#' @return a `matrix` containing the distances
#' @export
#'
#' @examples
#'
dist_sf = function(sf_foc, sf_comp){
  # Convert to centroid
    pt_foc = sf::st_centroid(sf::st_geometry(sf_foc))
    pt_comp = sf::st_centroid(sf::st_geometry(sf_comp))
  # Return distance matrix
    pt_d = sf::st_distance(pt_foc, pt_comp)
  # Convert to numerical object
    return(pt_d)
}

#' Calculate cosine similarity
#'
#' @param df A `data.frame` of variables upon which cosine similarity is to be
#'   calculated for each row.
#' @param foc_inds the indices of `df` which are considered to be the focal rows
#'   for which you want the similarity (with all rows) returned. The default
#'   `foc_inds = NULL` will return the full similarity matrix, othewise the
#'   results will be filtered to rows for the rows specified in `foc_ind`.
#'
#' @return a matrix of similarity scores
#'
#' @examples
cosine_sim = function(df, foc_inds = NULL){
# Modified from post in https://stats.stackexchange.com/questions/31565/compute-a-cosine-dissimilarity-matrix-in-r
  df_mat = as.matrix(df)
  sim = df_mat / sqrt(rowSums(df_mat*df_mat))
  if(!is.null(foc_inds)){
    sim = sim[foc_inds,] %*% t(sim)
  } else {
    sim = sim %*% t(sim)
  }
  return(sim)
}

#' Calculate frescalo weights values for a neighbourhood based on distance and
#' similarity scores
#'
#' @inheritParams calc_neigh_wts
#' @param nr A `data.frame` containing location identifiers, distances and
#'   similarity scores of all locations in the neighbourhood for a focal
#'   location. The columns will be named `location`,`d` and `s` respectively.
#' @param sort Logical specifying whether the `data.frame` that is returned is
#'   to be sorted based on the weights (in descending order).
#'
#' @return a `data.frame` containing two columns; `location2` and `w` containing
#'   the location identifiers for all locations in the neighbourhood and their
#'   corresponding weights
#'
#' @examples
#'
#'
wts_fun = function(nr,k,n, sort = TRUE){
  # Formula on page 3 of Hill (2012)
  # Create vectors
    rs = w = rep(NA,nrow(nr))
  # Rank distances
    rd = rank(nr$d, ties.method = "first")
  # Rank similarities
    rs = rank(-nr$s, ties.method = "first")
  # Determine which of the k nearest distances are in the top n in terms of similarity
  # Note sort to keep original order
    top_s = sort(order(nr$s,decreasing = TRUE)[1:n])
  # subset rd and rs only keeping ones in top_s
    rd = rd[top_s]
    rs = rs[top_s]
  # Calc wts based on ranks
    w = (1 - (rd-1)^2/k^2)^4 * (1 - (rs-1)^2/n^2)^4
  # Build data.frame to return with indices of d/s that correspond to top_s
    out_obj = data.frame(location2 = nr$location[top_s], w = w, row.names = NULL)
  # If sort == TRUE Reorder by weights
  if(sort){
    out_obj = out_obj[order(out_obj$w, decreasing = TRUE),]
    row.names(out_obj) = NULL
  }
  # Return weights
    return(out_obj)
}



#' Determine neighbourhood weights for a subset of locations
#'
#' @inheritParams calc_neigh_wts
#' @param foc_inds Indices of rows in `comp_param` that specify the locations in
#'   the dataset for which neighbourhood weights are to be calculated
#'
#' @return a `data.frame` with the neighbourhood locations and corresponding
#'   weights for each focal location (i.e. those in rows `foc_ind` in
#'   `comp_param`. The `data.frame` has following three columns `location1`,
#'   `location2` and `w` which contain the identifier of the focal location, the
#'   identifier of locations in it's neighbourhood and the corresponding
#'   weighting respectively.
#'
#' @export
#'
#' @description This function calculates the neighbourhood weights for a subset
#'   of locations. Typically user will instead use the convenient wrapper
#'   function `calc_neigh_wts()` to calculate the neighbourhood weights for all
#'   locations in `comp_param`.
#'
#'
#' @seealso [calc_neigh_wts()] for a wrapper function that uses the current
#'   function to calculates the neighbourhood weights for all locations in
#'   `comp_param`. calling this function.
#'
#' @examples
#'
#'
batch_neigh_wts = function(foc_inds,comp_param, k = 200, n = 100){
  # If k is >= to rows in comp_param then stop
    if(k >= nrow(comp_param)) warning("specified k is >= to number of row/locations in 'comp_param'")
  # Determine distances between foc_loc and comp_locs
    dist_mat = dist_sf(comp_param[foc_inds,],comp_param)
  # Determine columns in comp_param
    par_cnames = names(sf::st_drop_geometry(comp_param))
    par_cnames = par_cnames[!grepl("location",par_cnames)]
  # Setup list to hold output
    temp_wts = vector("list",length(foc_inds))

  # Loop through foc_inds and determine k nearest sites (including self)
  for(i in 1:length(foc_inds)){
    # Determine indices of k nearest locs
      nr_locs = order(dist_mat[i,])[1:k]
    # Create data.frame to store neigh_stats
      nr_stats = data.frame(location = comp_param$location[nr_locs], d = NA, s = NA)
    # Extract distances for these
      nr_stats[,"d"] = as.numeric(dist_mat[i,nr_locs])
      #nr_d = dist_mat[i,nr_locs]
    # Now calculate similarities between these
      nr_stats[,"s"] = cosine_sim(
        sf::st_drop_geometry(comp_param[nr_locs,par_cnames]),which(nr_locs == foc_inds[i])
      )[1,] # First element of nr_locs should be foc_inds[i] but use which to be sure
      #nr_sim = cosine_sim(sf::st_drop_geometry(comp_param[nr_locs,par_cnames]),foc_inds[i])[1,]
    # Determine weights & add focal location to produce output in weights format
      temp_wts[[i]] = data.frame(
        location1 = comp_param$location[foc_inds[i]],
        wts_fun(nr_stats, k, n),
        row.names = NULL
      )
  }

  # Now collapse temp_wts list to a single data.frame
    out_wts = do.call("rbind",temp_wts)
  # Return weights
    return(out_wts)
}

#' Calculate Neighbourhood Weights for frescaloo based on distance and
#' environmental or biological variables
#'
#' @param comp_param The set of variables upon which the similarity between
#'   locations will be determined. The object need to be an `sf` data.frame with
#'   locations represented as point or polygon geometries and the variables used
#'   for the similarity calculation being columns in the associated attribute
#'   table. The first column of the attribute table should be titled `location`
#'   and should give the location identifiers that correspond to the occupancy
#'   data you wish to analyse with frescalo.
#' @param k The number of locations closest in distance to the focal location
#'   for each neighbourhood that will be used in the next stage of neighbourhood
#'   selection (i.e. similarity testing). The default `k = 200` will use the 200
#'   closest locations to the focal location in the similarity calculations for
#'   each neighbourhood.
#' @param n The number of locations that are most similar to the focal location
#'   selected to form the final neighbourhood. These are selected from the `k`
#'   closest locations.
#' @param in_parallel Logical determining whether the function will run in
#'   parallel using foreach across a series of workers or will run sequentially.
#'   The default `in_parallel = TRUE` will run the function in parallel using a
#'   previously registered cluster or creating a new one as required.
#' @param n_cores The number of cores/worker processes to spawn when running in
#'   parallel. This only has an effect if a new cluster is being created by the
#'   function otherwise it will use the settings from the previously registered
#'   cluster.
#' @param chunkSize The number of locations in each chunk or batch run by a
#'   worker process.
#'
#' @return a `data.frame` with the neighbourhood locations and corresponding
#'   weights for each focal location. The `data.frame` has following three
#'   columns `location1`, `location2` and `w` which contain the identifier of
#'   the focal location, the identifier of locations in it's neighbourhood and
#'   the corresponding weighting respectively.
#'
#' @export
#'
#' @description A function to produce a neighbourhood weights dataset that can
#'   be used in frescalo analyses. For each location a neighbourhood will be
#'   defined as nearby locations that are most similar based on values from a
#'   set of variables. The process is two stage, first the `k` nearest locations
#'   to a focal location, including the focal location itself, are determined
#'   and then from these `K` nearest location the `n` most similar are selected
#'   as the neighbourhood and will have weights calculated based on the
#'   distances and similarities.
#'
#'   The similarities are determined based on a set of variables provided with
#'   values for each location. It may be worth normalising each variable in
#'   `comp_param`.
#'
#'
#' @seealso [batch_neigh_wts()] for the function used internally to calculate
#'   weights for a subsets of locations, [cosine_sim()] for the function used to
#'   calculate similarity, [wts_fun()] for the function to calculating weights
#'   from distances & similarity measures, [frescalo()] for function to preform
#'   the frescalo analysis.
#'
#' @examples
#'
calc_neigh_wts = function(
    comp_param,
    k = 200,
    n= 100,
    in_parallel = TRUE,
    n_cores = parallel::detectCores(logical = FALSE)-1,
    chunkSize = 250
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

  # Split comp_param into chunks/batches
  # Determine number locations in comp_param
    n_locs = length(unique(comp_param$location))
  # Determine number of chunks/batches needed to ensure all locs are analysed
    n_chunk = ceiling(n_locs/chunkSize)
  # Build group vector allocating locations to chunks/batches
    loc_grps = rep(1:n_chunk, each = chunkSize)[1:n_locs]
  # Produce list of row indices for each group (each element vector of IDs for that chunk)
    cnk_locs = split(1:n_locs,loc_grps)

  # using foreach in parallel or sequentially to create neighbourhood weightings for each chunk
  if(in_parallel){
    wts_out = foreach(i = 1:n_chunk, .inorder=T, .combine='rbind', .multicombine=TRUE) %dopar% {
      devtools::load_all() # lines need for local testing remove for installed package
      batch_neigh_wts(foc_inds = cnk_locs[[i]], comp_param = comp_param, k = k, n = n)
    }
  } else {
    wts_out = foreach(i = 1:n_chunk, .inorder=T, .combine='rbind', .multicombine=TRUE) %do% {
      batch_neigh_wts(foc_inds = cnk_locs[[i]], comp_param = comp_param, k = k, n = n)
    }
  }

  # Now return
    return(wts_out)
}

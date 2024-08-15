dist_sf = function(sf_foc, sf_comp){
  # Convert to centroid
    pt_foc = sf::st_centroid(sf::st_geometry(sf_foc))
    pt_comp = sf::st_centroid(sf::st_geometry(sf_comp))
  # Return distance matrix
    pt_d = sf::st_distance(pt_foc, pt_comp)
  # Convert to numerical object
    return(pt_d)
}


# Modified from post in https://stats.stackexchange.com/questions/31565/compute-a-cosine-dissimilarity-matrix-in-r
cosine_sim = function(df, foc_inds = NULL){
  df_mat = as.matrix(df)
  sim = df_mat / sqrt(rowSums(df_mat*df_mat))
  if(!is.null(foc_inds)){
    sim = sim[foc_inds,] %*% t(sim)
  } else {
    sim = sim %*% t(sim)
  }
  return(sim)
}

wts_fun = function(nr,k,n){
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
  # subset rd and rsonly keeping ones in top_s
    rd = rd[top_s]
    rs = rs[top_s]
  # Calc wts based on ranks
    w = (1 - (rd-1)^2/k^2)^4 * (1 - (rs-1)^2/n^2)^4
  # Build data.frame to return with indices of d/s that correspond to top_s
    out_obj = data.frame(location2 = nr$location[top_s], w = w, row.names = NULL)
  # Return weights
    return(out_obj)
}


# comp_param is sf points or polyg data.frame with columns: location, geometry and the rest being params sim to be based on
# Note due to using order with ties it is possible that slighly different neighbourhoods will be returned if order of locations is changed (perhaps implement a sort to standardise, on coords??)
# Calculate neighbourhood & wts
calc_neigh_wts = function(foc_inds,comp_param, k = 200, n = 100){
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

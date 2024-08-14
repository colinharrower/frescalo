dist_sf = function(sf_foc, sf_comp){
  # Convert to centroid
    pt_foc = sf::st_centroid(sf_foc)
    pt_comp = sf::st_centroid(sf_comp)
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
    # Now calculate similarities between these
      nr_sim = cosine_sim(sf::st_drop_geometry(comp_param[nr_locs,par_cnames]),foc_inds[i])[1,]
    # Now build data.frame of results
      temp_wts[[i]] = data.frame(location1 = comp_param$location[foc_inds[i]], location2 = comp_param$location[nr_locs], w = nr_sim,row.names = NULL)
  }

}

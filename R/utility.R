speciesListFun = function(spList, species){
  # A function that returns the presences and absences of species per region
  # spList should be a data frame with columns location, species and time

  # Create list of all locations
  locationList = unique(spList$location)

  # Initialize output list
  out = lapply(c(1:length(locationList)), function(i) {i})

  for (i in 1:length(locationList)) {
    spListSub = subset(spList, location==locationList[i])
    out[[i]] = as.integer(species %in% spListSub$species)
  }
  return(out)
}

cfun = function(...) {
  # Bespoke function to combine the output from the frescalo() function
  input_list <- list(...)
  return(list(locs=Reduce('rbind',Map(function(x){x[[1]]},input_list)),
              freq=Reduce('rbind',Map(function(x){x[[2]]},input_list))))
}

cfunTrend = function(...) {
  # Another bespoke function to combine the output from the trend() function
  input_list <- list(...)
  return(list(trend=Reduce('rbind',Map(function(x){x[[1]]},input_list)),
              site_time=Reduce('rbind',Map(function(x){x[[2]]},input_list))))
}

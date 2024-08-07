# Can this function work or be modified to be parallelised over more combinations/nodes than just the number of time periods

# Addtional function for creating P_ijt values
predict_occ <- function(freq, trend){
  list_data <- by(data = freq, as.factor(freq$location), function(x,trend){(1-(1-x$freq_1)^trend[trend$species==x$species,]$tFactor)}, trend = trend)
  # Initialize variables
  species <- vector("list", length(list_data))
  IDs <- vector("list", length(list_data))
  Values <- vector("list", length(list_data))
  Labels <- vector("list", length(list_data))
  labels <- unique(trend$time) # Label vector
  # Populate the variables with data from the list (can this be done as vector calcs rather than loop?)
  for (i in seq_along(list_data)) {
    species[[i]] <- rep(unique(freq$species), times = length(list_data[[i]]))
    IDs[[i]] <- rep(names(list_data)[i], times = length(list_data[[i]]))
    Values[[i]] <- list_data[[i]]
    Labels[[i]] <- labels
  }
  # Combine all elements into vectors
  species <- unlist(species)
  IDs <- unlist(IDs)
  Values <- unlist(Values)
  Labels <- unlist(Labels)
  df <- data.frame(species = species, location = IDs, time = Labels, P_ijt = Values)
  return(df)
}

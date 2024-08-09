#' Estimate frescalo tFactors for occupancy data based on corrected
#' neighbourhood frequencies
#'
#' @param s_data occupancy dataset. Supplied as a data.frame with the following
#'   named columns; `time`, `location`, `species`. The dataset should contain a
#'   row corresponding to the recording or presence of the specified species at
#'   the specified location during the specified time period.
#' @param freq corrected neighbourhood frequencies returned by `fresc_scaling`.
#'   A data.frame containing the following columns; `location` `species`,
#'   `freq_1`,`rank_` and `benchmark`.
#' @param all_spp character vector containing all species identifiers. The
#'   default `all_spp = unique(s_data)` will determine this from `s_data`. If
#'   running in batches, i.e. in parallel, the list from the full dataset can
#'   instead be supplied overiding this behaviour
#'
#' @return a list with two elements each containing a data.frame. The
#'   first element `trend` cotains a data.frame with the tFactors for each
#'   species. This table is equiavlent to the `trend.out` file produced by the
#'   original fortan frescalo program. The second element `site_time` contains
#'   site x time period recording effort information used to make site/time
#'   period effort adjustments.
#'
#' @export
#'
#' @description Fits the second part of frescalo analyses which utilises the
#'   corrected neighbourhood frequencies to determine differences in frequency
#'   between the time periods.
#'
#'#' @seealso [frescalo()] for a convenience wrapper function that applies both
#' stages of the frescalo analyses or [fresc_scaling()] for the function to determine
#' rescaling factors and corrected frequencies required by `fresc_trends()`.
#'
#' @examples
#'
#' ### First part of this example is taken from the example for fresc_scaling()
#'
#' ## Extract results for a couple of focal locations from included
#' ## unicorns example dataset
#'
#' # Get subset of 250 locations from this dataset
#'   test_locs = unique(s$location)[1:250]
#'   test_spp = unique(s$species)
#'
#' # Extract data from weights and occ data
#'   test_wts = d[which(d$location1 %in% test_locs),]
#'   test_occ = s[which(s$location %in% test_locs),]
#' # Convert to list of prescence/absence
#'   pa_lt = frescalo:::speciesListFun(test_occ, species = test_spp)
#'
#' ## Now determine rescaling parameter and adjust frequencies using
#' ## using fresc_scaling()
#'
#'  out_freq = fresc_scaling(
#'   test_wts,pa_lt,
#'   all_loc = test_locs,
#'   all_spp = test_spp
#'  )
#'
#'  ## Now determine tFactor for first time period
#'  out_trend = fresc_trend(
#'   test_occ[which(test_occ$time == 1),],
#'   out_freq$freq,
#'   test_spp
#'  )
#'
#'  ## View results
#'  head(out_trend[["trend"]])
#'  head(out_trend[["site_time"]])
#'
fresc_trend = function(s_data, freq, all_spp = unique(s_data)) {
  # Frescalo trend analysis
  # Function to calculate the time scaling factor for each species in a single year bin
  #
  # s_data are the records from a single time period. Each row is a reocrd of a species.
  # The data frame needs at least 3 columns:
  #      location  an identifier for the location of the record
  #      species   an identifier for the species
  #      time      an identifier for the time bin
  #
  # freq is a data frame with the corrected neighbourhood frequencies from a FRESCALO analysis (see the frescalo function).
  # This data frame should have columns:
  #      location  an identifier for the location
  #      species   an identifier for the species
  #      freq_1    the corrected neighbourhood frequency
  #      rank_1    the rank of the corrected neighbourhood frequency
  #      benchmark =1 if this species should be used as a benchmark species

  # Check that the input data (s_data) is from a single time bin
  timeBin = unique(s_data$time)
  if (length(timeBin)>1) {
    warning('More than one time bin supplied to trend() function')
  }
  #freq <- output$freq
  #sSplit2 = split(s, as.factor(s$time))
  #s_data = sSplit2[[1]]
  #Species 5 in time period has zero TF, check how I cna keep this
  locationList = as.character(unique(freq$location))
  #spList = unique(s_data$species)
  # OLP: use full list of species so that time-period specific zeros are kept
  #spList = unique(s$species) # TODO fix this so that it is not using objects from the  not passed directly to the f

  # Calculate the proportion of benchmark species in each hectad (for this time bin)
  focal_s = split(s_data, factor(s_data$location, levels=locationList))
  focal_bench = split(freq, factor(freq$location, levels=locationList))
  s_it = mapply(FUN=function(x,y) {sum(x$benchmark[x$species %in% y$species])/sum(x$benchmark)},
                x=focal_bench, y=focal_s, SIMPLIFY=T, USE.NAMES=F)
  s_it[is.na(s_it)] = 1.0E-7 # added by OL Pescott Apr 2024, as apparently sometimes x$benchmark can be zero (and so s_it is NaN)
  # Calculate weights to downweight infrequenctly sampled locations
  # i.e. seeing fewer than 9.95% of the benchmark species (s_it<0.0995)
  w = rep(1, each=length(locationList))
  w[s_it<0.0995] = 10*s_it[s_it<0.0995]+0.005
  # Add site x time period recording effort info to output
  site_time = data.frame(
    location = locationList,
    time = rep(unique(s_data$time), times = length(locationList)),
    s_it = s_it,
    w = w
  )

  focal_s2 = split(s_data, factor(s_data$species, levels=all_spp))
  sumP_ijtw =  unlist(lapply(FUN=function(X){sum(w[locationList%in%X$location])}, X=focal_s2), use.names=F)

  focal_f = split(freq, factor(freq$species, levels=all_spp))
  sf = lapply(FUN=function(x) {x$freq_1[match(locationList,x$location)]*s_it}, focal_f)

  # Calculate Q_ijt. If sf>0.98 Set Q_ijt=-log(1-0.98)=3.912023
  # P_ijt = 1-exp(-Q_ijt x_jt) where x_jt is the time factor for species j at time t
  # and P_ijt is prob of observing species j in hectad i at time t
  Q_ijt = lapply(FUN=function(x){y=-log1p(-x); y[y>3.912023]=3.912023; return(y)}, X=sf)

  StDev = x = xSD = rep(NA, times=length(all_spp))   # Vector to contain the time factors
  sptot1 = estvar = rep(0, times=length(all_spp))   # Vector to contain the time factor standard deviations
  for (i in 1:length(all_spp)) {
    # Rescale frequencies by effort for all hectads
    # tmp = subset(freq, species==all_spp[s])
    # sf_tmp = tmp$freq_1[match(locationList,tmp$location)]*s_it

    # # Calculate Q_ijt
    # # P_ijt = 1-exp(-Q_ijt x_jt) where x_jt is the time factor for species j at time t
    # # and P_ijt is prob of observing species j in hectad i at time t
    # Q_ijt = rep(-log(1-0.98), each=length(locationList))
    # Q_ijt[sf[[s]]<0.98] = -log(1-sf[[s]][sf[[s]]<0.98])

    # Select a x max that ensures a sign change in min_trend_fun
    if (any(Q_ijt[[i]]>0)) {
      xMax = 5 # Yearsley original
      #xMax = 1
      while(min_trend_fun(xMax, Q_ijt[[i]], w, sumP_ijtw[i])<0) {xMax = xMax+5}
      sol = uniroot(min_trend_fun,interval=c(0,xMax), tol=0.0005, Q_ijt[[i]], w, sumP_ijtw[i])
      x[i] = sol$root
      #x = sol$root
      #i = 52
      # Calculate time factor SDs
      estvar[i] = sum((1-exp(-Q_ijt[[i]]*x[i])) * (1-(1-exp(-Q_ijt[[i]]*x[i]))) * w * w)
      sptot1[i] = sumP_ijtw[i] + sqrt(estvar[i])
      xSD[i] = sptot1[i]/(sumP_ijtw[i]+0.0000001)
      while(min_trend_fun(xSD[i], Q_ijt[[i]], w, sptot1[i])<0) {xSD[i]  = xSD[i] + 1}
      # Small value added to xSD[i] to accommodate species with no data in a time bin (stops uniroot failing)
      sol2 = uniroot(min_trend_fun,interval=c(0,xSD[i]+1.0E-10), tol=0.0005, Q_ijt[[i]], w, sptot1[i])
      xSD[i] = sol2$root
      StDev[i] = abs(xSD[i] - x[i])
    }
  }
  df = data.frame(species=all_spp, time=timeBin, tFactor=x, StDev = StDev, estvar=estvar, sptot1=sptot1)
  df = df[order(df$species, df$time),]
  # Reset row.names
  row.names(df) = NULL
  return(list(trend = df, site_time = site_time))
}


min_trend_fun = function(x, Q, w, sumPw) {
  # A function to solve sum_i P_ijt * w_i = sum_i (1-exp(-Q_ijt*x))*w_i in the Frescalo trend analysis Hill (2012)
  # P_ijt = 1-exp(-Q_ijt*x)
  # w_i is a weighting factor to downweight area with low fraction of benchmark species (s_it<0.095)
  # P_ijt = probability of observing species j in region i at time t
  return(sum((1-exp(-Q*x))*w) - sumPw)
}

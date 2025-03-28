#' Determine Frequency Scaling Factor Based on Local Occupancy
#'
#' @param neigh_wts neighbourhoods weighs data. A data.frame with the following
#'   named columns `location1`, `location2`, `w` containing 1) location
#'   identifier for the focal location of the neighbourhood, 2) location
#'   identifiers for locations in the neighbourhood for the current focal
#'   location, 3) weighting values specifying the similarity between the focal
#'   location and the specific location within it's neighbourhood.
#' @param spp_pa list specifying presence or absence of species at each
#'   location. The list should correspond to `all_loc` where the ith element of
#'   `spp_pa` list is the presence/absence data for ith location in `all_loc`.
#'   Each element of `spp_pa` should contain a vector of 0 or 1 values of length
#'   = `length(all_spp)` where 0 represents absence and 1 presence of that given
#'   species in the given location. The order of the presence and absence data
#'   should match the order of species in `all_spp`
#' @param all_loc Vector containing all unique locations in the occupancy data.
#'   The order should match the order of locations in `spp_pa`.
#' @param all_spp vector containing all unique taxa or species names in the
#'   occupancy data. The order should match the order of species presence and
#'   absence data in `spp_pa`.
#' @param Phi The target neighbourhood frequency. The default is `Phi = 0.74`.
#' @param R_star The benchmarking species rank proportion. The default `R_star =
#'   0.2703` means that after correction for recording effort the top 27% of
#'   species are used as benchmarks.
#' @param filter_wts Logical value determining whether neigh_wts should be
#'   filtered to only use neighbourhood info for locations in the occupancy
#'   data. The default `filter_wts = FALSE` will attempt to determine statistics
#'   for every focal location in `neigh_wts` while `TRUE` will exclude any focal
#'   locations in `neigh_wts` that are not appear in the occupancy data.
#' @param missing_data action to take when locations in neighbourhoods are
#'   missing from the occupancy dataset. The default `missing_data = 2` assumes
#'   that all species are absent or not recorded from the missing locations,
#'   which . The alternative `missing_data = 1` will exclude any neighbourhoods
#'   that contain locations missing from the occupancy data.
#'
#' @return A list comprised of two elements each containing a data.frame. The
#'   first element `locs` contains location metrics while the the second `freq`
#'   contains the original and re-scaled species frequencies. These two outputs
#'   are equivalent to `samples.txt` and `frequencies.out` output files produced
#'   by Mark Hill's original frescalo fortran program.
#'
#'   The data.frame in `locs` has a row for each location in `spp_pa` with the
#'   following columns: \item{location}{the location name(s) or identifier(s)}
#'   \item{nSpecies}{the total number of species or taxa observed}
#'   \item{phi_in}{the `phi_i`, i.e. local frequency-weighted mean frequency,
#'   for the local neighbourhood for the given location}
#'   \item{alpha}{the
#'   rescaling factor used to adjust the resepective neighbourhood `phi` so that
#'   it matches the target `phi` specied in the `Phi` argument to the function}
#'   \item{phi_out}{the final `phi` after rescaling, which should match the
#'   value supplied to the input argument `Phi `}
#'   \item{spnum_in}{}
#'   \item{iter}{the number of iterations required to determine the scaling
#'   factor `alpha` required}
#'
#'   The data.frame in `freq` has rows for each location taxa combination with
#'   the following columns: \item{location}{the location name(s) or
#'   identifier(s)} \item{species}{the species/taxa name(s) or identifier(s)}
#'   \item{pres}{logical determining prescence (`1`) or absence (`0`) of the
#'   species/taxa at that location}
#'   \item{freq}{the localisied frequency for that species & neighbourhood
#'   before rescaling}
#'   \item{freq_1}{the local frequency for that species &
#'   neighbourhood after rescaling}
#'   \item{rank}{the original frequency rank} \item{rank_1}{the frequency rank
#'   after rescaling}
#'   \item{benchmark}{indicates whether the species was a benchmark
#'   species for the local neighbourhood (`1` = benchmark species/taxa)}
#'
#' @export
#'
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#'
#' @description Determines and applies frequency scaling factors in order to
#'   attempt to correct for variation in recording effort in a biological
#'   records dataset.
#'
#' @seealso [frescalo()] for a convenience wrapper function that applies both
#' stages of the frescalo analyses or [fresc_trends()] for the function to determine
#' time factors for a given time period using the outputs of `fresc_scaling()`.
#'
#' @examples
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
#' ## View outputs
#'   head(out_freq[["locs"]])
#'   head(out_freq[["freq"]])
#'
#'
#'
fresc_scaling = function(neigh_wts, spp_pa, all_loc, all_spp, Phi=0.74, R_star=0.2703, filter_wts = FALSE, missing_data = 2){
  # Determine all unique focal locations in neigh_wts and filter to locations in spp_pa if necessary
    foc_locs = unique(neigh_wts$location1)
    n_foc = length(foc_locs)
    # Determine how many of these are not in the occupancy data
    nd_inds = foc_locs[which(!foc_locs %in% all_loc)]
    n_nd = length(nd_inds)
    # If filter_wts TRUE the reduce neigh_wts to only rows where location1 is in all_loc
    if(filter_wts){
      neigh_wts = neigh_wts[which(neigh_wts$location1 %in% all_loc),]
      # Redetermine foc_locs and n_foc after removal
      foc_locs = unique(neigh_wts$location1)
      n_foc = length(foc_locs)
    } else {
      if(n_nd > 10 & (n_nd/n_foc) >= 0.5)
      warning(paste0("A large proportion of focal locations in 'neight_wts' are not present in 'spp_pa'\n"))
    }
    # Cleanup
    rm(nd_inds, n_nd)

  # Calculate no. of spp for each focal location
    nSpecies = sapply(spp_pa[match(foc_locs,all_loc,nomatch=length(all_loc)+1)], FUN=sum, simplify=T) # CH length(all_loc)+1 for nomatch so now if foc_loc not in all_loc then nSpecies = 0 should be returned (though due to filtering of weights dataset I don't think this should ever actually happen)
  # Create data.frame to hold output location summary statistics for focal locs
    out_loc = data.frame(
      location=foc_locs,
      nSpecies=nSpecies,
      phi_in=NA, alpha=NA, phi_out=Phi, spnum_in=NA, spnum_out=NA, iter=NA
    )

  # Setup temporary object to hold location species frequency statistics (use temp list and collapse afterwards to avoid growing object)
  temp_freq = vector("list",n_foc)
  #freq.out = data.frame(location=c(), species=c(), pres=c(), freq=c(), rank=c(), rank_1=c(), benchmark=c())

  for(i_f in 1:n_foc){
    focal = foc_locs[i_f]
    focal_d = neigh_wts[which(neigh_wts$location1==focal),]
    #cat("DEBUG: i_f = ",i_f," - ",focal,"\n",sep="")

    # Identify species in neighbourhood of focal region
    # If data missing for any neighbourhood locations for current focal location assign (e.g. nomatch) give it a index 1 greater than length of (spp_pa) so it will have a NULL value in speciesRegional
    neighbourhood = match(focal_d$location2, all_loc, nomatch=length(spp_pa)+1) # empty neighbourhood sites are NULL
    speciesRegional = spp_pa[neighbourhood]
    # Find any locations in the neighbourhood that are not in spp_pa (should be NULL in speciesRegional) and create a spp_pa like vector for them with zeros for every species
    speciesRegional = lapply(speciesRegional, function(x) if(is.null(x)) rep(0, times = length(all_spp)) else x)

    missingData = neighbourhood==length(spp_pa)+1 # length(spp_pa) is number of locations; missing data is neighbourhoods with
    if (all(missingData) | (any(missingData) & missing_data==1)){
      if(all(missingData)){
        # Species data missing from a location in the neighbourhood. Ignore this focal location
        warning(paste('Removing location ',focal,'. No data for any location in the neighbourhood.', sep=''))
      } else {
        # Species data missing from a location in the neighbourhood. Ignore this focal location
        warning(paste('Removing location ',focal,'. Missing data for locations in the neighbourhood.', sep=''))
      }

      out_loc$alpha[i_f] = NA
      out_loc$iter[i_f] = NA
      out_loc$phi_in[i_f]= NA
    } else {
      # Calculate weights of locations in the neighbourhood
      weights = focal_d$w/(sum(focal_d$w)+1.0E-10)

      # Create weighted neighbourhood frequencies (checked against Frescalo)
      frequency = Reduce('+',Map('*',as.list(weights), speciesRegional))
      # Added by Oli Pescott, Apr 2024, to avoid zeros and match Hill (2012) outputs
      frequency = ifelse(frequency==0, frequency+1.0E-10, frequency)
      phi_in = sum(frequency^2) / sum(frequency)

      # Calculate the multiplier (alpha) that equalises recording effort
      # Modified from Jon/Oli's original code to try and deal with situations with sparse data
      # and unfiltered neight_wts object were resulting in infinite while loops. There is probably
      # a better way to deal with this but for now seems to work (batches of 10 and break if no improvement)
      alpha_min = 1  # Minimum alpha ( =1 means no correction required)
      alpha_max = 5

      cur_alpha_rng = seq(alpha_max,by = 5,length.out = 10)
      cur_min_f = -1
      # Increase alpha_max until freq_min_f() becomes positive (i.e. ensure there is a zero)
      while (all(cur_min_f < 0)) {
        cur_min_f = sapply(cur_alpha_rng, freq_min_f, fij = frequency, Phi = Phi)
        if(any(cur_min_f > 0)){
          alpha_max = cur_alpha_rng[which.max(cur_min_f > 0)]
        } else {
          if(length(unique(cur_min_f)) == 1){
            break
          }
          cur_alpha_rng = seq(max(cur_alpha_rng)+5,by = 5,length.out = 10)
        }
      }
      # Look for cur_min_f all being the same value
      if(length(unique(cur_min_f)) == 1){
        warning(paste('Removing location ',focal,'. Changing alpha_max had no effect.', sep=''))

        # Set values to NA
        out_loc$alpha[i_f] = NA
        out_loc$iter[i_f] = NA
        out_loc$phi_in[i_f]= NA

        # Skip
        next
      }

      cur_alpha_rng = rep(alpha_min,10)*0.5^(0:9)
      cur_min_f = 1
      while (all(cur_min_f > 0)) {
        cur_min_f = sapply(cur_alpha_rng, freq_min_f, fij = frequency, Phi = Phi)
        if(any(cur_min_f < 0)){
          alpha_min = cur_alpha_rng[which.max(cur_min_f < 0)]
        }
        if(length(unique(cur_min_f)) == 1){
          break
        }
        cur_alpha_rng = rep(min(cur_alpha_rng)/2,10)*0.5^(0:9)
      }
      if(length(unique(cur_min_f)) == 1){

        warning(paste('Removing location ',focal,'. Changing alpha_min had no effect.', sep=''))
        # Set values to NA
        out_loc$alpha[i_f] = NA
        out_loc$iter[i_f] = NA
        out_loc$phi_in[i_f]= NA

        # Skip
        next
      }
      #while (freq_min_f(alpha_min, frequency, Phi)>0) { alpha_min = alpha_min/2}

      # Find sampling-effort multiplier
      sol=uniroot(freq_min_f,interval=c(alpha_min,alpha_max), tol=0.0003, frequency, Phi)
      out_loc$alpha[i_f] = sol$root
      out_loc$phi_in[i_f] = phi_in
      out_loc$iter[i_f] = sol$iter
      # Expected species richness before recorder effort correct
      out_loc$spnum_in[i_f] = sum(frequency)
      # Expected species richness after recorder effort correct
      out_loc$spnum_out[i_f] = sum(1-exp(sol$root*log(1-frequency)))

      # Create the data frame with the local species frequencies after correction
      freq_ord = order(frequency, decreasing=T) # no explicit method to deal with ties, so any ties just arranged using original name order
      foc_ind = match(focal,all_loc,nomatch=length(all_loc))

      # Pick out benchmark species (assumes that species are ordered by rank)
      benchmarkSpecies = rep(0, times=length(frequency))
      R_prime = c(1:length(spp_pa[[foc_ind]]))/out_loc$spnum_out[i_f] # rescaled rank
      # Benchmark species are either those where rescaled rank < R_star, or are ranked number 1 even tho all R_prime > R_star
      benchmarkSpecies[R_prime<R_star | c(1:length(spp_pa[[foc_ind]]))==1 ] = 1
      #benchmarkSpecies[R_prime<R_star] = 1 ## Yearsley original only had one condition. Adding second increases correlation with Hill fortran marginally

      # freq_1 is the corrected neighbourhood frequencies
      temp_freq[[i_f]] = data.frame(
        location=focal,
        species=all_spp[freq_ord],
        pres=spp_pa[[foc_ind]][freq_ord],
        freq=frequency[freq_ord],
        freq_1=1-exp(sol$root*log(1-frequency[freq_ord])),
        rank=c(1:length(spp_pa[[foc_ind]])),
        rank_1=R_prime,
        benchmark=benchmarkSpecies
      )
    }
  }
  # Collapse temp_freq to single data.frame
    out_freq = do.call("rbind",temp_freq)
  # Return outputs as list
  return(list(locs=out_loc, freq=out_freq))
}

#' Determine the difference between the rescaled neighbourhod frequency and the target Phi
#' for a given value of alpha
#'
#' @param alpha the value of alpha
#' @param fij the frequencies of each species in a given neighbourhood
#' @param Phi The target Phi value
#'
#' @return The difference between the target Phi and the current rescaled Neighbourhood
#'   frequency. The alpha value at which this equals zero indicates the
#'   scaling factor required for that neighbourhood to reach the target Phi.
#'
#' @description
#' This function is iteratively called within `fresc_scaling` when attempting to determine the appropriate alpha values for each neighbourhood.
#' The goal is to find the alpha value that minimises the value returned by this function (i.e. where it reaches zero).
#'
#'
#'
#' @examples
#'
#' ## Create a set of frequency values
#'   test_freq = runif(1000,0,1)
#'
#' ## Run freq_min_f on these test_frequencies
#'   freq_min_f(5,test_freq,0.74)
#'
#'
freq_min_f <- function(alpha,fij,Phi){
  # The function to minimise to fit Frescalo model to the data
  # Page 5 column 2 in Hill (2011)
  fr=1-(1-fij)^alpha
  #fr[abs(fij-1)<1e-10] = 1 # Yearsley
  fr[abs(fij-1)<1.0E-10] = 0.99999 # to match Hill (2012), although not needed for Yearlsey calc strictly
  fr[abs(fij-1)>0.99999] = 1.0E-10
  out = sum(fr^2)/sum(fr) - Phi
  return(out)
}

#' Predict Species Occupancies across locations and time
#'
#' @param freq the frequency results `data.frame` returned from `frescalo()`
#' @param trend the trend results `data.frame` returned from `frescalo()`
#'
#' @return a data.frame containing the predicted occupancies (`p_occ`) for the
#'   location, species, time combinations
#' @export
#'
#' @description The outputs from the frescalo model are used to predict
#' occupancies for each species, location, time combination.
#'
#'
#' @examples
#' \dontrun{
#' ## Run frescalo on the included example dataset
#'   fres_res = frescalo(s,d,filter_wts = TRUE)
#'
#' ## Predict occupancies for all species,location,time combinations
#'   p_occ = predict_occ(fres_res$freq, fres_res$trend)
#' }
#'
#'
predict_occ <- function(freq, trend){
  list_data <- by(
    data = freq, freq[,c("location","species")],
    function(x,trend){
      cur_t = trend[trend$species==x$species,]
      cur_p = (1-(1-x$freq_1)^cur_t$tFactor)
      ret_obj = data.frame(species = x$species, location = x$location, time = cur_t$time, p_occ = cur_p)
      return(ret_obj)
    }, trend = trend, simplify = FALSE)
  # Collapse list
    df = do.call("rbind",list_data)
  # Now return
  return(df)
}

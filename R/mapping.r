#' Map output from frescalo model
#'
#' @param fres_sf A frescalo results table as an `sf` object where the locations
#'   have been assigned with sf geometries (e.g. `POLYGON`, `POINT`). If the
#'   object passed to `fres_sf` has a single non-geometry column that is the
#'   column that will be used to colour the points/polygons, otherwise you will
#'   need to specify the name of the specific column to be used for this purpose
#'   in `zcol`
#' @param zcol The name of the column in `fres_sf` to be used for the plotting.
#'   The default `zcol = NULL` assumes that `fres_sf` will have a single
#'   non-geometry column that will be used, otherwise an error will be returned
#'   unless you specify a column name.
#' @param over_bd An `sf` geometry only object that you wish to overlay on top
#'   of the map (e.g. a country boundary). The default `over_bd = NULL` assumes
#'   no overlay.
#' @param over_args A list of additional arguments that will be passed to the
#'   `sf plot()` function used to add `over_bd` to the map (if `over_bd` is not
#'   `NULL`). The default `over_args = NULL` assumes no additional arguments
#' @param ... Additional arguments that will be passed to the `sf plot` function
#'   used to plot `fres_sf`. For instance `breaks = "quantile"` can be used to
#'   changes the breaks style to quantile based breaks.
#'
#' @return NULL
#' @export
#'
#' @description This plotting fun
#'
#' @seealso [sf][plot.sf()]
#'
#' @examples
fres_map = function(fres_sf,zcol = NULL, over_bd = NULL, over_args = NULL, ...){
  # If zcol specified only keep that column otherwise check is only 1 non-geometry column
  if(is.null(zcol)){
    stopifnot(ncol(sf::st_drop_geometry(sf_loc[,"spnum_out"])) == 1)
  } else {
    fres_sf = fres_sf[,zcol]
  }

  # Plot using plot.sf with reset = FALSE so that things can be added after the fact
    plot(fres_sf, reset = FALSE, ...)

  # If over_bd specified then add to e
  if(!is.null(over_bd)){
    do.call(plot, args = c(list(x = over_bd, add = TRUE), over_args))
  }
}

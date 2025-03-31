#' Plot tFactors for a single species
#'
#' @param sp_trend results from `trend` output of frescalo for a single species.
#' @param ylab The y axis label
#' @param xlab The x axis label
#' @param ylim The limits for the y axis. The default ylim = NULL will set these
#'   based on the range of the TFactors +- SD
#' @param ... Any other arguments to be passed to `plot()` and or `arrows()`
#'   functions used to create the plot or add the error bars
#'
#' @description Produces a time factor (tFactor) plot for a given species. The
#'   time factor plot shows the tFactor estimated for each time period. Error
#'   bars are used to indicate the Standard Deviation around each tFactor.
#'
#'   Note this function expects trend data for a single species and will produce
#'   a single plot for that species. If you want to plot tFactor plots for all
#'   species you need write code to open a plotting device and then loop or run
#'   through all species subsetting the data and fitting this function then
#'   closing the plotting device.
#'
#' @return NULL
#' @export
#'
#' @examples
#' # Fit frescalo
#' out_fres = frescalo(s,d,in_parallel = FALSE,filter_wts = TRUE)
#'
#' # Plot tFactor for species "Species 1"
#' plot_tfactor(out_fres$trend[which(out_fres$trend$species == "Species 1"),])
#'
plot_tfactor = function(sp_trend, ylab = "Relative Frequency", xlab = "Time", ylim = NULL, ...){
  if(is.null(ylim)){
    ylim = range(c(sp_trend$tFactor-sp_trend$StDev,sp_trend$tFactor+sp_trend$StDev))
  }
  plot(sp_trend[,c("time","tFactor")],ylim = ylim, main = sp_trend$species[1], ylab = ylab, xlab = xlab, xaxt = "n",...)
  axis(1, at = sort(unique(sp_trend$time)))
  arrows(x0=sp_trend$time, x1=sp_trend$time, y0=sp_trend$tFactor-sp_trend$StDev, y1=sp_trend$tFactor+sp_trend$StDev, code = 3, angle = 90,length = 0.05,...)
}

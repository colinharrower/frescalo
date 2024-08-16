#' Plot tFactor for a species
#'
#' @param sp_trend results from `trend` output of frescalo for a single species.
#' @param ylab The y axis label
#' @param xlab The x axis label
#' @param ... Any other arguments to be passed to `plot()` and or `arrows()`
#'   functions used to create the plot or add the error bars
#'
#' @return NULL
#' @export
#'
#' @examples
plot_tfactor = function(sp_trend, ylab = "Relative Frequency", xlab = "Time", ...){
  y_lims = range(c(sp_trend$tFactor-sp_trend$StDev,sp_trend$tFactor+sp_trend$StDev))
  plot(sp_trend[,c("time","tFactor")],ylim = y_lims, main = sp_trend$species[1], ylab = ylab, xlab = xlab, xaxt = "n",...)
  axis(1, at = sort(unique(sp_trend$time)))
  arrows(x0=sp_trend$time, x1=sp_trend$time, y0=sp_trend$tFactor-sp_trend$StDev, y1=sp_trend$tFactor+sp_trend$StDev, code = 3, angle = 90,length = 0.05,...)
}

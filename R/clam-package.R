#' clam
#'
#' Classical (non-Bayesian) age-depth modelling.
#'
#' @docType package
#' @aliases clam-package
"_PACKAGE"
#' @author Maarten Blaauw <maarten.blaauw@qub.ac.uk>
#' @description Performs 'classical' age-depth modelling of dated sediment deposits - prior to applying more sophisticated techniques such as Bayesian age-depth modelling. Any radiocarbon dated depths are calibrated. Age-depth models are constructed by sampling repeatedly from the dated levels, each time drawing age-depth curves. Model types include linear interpolation, linear or polynomial regression, and a range of splines. See Blaauw (2010). <doi:10.1016/j.quageo.2010.01.002>.
#' @importFrom grDevices dev.off grey pdf png rgb
#' @importFrom graphics abline image layout legend lines par points polygon rect plot
#' @importFrom stats approx density dnorm lm loess pnorm predict qnorm quantile rnorm runif smooth.spline spline weighted.mean
#' @importFrom utils packageName read.csv read.table write.table
#' @importFrom data.table fread fwrite
#' @importFrom rice hpd
#' @name clam
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL


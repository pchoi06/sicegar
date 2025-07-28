# parameterCalculation_h0 Double Sigmoidal Helper Functions
#' @title Log‐transformed Double‐Sigmoidal Core Function
#' @description Computes the log of the core double‐sigmoidal form used for
#'   optimizing the peak time.  This function is internal to the SPT parameter
#'   calculation.
#' @param x Numeric vector of time points.
#' @param B1 Numeric; slope of the rising phase.
#' @param M1 Numeric; midpoint of the rising phase.
#' @param B2 Numeric; absolute slope of the declining phase.
#' @param L  Numeric; distance between the two steepest points.
#' @return Numeric vector of \code{log(1/((1+e^{-B1(x-M1)})*(1+e^{B2(x-(M1+L))})))}.
#' @keywords internal
f1 <- function (x, B1, M1, B2, L) # untouched
{
  log((1/((1 + exp(-B1 * (x - M1))) * (1 + exp(B2 * (x - (M1 +
                                                            L)))))))
}


#' Half-Maximum Midpoint Function for Double-Sigmoidal Curve
#'
#' Given the parameters of the double-sigmoidal backbone \(f_0(x)\) and the peak
#' value `const`, returns \(f_0(x) - \tfrac{1}{2}\,const\).  This is used with
#' `uniroot()` to find the time points where the curve crosses half its maximum.
#'
#' @param x Numeric vector of time points at which to evaluate.
#' @param B1 Numeric; the first slope parameter of the double-sigmoidal backbone.
#' @param M1 Numeric; the first midpoint parameter (rise center).
#' @param B2 Numeric; the second slope parameter of the double-sigmoidal backbone.
#' @param L Numeric; the distance between the two midpoints.
#' @param const Numeric; the peak (maximum) value of \(f_0(x)\).
#' @return Numeric vector of the same length as `x`, giving \(\,f_0(x) - \tfrac{1}{2}\,const\).
#' @noRd
f0mid <- function (x, B1, M1, B2, L,const) {
  -const / 2 + 1 / ( ( 1 + exp(-B1 * (x-M1)) ) * ( 1 + exp(B2*(x-(M1+L))) ) )
}


# Argmax
#' @title Argmax for Double‐Sigmoidal Peak Time
#' @description Finds the time (within the scaled range) at which the fitted
#'   double‐sigmoidal curve reaches its maximum.
#' @param parameterVector A named list or data frame containing:
#'   \describe{
#'     \item{\code{slope1Param_N_Estimate}}{Rising slope (raw).}
#'     \item{\code{slope2Param_N_Estimate}}{Declining slope (raw).}
#'     \item{\code{midPoint1Param_N_Estimate}}{Rising midpoint (raw).}
#'     \item{\code{midPointDistanceParam_N_Estimate}}{Distance between midpoints (raw).}
#'     \item{\code{dataScalingParameters.timeRange}}{Total time range scaling factor.}
#'   }
#' @return Numeric scalar: the time (in original units) of the curve’s peak.
#' @keywords internal
f_argmax_doublesigmoidal_h0 <- function (parameterVector)
{
  slope1Param <- parameterVector$slope1Param_N_Estimate
  slope2Param <- parameterVector$slope2Param_N_Estimate
  midPointDistanceParam <- parameterVector$midPointDistanceParam_N_Estimate
  finalAsymptoteIntensityRatio <- parameterVector$finalAsymptoteIntensityRatio_N_Estimate
  maximum <- parameterVector$maximum_N_Estimate
  midPoint1Param <- parameterVector$midPoint1Param_N_Estimate
  timeRange <- parameterVector$dataScalingParameters.timeRange
  if (slope1Param < 0) {
    stop("slope1Param should be a positive number")
  }
  if (slope2Param < 0) {
    stop("slope2Param should be a positive number. It is the absolute value of the slope2Param")
  }
  if (midPointDistanceParam < 0) {
    stop("midPointDistanceParam should be a positive number. It is the distance between two steppest points of exponential phase and lysis")
  }
  #if (finalAsymptoteIntensityRatio < 0 | finalAsymptoteIntensityRatio >
  #  1) {
  #  stop("finalAsymptoteIntensityRatio should be a number between 0 and 1")
  #}
  #if (maximum < 0) {
  #  stop("maximum should be a positive number")
  #}
  optimizeIntervalMin <- midPoint1Param - 2 * midPointDistanceParam
  optimizeIntervalMax <- midPoint1Param + 3 * midPointDistanceParam
  xmax <- stats::optimize(f1, c(optimizeIntervalMin, optimizeIntervalMax),
                          tol = 1e-04, B1 = slope1Param, M1 = midPoint1Param,
                          B2 = slope2Param, L = midPointDistanceParam, maximum = TRUE)
  argumentt <- xmax$maximum * timeRange
  return(argumentt)
}


# Mid1
#' Find First Half-Maximum Point of Double-Sigmoidal Fit
#'
#' Given a parameter data frame from a double-sigmoidal fit, finds the
#' first time point (x) at which the model reaches half of its maximum
#' response by solving for the root of the mid-point equation.
#'
#' @param parameterDf A data.frame or list containing at least:
#'   - `dataScalingParameters.timeRange`
#'   - `slope1Param_Estimate`
#'   - `midPoint1Param_Estimate`
#'   - `slope2Param_Estimate`
#'   - `midPointDistanceParam_Estimate`
#' @return A numeric scalar giving the time of the first half-maximum point.
#' @importFrom stats optimize uniroot
#' @export
f_mid1_doublesigmoidal_h0 <- function (parameterDf)
{
  max_x <- parameterDf$dataScalingParameters.timeRange
  xmax <- stats::optimize(f1, interval = c(-1.125 * max_x,
                                           max_x * 3), tol = 1e-04, B1 = parameterDf$slope1Param_Estimate,
                          M1 = parameterDf$midPoint1Param_Estimate, B2 = parameterDf$slope2Param_Estimate,
                          L = parameterDf$midPointDistanceParam_Estimate, maximum = TRUE)
  argumentt <- xmax$maximum
  constt <- f0(argumentt, B1 = parameterDf$slope1Param_Estimate,
               M1 = parameterDf$midPoint1Param_Estimate, B2 = parameterDf$slope2Param_Estimate,
               L = parameterDf$midPointDistanceParam_Estimate)
  mid1x <- try(expr = stats::uniroot(f0mid, interval = c(-1.125 *
                                                           max_x, argumentt), tol = 1e-04, B1 = parameterDf$slope1Param_Estimate,
                                     M1 = parameterDf$midPoint1Param_Estimate, B2 = parameterDf$slope2Param_Estimate,
                                     L = parameterDf$midPointDistanceParam_Estimate, const = constt),
               silent = TRUE)
  if (!inherits(mid1x, "try-error"))  {
    rangeList <- seq(-1, 30)
    rangeListCounter <- -1
    while (class(mid1x) == "try-error" & rangeListCounter <=
           30) {
      rangeListSelection <- rangeList[rangeListCounter +
                                        2]
      mid1x <- try(expr = stats::uniroot(f0mid, interval = c((-1) *
                                                               (2^rangeListSelection) * max_x, argumentt),
                                         tol = 1e-04, B1 = parameterDf$slope1Param_Estimate,
                                         M1 = parameterDf$midPoint1Param_Estimate, B2 = parameterDf$slope2Param_Estimate,
                                         L = parameterDf$midPointDistanceParam_Estimate,
                                         const = constt), silent = TRUE)
      rangeListCounter <- rangeListCounter + 1
    }
  }
  if (inherits(mid1x, "try-error")) {
    stop("f_mid1_doublesigmoidal_h0(): unable to find first half-max root with uniroot()")
  }
  return(mid1x$root)
}


# Mid2
#' Find Second Half-Maximum Point of Double-Sigmoidal Fit
#'
#' For a double-sigmoidal fit, locates the time at which the response
#' declines back to half-maximum, by root-finding on the mid-point function
#' over the interval from the peak to beyond the original time range.
#'
#' @param parameterDf A data.frame or list containing the same fields as
#'   \code{f_mid1_doublesigmoidal_h0()}, including \code{dataScalingParameters.timeRange}.
#' @return A numeric scalar giving the time of the second half-maximum point.
#' @importFrom stats optimize uniroot
#' @export
f_mid2_doublesigmoidal_h0 <- function (parameterDf)
{
  max_x <- parameterDf$dataScalingParameters.timeRange
  xmax <- stats::optimize(f1, interval = c(-1.125 * max_x,
                                           max_x * 3), tol = 1e-04, B1 = parameterDf$slope1Param_Estimate,
                          M1 = parameterDf$midPoint1Param_Estimate, B2 = parameterDf$slope2Param_Estimate,
                          L = parameterDf$midPointDistanceParam_Estimate, maximum = TRUE)
  argumentt <- xmax$maximum
  constt <- f0(argumentt, B1 = parameterDf$slope1Param_Estimate,
               M1 = parameterDf$midPoint1Param_Estimate, B2 = parameterDf$slope2Param_Estimate,
               L = parameterDf$midPointDistanceParam_Estimate)
  mid2x <- try(stats::uniroot(f0mid, interval = c(argumentt,
                                                  max_x * (1.25)), tol = 1e-04, B1 = parameterDf$slope1Param_Estimate,
                              M1 = parameterDf$midPoint1Param_Estimate, B2 = parameterDf$slope2Param_Estimate,
                              L = parameterDf$midPointDistanceParam_Estimate, const = constt),
               silent = TRUE)
  if (!inherits(mid2x, "try-error")) {
    rangeList <- seq(-1, 30)
    rangeListCounter <- -1
    while (class(mid2x) == "try-error" & rangeListCounter <=
           30) {
      rangeListSelection <- rangeList[rangeListCounter +
                                        2]
      mid2x <- try(stats::uniroot(f0mid, interval = c(argumentt,
                                                      max_x * (2^rangeListSelection)), tol = 1e-04,
                                  B1 = parameterDf$slope1Param_Estimate, M1 = parameterDf$midPoint1Param_Estimate,
                                  B2 = parameterDf$slope2Param_Estimate, L = parameterDf$midPointDistanceParam_Estimate,
                                  const = constt), silent = TRUE)
      rangeListCounter <- rangeListCounter + 1
    }
  }
  if (inherits(mid2x, "try-error")) {
    stop("f_mid2_doublesigmoidal_h0(): unable to find second half-max root with uniroot()")
  }
  return(mid2x$root)
}


# Slope
#' @title Numerical Slope for Double‐Sigmoidal at a Given Time
#' @description Approximates the first derivative of the fitted
#'   double‐sigmoidal curve at time \code{x} using a 5‐point stencil.
#' @param x Numeric scalar; time point at which to compute the slope.
#' @param parameterDf A named list or data frame of fitted parameters, including:
#'   \describe{
#'     \item{\code{finalAsymptoteIntensityRatio_Estimate}}{Final ratio.}
#'     \item{\code{maximum_Estimate}}{Maximum intensity.}
#'     \item{\code{slope1Param_Estimate}}{Rising slope parameter.}
#'     \item{\code{midPoint1Param_Estimate}}{Rising midpoint.}
#'     \item{\code{slope2Param_Estimate}}{Declining slope parameter.}
#'     \item{\code{midPointDistanceParam_Estimate}}{Distance between midpoints.}
#'     \item{\code{h0_Estimate}}{Baseline intensity.}
#'   }
#' @param timeStep Numeric; small increment used for finite differences. Default 1e-05.
#' @return Numeric scalar: approximate slope (dY/dX) at \code{x}.
#' @keywords internal
f_slope_doublesigmoidal_h0 <- function (x, parameterDf, timeStep = 1e-05)
{
  fxp2h <- doubleSigmoidalFitFormula_h0(x = x + 2 * timeStep,
                                        finalAsymptoteIntensityRatio = parameterDf$finalAsymptoteIntensityRatio_Estimate,
                                        maximum = parameterDf$maximum_Estimate, slope1Param = parameterDf$slope1Param_Estimate,
                                        midPoint1Param = parameterDf$midPoint1Param_Estimate,
                                        slope2Param = parameterDf$slope2Param_Estimate, midPointDistanceParam = parameterDf$midPointDistanceParam_Estimate, h0 = parameterDf$h0_Estimate)
  fxph <- doubleSigmoidalFitFormula_h0(x = x + timeStep, finalAsymptoteIntensityRatio = parameterDf$finalAsymptoteIntensityRatio_Estimate,
                                       maximum = parameterDf$maximum_Estimate, slope1Param = parameterDf$slope1Param_Estimate,
                                       midPoint1Param = parameterDf$midPoint1Param_Estimate,
                                       slope2Param = parameterDf$slope2Param_Estimate, midPointDistanceParam = parameterDf$midPointDistanceParam_Estimate, h0 = parameterDf$h0_Estimate)
  fxmh <- doubleSigmoidalFitFormula_h0(x = x - timeStep, finalAsymptoteIntensityRatio = parameterDf$finalAsymptoteIntensityRatio_Estimate,
                                       maximum = parameterDf$maximum_Estimate, slope1Param = parameterDf$slope1Param_Estimate,
                                       midPoint1Param = parameterDf$midPoint1Param_Estimate,
                                       slope2Param = parameterDf$slope2Param_Estimate, midPointDistanceParam = parameterDf$midPointDistanceParam_Estimate, h0 = parameterDf$h0_Estimate)
  fxm2h <- doubleSigmoidalFitFormula_h0(x = x - 2 * timeStep,
                                        finalAsymptoteIntensityRatio = parameterDf$finalAsymptoteIntensityRatio_Estimate,
                                        maximum = parameterDf$maximum_Estimate, slope1Param = parameterDf$slope1Param_Estimate,
                                        midPoint1Param = parameterDf$midPoint1Param_Estimate,
                                        slope2Param = parameterDf$slope2Param_Estimate, midPointDistanceParam = parameterDf$midPointDistanceParam_Estimate, h0 = parameterDf$h0_Estimate)
  der <- (-1 * fxp2h + 8 * fxph - 8 * fxmh + 1 * fxm2h)/(12 *
                                                           timeStep)
  return(der)
}


#' @title Compute Intuitively Meaningful Curve Parameters
#' @description From the raw fit vector of either a sigmoidal or double‐sigmoidal
#'   model, calculates intuitively interpretable quantities (e.g. peak times,
#'   slopes at midpoints, plateau heights) and appends them to the parameter list.
#' @param parameterVector A named list or data frame produced by
#'   \code{\link{sigFit}} or \code{\link{sigFit2}}, containing at minimum:
#'   \describe{
#'     \item{\code{model}}{Either \dQuote{sigmoidal} or \dQuote{doublesigmoidal}.}
#'     \item{\code{*_Estimate}}{Raw parameter estimates (renormalized back to the original scale).}
#'     \item{\code{dataScalingParameters.timeRange},
#'           \code{dataScalingParameters.intensityRange},
#'           \code{dataScalingParameters.intensityMin}}{
#'       Scaling parameters used to normalize/de‐normalize the data.}
#'   }
#' @param stepSize Numeric; finite‐difference increment for slope calculation. Default is \code{1e-05}.
#' @return The input \code{parameterVector}, augmented with new fields:
#'   for sigmoidal: \code{midPoint_x}, \code{midPoint_y}, \code{slope}, etc.;
#'   for double‐sigmoidal: \code{maximum_x}, \code{midPoint1_x}, \code{slope1},
#'   \code{finalAsymptoteIntensity}, etc.
#' @export
parameterCalculation_h0 <- function (parameterVector, stepSize = 1e-05)
{
  if (parameterVector$model == "sigmoidal") {
    parameterVector$maximum_x <- NA
    parameterVector$maximum_y <- parameterVector$maximum_Estimate
    parameterVector$midPoint_x <- parameterVector$midPoint_Estimate
    parameterVector$midPoint_y <- parameterVector$maximum_y/2
    parameterVector$h0_y      <- parameterVector$h0_Estimate # changed
    parameterVector$slope <- parameterVector$slopeParam_Estimate *
      parameterVector$maximum_y * (1/4)
    parameterVector$incrementTime <- parameterVector$maximum_y/parameterVector$slope
    parameterVector$startPoint_x <- parameterVector$midPoint_x -
      (parameterVector$incrementTime/2)
    parameterVector$startPoint_y <- parameterVector$h0_y # changed
    parameterVector$reachMaximum_x <- parameterVector$midPoint_x +
      (parameterVector$incrementTime/2)
    parameterVector$reachMaximum_y <- parameterVector$maximum_y
    parameterVector$additionalParameters <- TRUE
  }
  if (parameterVector$model == "doublesigmoidal") {
    parameterVector$maximum_x <- f_argmax_doublesigmoidal_h0(parameterVector)
    parameterVector$maximum_y <- parameterVector$maximum_Estimate
    parameterVector$midPoint1_x <- f_mid1_doublesigmoidal_h0(parameterVector)
    parameterVector$midPoint1_y <- doubleSigmoidalFitFormula_h0(x = parameterVector$midPoint1_x,
                                                                finalAsymptoteIntensityRatio = parameterVector$finalAsymptoteIntensityRatio_Estimate,
                                                                maximum = parameterVector$maximum_y, slope1Param = parameterVector$slope1Param_Estimate,
                                                                midPoint1Param = parameterVector$midPoint1Param_Estimate,
                                                                slope2Param = parameterVector$slope2Param_Estimate,
                                                                midPointDistanceParam = parameterVector$midPointDistanceParam_Estimate,
                                                                h0 = parameterVector$h0_Estimate)
    parameterVector$midPoint2_x <- f_mid2_doublesigmoidal_h0(parameterVector)
    parameterVector$midPoint2_y <- doubleSigmoidalFitFormula_h0(x = parameterVector$midPoint2_x,
                                                                finalAsymptoteIntensityRatio = parameterVector$finalAsymptoteIntensityRatio_Estimate,
                                                                maximum = parameterVector$maximum_y, slope1Param = parameterVector$slope1Param_Estimate,
                                                                midPoint1Param = parameterVector$midPoint1Param_Estimate,
                                                                slope2Param = parameterVector$slope2Param_Estimate,
                                                                midPointDistanceParam = parameterVector$midPointDistanceParam_Estimate,
                                                                h0 = parameterVector$h0_Estimate)
    parameterVector$slope1 <- f_slope_doublesigmoidal_h0(parameterVector$midPoint1_x, # changed
                                                         parameterVector, timeStep = stepSize)
    parameterVector$slope2 <- f_slope_doublesigmoidal_h0(parameterVector$midPoint2_x, # changed
                                                         parameterVector, timeStep = stepSize)
    parameterVector$finalAsymptoteIntensity <- parameterVector$finalAsymptoteIntensityRatio_Estimate *
      parameterVector$maximum_y
    parameterVector$incrementTime <- parameterVector$maximum_y/parameterVector$slope1
    parameterVector$startPoint_x <- parameterVector$midPoint1_x -
      (parameterVector$incrementTime/2)
    parameterVector$startPoint_y <- parameterVector$h0_Estimate # changed
    parameterVector$reachMaximum_x <- parameterVector$midPoint1_x +
      (parameterVector$incrementTime/2)
    if (parameterVector$reachMaximum_x > parameterVector$maximum_x) {
      parameterVector$reachMaximum_x <- parameterVector$maximum_x
      parameterVector$warning.reachMaximum_cor = TRUE
    }
    parameterVector$reachMaximum_y <- parameterVector$maximum_y
    parameterVector$decrementTime <- (parameterVector$maximum_y -
                                        parameterVector$finalAsymptoteIntensity)/(-parameterVector$slope2)
    parameterVector$startDeclinePoint_x <- parameterVector$midPoint2_x -
      (parameterVector$decrementTime/2)
    if (parameterVector$startDeclinePoint_x < parameterVector$maximum_x) {
      parameterVector$startDeclinePoint_x <- parameterVector$maximum_x
      parameterVector$warning.startDeclinePoint_cor = TRUE
    }
    parameterVector$startDeclinePoint_y <- parameterVector$maximum_y
    parameterVector$endDeclinePoint_x <- parameterVector$midPoint2_x +
      (parameterVector$decrementTime/2)
    parameterVector$endDeclinePoint_y <- parameterVector$finalAsymptoteIntensity
    parameterVector$additionalParameters = TRUE
  }
  return(parameterVector)
}

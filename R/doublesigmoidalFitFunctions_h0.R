#' Double‐sigmoidal base function
#'
#' A helper “base” function used internally by the double‐sigmoidal fit formula.
#'
#' @param x Numeric vector of time points.
#' @param B1 Positive slope parameter for the first rise.
#' @param M1 Midpoint of the first rise.
#' @param B2 Positive slope parameter for the second fall.
#' @param L  Distance between the two midpoints.
#' @return Numeric vector of the same length as `x`.
#' @keywords internal
f0 <- function (x, B1, M1, B2, L) #unchanged
{
  (1/((1 + exp(-B1 * (x - M1))) * (1 + exp(B2 * (x - (M1 + L))))))
}

#' Heaviside Piecewise double‐sigmoidal
#'
#' Internal helper that stitches together two scaled copies of `f0()` at a
#' switch point.
#'
#' @param x        Numeric vector of time points.
#' @param A2       Final asymptote ratio.
#' @param Ka       Maximum plateau height.
#' @param B1       First slope (see `f0()`).
#' @param M1       First midpoint.
#' @param B2       Second slope.
#' @param L        Midpoint distance.
#' @param const    Normalizing constant (value of `f0()` at its peak).
#' @param argument Switch point on the x‐axis.
#' @param h0       Baseline intensity.
#' @return Numeric vector of same length as `x`.
#' @keywords internal
f2_h0 <- function (x, A2, Ka, B1, M1, B2, L, const, argument, h0)
{
  fBasics::Heaviside(x - argument) * (f0(x, B1, M1, B2, L) * ((Ka - A2 * Ka)/(const)) + A2 * Ka) +
    (1 - fBasics::Heaviside(x - argument)) * (f0(x, B1, M1, B2, L) * ((Ka - h0)/(const)) + h0)
}

#' Convert normalized min.pack estimated parameters back to raw scale (with h0)
#'
#' Takes a one‐row parameter data.frame and “unscales” its
#' time‐normalized and intensity‐normalized columns.
#'
#' @param parameterDF Data.frame with columns `*_N_Estimate` and
#'   `dataScalingParameters.*`.
#' @param isalist     Logical; if `TRUE` input came from a list,
#'   otherwise already raw.
#' @return A one‐row data.frame with `*_Estimate` columns on the raw
#'   scale (including `h0_Estimate`).
#' @export
#' @examples
#' df <- data.frame(
#'   model="doublesigmoidal",
#'   maximum_N_Estimate=1, finalAsymptoteIntensityRatio_N_Estimate = .5,
#'   slope1Param_N_Estimate=1, midPoint1Param_N_Estimate = .3,
#'   slope2Param_N_Estimate=1, midPointDistanceParam_N_Estimate = .2,
#'   h0_N_Estimate = .1,
#'   dataScalingParameters.timeRange = 100,
#'   dataScalingParameters.intensityRange = 10,
#'   dataScalingParameters.intensityMin = 0
#' )
#' renormalizeParameters_h0(df, isalist=TRUE)
renormalizeParameters_h0 <- function(parameterDF, isalist) {
  model         <- parameterDF$model
  dataInputName <- parameterDF$dataInputName
  if (isalist) {
    # raw plateau heights
    h1_raw <- parameterDF$maximum_N_Estimate  * parameterDF$dataScalingParameters.intensityRange + parameterDF$dataScalingParameters.intensityMin
    h2_raw <- (fAIR_n <- parameterDF$finalAsymptoteIntensityRatio_N_Estimate * parameterDF$maximum_N_Estimate) * parameterDF$dataScalingParameters.intensityRange + parameterDF$dataScalingParameters.intensityMin
    parameterDF$finalAsymptoteIntensityRatio_Estimate <-
      h2_raw / h1_raw
    #parameterDF$finalAsymptoteIntensityRatio_Estimate <-
    #  parameterDF$finalAsymptoteIntensityRatio_N_Estimate
    parameterDF$maximum_Estimate <-
      parameterDF$maximum_N_Estimate *
      parameterDF$dataScalingParameters.intensityRange +
      parameterDF$dataScalingParameters.intensityMin
    parameterDF$slope1Param_Estimate <-
      parameterDF$slope1Param_N_Estimate /
      parameterDF$dataScalingParameters.timeRange
    parameterDF$midPoint1Param_Estimate <-
      parameterDF$midPoint1Param_N_Estimate *
      parameterDF$dataScalingParameters.timeRange
    parameterDF$slope2Param_Estimate <-
      parameterDF$slope2Param_N_Estimate /
      parameterDF$dataScalingParameters.timeRange
    parameterDF$midPointDistanceParam_Estimate <-
      parameterDF$midPointDistanceParam_N_Estimate *
      parameterDF$dataScalingParameters.timeRange
    parameterDF$h0_Estimate <-
      parameterDF$h0_N_Estimate *
      parameterDF$dataScalingParameters.intensityRange +
      parameterDF$dataScalingParameters.intensityMin
  }
  if (!isalist) {
    # if already raw, just copy over
    parameterDF$finalAsymptoteIntensityRatio_Estimate <-
      parameterDF$finalAsymptoteIntensityRatio_N_Estimate
    parameterDF$maximum_Estimate <-
      parameterDF$maximum_N_Estimate
    parameterDF$slope1Param_Estimate <-
      parameterDF$slope1Param_N_Estimate
    parameterDF$midPoint1Param_Estimate <-
      parameterDF$midPoint1Param_N_Estimate
    parameterDF$slope2Param_Estimate <-
      parameterDF$slope2Param_N_Estimate
    parameterDF$midPointDistanceParam_Estimate <-
      parameterDF$midPointDistanceParam_N_Estimate
    parameterDF$h0_Estimate <-
      parameterDF$h0_N_Estimate
  }
  parameterDF$model         <- model
  parameterDF$dataInputName <- dataInputName
  return(parameterDF)
}

#' Double‐sigmoidal fit formula with baseline h0
#'
#' Evaluates the full 5‐parameter “rise+fall” model at `x`.
#'
#' @param x                         Numeric vector.
#' @param finalAsymptoteIntensityRatio Ratio of second plateau to first.
#' @param maximum                   Maximum plateau height.
#' @param slope1Param               First slope (>0).
#' @param midPoint1Param            First midpoint (>0).
#' @param slope2Param               Second slope (>0).
#' @param midPointDistanceParam     Distance to second midpoint (>0).
#' @param h0                        Baseline when x is negative infinity.
#' @return Numeric vector of same length as `x`.
#' @export
doubleSigmoidalFitFormula_h0 <- function (x, finalAsymptoteIntensityRatio, maximum, slope1Param,
                                          midPoint1Param, slope2Param, midPointDistanceParam, h0)
{
  if (slope1Param < 0) {
    stop("slope1Param should be a positive number")
  }
  if (slope2Param < 0) {
    stop("slope2Param should be a positive number. It is the absolute value of the second slopeParam")
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
  xmax <- stats::optimize(
    f1,
    c(optimizeIntervalMin, optimizeIntervalMax),
    tol = 1e-4,
    B1 = slope1Param,
    M1 = midPoint1Param,
    B2 = slope2Param,
    L  = midPointDistanceParam,
    maximum = TRUE
  )
  argumentt <- xmax$maximum
  constt <- f0(argumentt, slope1Param, midPoint1Param, slope2Param,
               midPointDistanceParam)
  y <- f2_h0(x, finalAsymptoteIntensityRatio, maximum, slope1Param,
             midPoint1Param, slope2Param, midPointDistanceParam,
             constt, argumentt, h0)
  return(y)
}

#' Double‐sigmoidal NLS fit (with baseline h0)
#'
#' Fits `doubleSigmoidalFitFormula_h0()` to a `data.frame(time, intensity)`
#' via `minpack.lm::nlsLM()`, appending h0 and AIC/BIC.
#'
#' @param dataInput      data.frame or list with `$timeIntensityData`.
#' @param tryCounter     Integer: which random‐start iteration this is.
#' @param startList      Named list of starting values (including h0).
#' @param lowerBounds    Named vector of lower bounds.
#' @param upperBounds    Named vector of upper bounds.
#' @param min_Factor     Passed to `nlsLM()`.
#' @param n_iterations   Max iterations for `nlsLM()`.
#' @return A one‐row data.frame of estimates, fit‐statistics, and
#'   `h0_Estimate`.
#' @export
#' @importFrom minpack.lm nlsLM
sigFit2 <- function (
    dataInput,
    tryCounter,
    startList = list(
      finalAsymptoteIntensityRatio = 0.5,
      maximum                      = 1,
      slope1Param                  = 1,
      midPoint1Param               = 0.33,
      slope2Param                  = 1,
      midPointDistanceParam        = 0.29,
      h0                           = 0.1     # ← added
    ),

    lowerBounds = c(
      finalAsymptoteIntensityRatio = 0,
      maximum                      = 0.3,
      slope1Param                  = 0.01,
      midPoint1Param               = -0.52,
      slope2Param                  = 0.01,
      midPointDistanceParam        = 0.04,
      h0                           = 0    # ← added
    ),

    upperBounds = c(
      finalAsymptoteIntensityRatio = 1.5,
      maximum                      = 1.5,
      slope1Param                  = 180,
      midPoint1Param               = 1.15,
      slope2Param                  = 180,
      midPointDistanceParam        = 0.63,
      h0                           = 0.3    # ← added
    ),
    min_Factor   = 1/2^20,
    n_iterations = 1000
)
{
  isalist     <- (is.list(dataInput) & !is.data.frame(dataInput))
  if (isalist) {
    dataFrameInput <- dataInput$timeIntensityData
  }
  isadataframe <- is.data.frame(dataInput)
  if (isadataframe) {
    dataFrameInput <- dataInput
  }
  if (tryCounter == 1) {
    counterDependentStartList <- startList
  }
  if (tryCounter != 1) {
    randomVector <- stats::runif(length(startList), 0, 1)
    names(randomVector) <- names(startList)    # now includes "h0"
    counterDependentStartVector <- randomVector * (upperBounds - lowerBounds) + lowerBounds
    counterDependentStartList <- as.list(counterDependentStartVector)
  }

  theFitResult <- try(
    minpack.lm::nlsLM(
      intensity ~ doubleSigmoidalFitFormula_h0(
        time,
        finalAsymptoteIntensityRatio,
        maximum,
        slope1Param,
        midPoint1Param,
        slope2Param,
        midPointDistanceParam,
        h0
      ),
      data    = dataFrameInput,
      start   = counterDependentStartList,
      control = list(maxiter = n_iterations, minFactor = min_Factor),
      lower   = lowerBounds,
      upper   = upperBounds,
      trace   = FALSE,
    ),
    silent = TRUE
  )

  if (!inherits(theFitResult, "try-error")) {
    parameterMatrix <- summary(theFitResult)$parameters
    colnames(parameterMatrix) <- c("Estimate", "Std_Error", "t_value", "Pr_t")
    parameterVector <- c(t(parameterMatrix))
    # append 4 h0 names at the end
    names(parameterVector) <- c(
      "finalAsymptoteIntensityRatio_N_Estimate", "finalAsymptoteIntensityRatio_Std_Error",
      "finalAsymptoteIntensityRatio_t_value",    "finalAsymptoteIntensityRatio_Pr_t",
      "maximum_N_Estimate", "maximum_Std_Error", "maximum_t_value", "maximum_Pr_t",
      "slope1Param_N_Estimate","slope1Param_Std_Error","slope1Param_t_value","slope1Param_Pr_t",
      "midPoint1Param_N_Estimate","midPoint1Param_Std_Error","midPoint1Param_t_value","midPoint1Param_Pr_t",
      "slope2Param_N_Estimate","slope2Param_Std_Error","slope2Param_t_value","slope2Param_Pr_t",
      "midPointDistanceParam_N_Estimate","midPointDistanceParam_Std_Error","midPointDistanceParam_t_value","midPointDistanceParam_Pr_t",
      "h0_N_Estimate","h0_Std_Error","h0_t_value","h0_Pr_t"   # ← appended
    )
    parameterVector <- c(
      parameterVector,
      residual_Sum_of_Squares = sum((as.vector(stats::resid(theFitResult)))^2)[1],
      log_likelihood         = as.vector(stats::logLik(theFitResult))[1],
      AIC_value              = as.vector(stats::AIC(theFitResult))[1],
      BIC_value              = as.vector(stats::BIC(theFitResult))[1]
    )
    parameterList <- as.list(parameterVector)
    parameterList$isThisaFit <- TRUE
    parameterList$startVector <- counterDependentStartList
    if (isalist) parameterList$dataScalingParameters <- as.list(dataInput$dataScalingParameters)
    parameterList$model <- "doublesigmoidal"
    parameterList$additionalParameters <- FALSE
    parameterDf <- as.data.frame(parameterList)
    parameterDf <- renormalizeParameters_h0(parameterDf, isalist)
  } else {
    parameterVector <- rep(NA, 24 + 4)
    names(parameterVector) <- c(
      "finalAsymptoteIntensityRatio_N_Estimate","finalAsymptoteIntensityRatio_Std_Error",
      "finalAsymptoteIntensityRatio_t_value","finalAsymptoteIntensityRatio_Pr_t",
      "maximum_N_Estimate","maximum_Std_Error","maximum_t_value","maximum_Pr_t",
      "slope1Param_N_Estimate","slope1Param_Std_Error","slope1Param_t_value","slope1Param_Pr_t",
      "midPoint1Param_N_Estimate","midPoint1Param_Std_Error","midPoint1Param_t_value","midPoint1Param_Pr_t",
      "slope2Param_N_Estimate","slope2Param_Std_Error","slope2Param_t_value","slope2Param_Pr_t",
      "midPointDistanceParam_N_Estimate","midPointDistanceParam_Std_Error","midPointDistanceParam_t_value","midPointDistanceParam_Pr_t",
      "h0_N_Estimate","h0_Std_Error","h0_t_value","h0_Pr_t"   # ← appended
    )
    parameterVector <- c(
      parameterVector,
      residual_Sum_of_Squares = Inf,
      log_likelihood         = NA,
      AIC_value              = NA,
      BIC_value              = NA
    )
    parameterList <- as.list(parameterVector)
    parameterList$isThisaFit <- FALSE
    parameterList$startVector     <- counterDependentStartList
    if (isalist) parameterList$dataScalingParameters <- as.list(dataInput$dataScalingParameters)
    parameterList$model <- "doublesigmoidal"
    parameterDf <- as.data.frame(parameterList)
    parameterDf <- renormalizeParameters_h0(parameterDf, isalist)
  }

  return(parameterDf)  # unchanged
}

#' @title Sigmoidal Fit Formula with Baseline
#' @description Compute the sigmoidal curve value at given time points,
#'   allowing a nonzero baseline (\code{h0}).
#' @param x Numeric vector of time points.
#' @param maximum Numeric; the upper asymptote of the sigmoid.
#' @param slopeParam Numeric; the slope parameter (rate of rise).
#' @param midPoint Numeric; the time at which the sigmoid reaches halfway between \code{h0} and \code{maximum}.
#' @param h0 Numeric; the lower asymptote (baseline) intensity.
#' @return Numeric vector of model‐predicted intensities.
#' @examples
#' time <- seq(0, 100, length.out = 50)
#' y <- sigmoidalFitFormula_h0(time, maximum = 10, slopeParam = 0.1, midPoint = 50, h0 = 2)
#' @export
sigmoidalFitFormula_h0 <- function (x, maximum, slopeParam, midPoint, h0)
{
  y = (h0 + (maximum - h0)/(1 + exp((-slopeParam) * (x - midPoint))))
  return(y)
}

#' @title Renormalize Sigmoidal Fit Parameters
#' @description Transform normalized parameter estimates back to the
#'   original data scale, or pass through raw estimates unchanged.
#' @param parameterDF Data frame of parameter estimates with columns
#'   \code{*_N_Estimate} and, if present, \code{dataScalingParameters.*}.
#' @param isalist Logical; \code{TRUE} if \code{parameterDF} was returned
#'   from a list‐based fit (needs renormalization), \code{FALSE} otherwise.
#' @return Data frame with added or overwritten columns
#'   \code{h0_Estimate}, \code{maximum_Estimate}, \code{slopeParam_Estimate},
#'   and \code{midPoint_Estimate} on the original scale.
#' @export
sigmoidalRenormalizeParameters_h0 <- function(parameterDF, isalist) {
  model         <- parameterDF$model
  dataInputName <- parameterDF$dataInputName
  if (isalist) {
    parameterDF$maximum_Estimate    <- (parameterDF$maximum_N_Estimate    *
                                          parameterDF$dataScalingParameters.intensityRange) +
      parameterDF$dataScalingParameters.intensityMin
    parameterDF$slopeParam_Estimate <-  parameterDF$slopeParam_N_Estimate /
      parameterDF$dataScalingParameters.timeRange
    parameterDF$midPoint_Estimate   <-  parameterDF$midPoint_N_Estimate   *
      parameterDF$dataScalingParameters.timeRange
    parameterDF$h0_Estimate         <-  parameterDF$h0_N_Estimate         *
      parameterDF$dataScalingParameters.intensityRange +
      parameterDF$dataScalingParameters.intensityMin
  }
  if (!isalist) {
    parameterDF$maximum_Estimate    <- parameterDF$maximum_N_Estimate
    parameterDF$slopeParam_Estimate <- parameterDF$slopeParam_N_Estimate
    parameterDF$midPoint_Estimate   <- parameterDF$midPoint_N_Estimate
    parameterDF$h0_Estimate         <- parameterDF$h0_N_Estimate
  }
  parameterDF$model         <- model
  parameterDF$dataInputName <- dataInputName
  return(parameterDF)
}

#' @title Fit a Sigmoidal Model via Non‐linear Least Squares
#' @description Attempt to fit a 4‐parameter sigmoidal curve with baseline
#'   using \code{nlsLM}. Multiple random restarts may be used internally.
#' @param dataInput Data frame or named list containing columns \code{time} and \code{intensity}.
#' @param tryCounter Integer; if greater than 1, randomize starting values.
#' @param startList Named list of initial parameter values (\code{h0}, \code{maximum}, \code{slopeParam}, \code{midPoint}).
#' @param lowerBounds Named numeric vector of lower bounds for each parameter.
#' @param upperBounds Named numeric vector of upper bounds for each parameter.
#' @param min_Factor Numeric; minimal step‐factor passed to the optimizer.
#' @param n_iterations Integer; maximum iterations for \code{nlsLM}.
#' @return Data frame of parameter estimates and fit diagnostics, including:
#'   \describe{
#'     \item{\code{*_N_Estimate}}{Raw normalized estimates.}
#'     \item{\code{*_Estimate}}{Renormalized estimates on the original scale.}
#'     \item{\code{residual_Sum_of_Squares}, \code{AIC_value}, \code{BIC_value}}{
#'       Fit diagnostics: sum of squared residuals, Akaike’s IC, Bayesian IC, etc.}
#'     \item{\code{isThisaFit}}{Logical indicating successful convergence.}
#'   }
#' @examples
#' df <- data.frame(time = seq(0, 100, 10),
#'                  intensity = runif(11, 0, 10))
#' fit <- sigFit(df, tryCounter = 1)
#' @export
sigFit <- function(
    dataInput,
    tryCounter,
    startList = list(h0 = 0, maximum = 1, slopeParam = 1, midPoint = 0.33),
    lowerBounds = c(h0 = 0, maximum = 0.3, slopeParam = 0.01, midPoint = -0.52),
    upperBounds = c(h0 = 0.3, maximum = 1.5, slopeParam = 180,  midPoint = 1.15),
    min_Factor = 1/2^20,
    n_iterations = 1000
) {
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
  } else {
    randomVector <- stats::runif(length(startList), 0, 1)
    names(randomVector) <- c("h0", "maximum", "slopeParam", "midPoint")
    counterDependentStartVector <- randomVector * (upperBounds - lowerBounds) + lowerBounds
    counterDependentStartList <- as.list(counterDependentStartVector)
  }

  theFitResult <- try(
    minpack.lm::nlsLM(
      intensity ~
        h0 + (maximum - h0) / (1 + exp(-slopeParam * (time - midPoint))),
      data    = dataFrameInput,
      start   = counterDependentStartList,
      control = list(maxiter = n_iterations, minFactor = min_Factor),
      lower   = lowerBounds,
      upper   = upperBounds,
      trace   = FALSE
    ),
    silent = TRUE
  )

  if (!inherits(theFitResult, "try-error")) {
    parameterMatrix <- summary(theFitResult)$parameters
    colnames(parameterMatrix) <- c("Estimate", "Std_Error", "t_value", "Pr_t")

    parameterVector <- c(t(parameterMatrix))
    names(parameterVector) <- c(
      "h0_N_Estimate",    "h0_Std_Error",    "h0_t_value",    "h0_Pr_t",    # ← added h0
      "maximum_N_Estimate",    "maximum_Std_Error",    "maximum_t_value",    "maximum_Pr_t",
      "slopeParam_N_Estimate", "slopeParam_Std_Error", "slopeParam_t_value", "slopeParam_Pr_t",
      "midPoint_N_Estimate",   "midPoint_Std_Error",   "midPoint_t_value",   "midPoint_Pr_t"
    )

    parameterVector <- c(
      parameterVector,
      residual_Sum_of_Squares = sum((as.vector(stats::resid(theFitResult)))^2)[1],
      log_likelihood          = as.vector(stats::logLik(theFitResult))[1],
      AIC_value               = as.vector(stats::AIC(theFitResult))[1],
      BIC_value               = as.vector(stats::BIC(theFitResult))[1]
    )

    parameterList <- as.list(parameterVector)
    parameterList$isThisaFit  <- TRUE
    parameterList$startVector <- counterDependentStartList
    if (isalist) {
      parameterList$dataScalingParameters <- as.list(dataInput$dataScalingParameters)
    }
    parameterList$model                <- "sigmoidal"
    parameterList$additionalParameters <- FALSE

    parameterDf <- as.data.frame(parameterList)
    parameterDf <- sigmoidalRenormalizeParameters_h0(parameterDf, isalist)

  } else {
    naNames <- c(
      "h0_N_Estimate",    "h0_Std_Error",    "h0_t_value",    "h0_Pr_t",    # ← added h0
      "maximum_N_Estimate",    "maximum_Std_Error",    "maximum_t_value",    "maximum_Pr_t",
      "slopeParam_N_Estimate", "slopeParam_Std_Error", "slopeParam_t_value", "slopeParam_Pr_t",
      "midPoint_N_Estimate",   "midPoint_Std_Error",   "midPoint_t_value",   "midPoint_Pr_t"
    )
    parameterVector <- rep(NA, length(naNames))
    names(parameterVector) <- naNames

    parameterVector <- c(
      parameterVector,
      residual_Sum_of_Squares = Inf,
      log_likelihood          = NA,
      AIC_value               = NA,
      BIC_value               = NA
    )

    parameterList <- as.list(parameterVector)
    parameterList$isThisaFit  <- FALSE
    parameterList$startVector <- counterDependentStartList
    if (isalist) {
      parameterList$dataScalingParameters <- as.list(dataInput$dataScalingParameters)
    }
    parameterList$model <- "sigmoidal"

    parameterDf <- as.data.frame(parameterList)
    parameterDf <- sigmoidalRenormalizeParameters_h0(parameterDf, isalist)
  }

  return(parameterDf)
}

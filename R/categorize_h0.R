#' Preliminary Categorization of Signal Presence
#'
#' Conducts the first screening step on normalized time–intensity data to
#' decide whether there is any signal worth fitting, based on intensity
#' range and maximum thresholds.
#'
#' @param normalizedInput A list returned by \code{\link{normalizeData_h0}},
#'   containing at least the fields:
#'   \describe{
#'     \item{\code{dataInputName}}{Original name of the data input.}
#'     \item{\code{dataScalingParameters}}{A list with elements
#'       \code{intensityMax} and \code{intensityRange}.}
#'   }
#' @param threshold_intensity_range Numeric scalar. Minimum fraction of the
#'   full intensity range that must be spanned to consider there is any
#'   dynamic change. Defaults to \code{0.1}.
#' @param threshold_minimum_for_intensity_maximum Numeric scalar. Minimum
#'   fraction of the maximum intensity that must be exceeded to consider
#'   a signal. Defaults to \code{0.3}.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{dataInputName}}{As in \code{normalizedInput}.}
#'     \item{\code{intensityMaximum}}{Observed maximum intensity of the data.}
#'     \item{\code{test.minimum_for_intensity_maximum}}{Logical; \code{TRUE}
#'       if \code{intensityMaximum} exceeds \code{threshold_minimum_for_intensity_maximum}.}
#'     \item{\code{intensityRange}}{Observed intensity range of the data.}
#'     \item{\code{test.intensity_range}}{Logical; \code{TRUE} if
#'       \code{intensityRange} exceeds \code{threshold_intensity_range}.}
#'     \item{\code{decision}}{Character; \code{"not_no_signal"} if both tests
#'       pass, otherwise \code{"no_signal"}.}
#'     \item{\code{decisionSteps}}{String summarizing which sub-steps were
#'       triggered (e.g. \code{"1a_1b_1c"}).}
#'   }
#'
#' @export
preCategorize_h0 <- function (normalizedInput, threshold_intensity_range = 0.1,
                              threshold_minimum_for_intensity_maximum = 0.3)
{
  decisionList <- list()
  decisionList$dataInputName <- normalizedInput$dataInputName
  decisionList$intensityMaximum <- normalizedInput$dataScalingParameters[["intensityMax"]]
  decisionList$threshold_minimum_for_intensity_maximum <- threshold_minimum_for_intensity_maximum
  decisionList$test.minimum_for_intensity_maximum <- threshold_minimum_for_intensity_maximum <
    decisionList$intensityMaximum
  decisionList$intensityRange <- normalizedInput$dataScalingParameters[["intensityRange"]]
  decisionList$threshold_intensity_range <- threshold_intensity_range
  decisionList$test.intensity_range <- threshold_intensity_range <
    decisionList$intensityRange
  choices <- c("no_signal", "sigmoidal", "double_sigmoidal",
               "ambiguous")
  decisonSteps <- c()
  if (!decisionList$test.minimum_for_intensity_maximum) {
    decisonSteps <- c(decisonSteps, "1a")
  }
  if (!decisionList$test.minimum_for_intensity_maximum) {
    decisonSteps <- c(decisonSteps, "1b")
  }
  if (!setequal(choices, c("no_signal"))) {
    decisonSteps <- c(decisonSteps, "1c")
  }
  decisionList$decisonSteps <- paste0(decisonSteps, collapse = "_")
  decisionList$decision <- ifelse(decisionList$test.minimum_for_intensity_maximum &
                                    decisionList$test.intensity_range, "not_no_signal",
                                  "no_signal")
  return(decisionList)
}

#' Decide between Sigmoidal and Double-Sigmoidal Fits
#'
#' Runs a series of logical and AIC‐based tests to choose the best model
#' (“sigmoidal”, “double_sigmoidal”, “no_signal” or “ambiguous”) given two
#' fitted parameter vectors (one sigmoidal, one double-sigmoidal).
#'
#' @param parameterVectorSigmoidal A named list or data frame of sigmoidal fit parameters,
#'   produced by \code{\link{sigFit}} and
#'   \code{\link{parameterCalculation_h0}}.
#' @param parameterVectorDoubleSigmoidal A named list or data.frame of double-sigmoidal fit
#'   parameters, produced by \code{\link{sigFit2}} and
#'   \code{\link{parameterCalculation_h0}}.
#' @param threshold_intensity_range Numeric; minimum required intensity range to consider
#'   any signal (default 0.1).
#' @param threshold_minimum_for_intensity_maximum Numeric; minimum maximum intensity
#'   required to consider a real signal (default 0.3).
#' @param threshold_bonus_sigmoidal_AIC Numeric; AIC bonus given to the sigmoidal model
#'   when comparing with the double-sigmoidal (default 0).
#' @param threshold_sm_tmax_IntensityRatio Numeric; minimum fraction of maximum
#'   intensity reached by sigmoidal at final time to allow sigmoidal (default 0.85).
#' @param threshold_dsm_tmax_IntensityRatio Numeric; maximum fraction of maximum
#'   intensity reached by double-sigmoidal at final time to allow double-sigmoidal (default 0.75).
#' @param threshold_AIC Numeric; maximum AIC value to allow either model (default –10).
#' @param threshold_t0_max_int Numeric; maximum allowed starting intensity ratio
#'   at time zero (default 0.05).
#' @param showDetails Logical; if TRUE, prints the internal decision list structure
#'   for debugging (default FALSE).
#'
#' @return A named list with components:
#'   \item{decisionList}{All intermediate TRUE/FALSE tests and thresholds.}
#'   \item{decision}{One of “no_signal”, “sigmoidal”, “double_sigmoidal” or “ambiguous”.}
#'   \item{decisonSteps}{Which numbered tests were applied (as a single string).}
#'
#' @export
Categorize_h0 <- function (parameterVectorSigmoidal, parameterVectorDoubleSigmoidal,
                           threshold_intensity_range = 0.1, threshold_minimum_for_intensity_maximum = 0.3,
                           threshold_bonus_sigmoidal_AIC = 0, threshold_sm_tmax_IntensityRatio = 0.75,
                           threshold_dsm_tmax_IntensityRatio = 0.75, threshold_AIC = -10,
                           threshold_t0_max_int = 0.05, showDetails = FALSE)
{
  decisionList <- list()
  if ((!is.na(parameterVectorSigmoidal$dataInputName)) & (!is.na(parameterVectorDoubleSigmoidal$dataInputName))) {
    decisionList$test.name <- parameterVectorSigmoidal$dataInputName ==
      parameterVectorDoubleSigmoidal$dataInputName
    decisionList$dataInputName <- as.vector(parameterVectorSigmoidal$dataInputName)
  }
  else {
    decisionList$test.name <- NA
    decisionList$dataInputName <- NA
  }
  decisionList$test.sm_modelCheck <- parameterVectorSigmoidal$model ==
    "sigmoidal"
  decisionList$test.dsm_modelCheck <- parameterVectorDoubleSigmoidal$model ==
    "doublesigmoidal"
  test.timeRange <- parameterVectorSigmoidal$dataScalingParameters.timeRange ==
    parameterVectorDoubleSigmoidal$dataScalingParameters.timeRange
  if (test.timeRange) {
    timeRange <- parameterVectorSigmoidal$dataScalingParameters.timeRange
  }
  test.intensityMin <- parameterVectorSigmoidal$dataScalingParameters.intensityMin ==
    parameterVectorDoubleSigmoidal$dataScalingParameters.intensityMin
  test.intensityMax <- parameterVectorSigmoidal$dataScalingParameters.intensityMax ==
    parameterVectorDoubleSigmoidal$dataScalingParameters.intensityMax
  if (test.intensityMax) {
    intensityMax <- parameterVectorSigmoidal$dataScalingParameters.intensityMax
  }
  test.intensityRange <- parameterVectorSigmoidal$dataScalingParameters.intensityRange ==
    parameterVectorDoubleSigmoidal$dataScalingParameters.intensityRange
  if (test.intensityRange) {
    intensityRange <- parameterVectorSigmoidal$dataScalingParameters.intensityRange
  }
  decisionList$test.dataScalingParameters <- test.timeRange &
    test.intensityMin & test.intensityMax & test.intensityRange
  decisionList$intensityMaximum <- intensityMax
  decisionList$threshold_minimum_for_intensity_maximum <- threshold_minimum_for_intensity_maximum
  decisionList$test.minimum_for_intensity_maximum <- threshold_minimum_for_intensity_maximum <
    intensityMax
  decisionList$intensityRange <- intensityRange
  decisionList$threshold_intensity_range <- threshold_intensity_range
  decisionList$test.intensity_range <- threshold_intensity_range <
    intensityRange
  decisionList$test.sigmoidalFit <- parameterVectorSigmoidal$isThisaFit
  decisionList$test.doublesigmoidalFit <- parameterVectorDoubleSigmoidal$isThisaFit
  decisionList$test.sigmoidalAdditionalParameters <- parameterVectorSigmoidal$additionalParameters
  decisionList$test.doublesigmoidalAdditionalParameters <- parameterVectorDoubleSigmoidal$additionalParameters
  decisionList$threshold_AIC <- threshold_AIC
  decisionList$sigmoidalAIC <- parameterVectorSigmoidal$AIC_value
  decisionList$test.sigmoidalAIC <- parameterVectorSigmoidal$AIC_value <
    threshold_AIC
  decisionList$doublesigmoidalAIC <- parameterVectorDoubleSigmoidal$AIC_value
  decisionList$test.doublesigmoidalAIC <- parameterVectorDoubleSigmoidal$AIC_value <
    threshold_AIC
  decisionList$dsm_maximum_x <- parameterVectorDoubleSigmoidal$maximum_x
  decisionList$timeRange <- timeRange
  decisionList$test.dsm_maximum_x <- decisionList$dsm_maximum_x <
    timeRange
  sm_intensity_at_tmax <- sigmoidalFitFormula_h0(x = timeRange,
                                                 maximum = parameterVectorSigmoidal$maximum_y, slopeParam = parameterVectorSigmoidal$slopeParam_Estimate,
                                                 midPoint = parameterVectorSigmoidal$midPoint_Estimate, h0 = parameterVectorSigmoidal$h0_Estimate) # changed but dont know if it works yet
  decisionList$sm_tmax_IntensityRatio <- sm_intensity_at_tmax/parameterVectorSigmoidal$maximum_y
  decisionList$threshold_sm_tmax_IntensityRatio <- threshold_sm_tmax_IntensityRatio
  decisionList$test.sm_tmax_IntensityRatio <- decisionList$sm_tmax_IntensityRatio >
    threshold_sm_tmax_IntensityRatio
  dsm_intensity_at_tmax <- doubleSigmoidalFitFormula_h0(x = timeRange,
                                                        finalAsymptoteIntensityRatio = parameterVectorDoubleSigmoidal$finalAsymptoteIntensityRatio_Estimate,
                                                        maximum = parameterVectorDoubleSigmoidal$maximum_y,
                                                        slope1Param = parameterVectorDoubleSigmoidal$slope1Param_Estimate,
                                                        midPoint1Param = parameterVectorDoubleSigmoidal$midPoint1Param_Estimate,
                                                        slope2Param = parameterVectorDoubleSigmoidal$slope2Param_Estimate,
                                                        midPointDistanceParam = parameterVectorDoubleSigmoidal$midPointDistanceParam_Estimate,
                                                        h0 = parameterVectorDoubleSigmoidal$h0_Estimate) # changed but don't know if it works yet
  decisionList$dsm_tmax_IntensityRatio <- dsm_intensity_at_tmax/parameterVectorDoubleSigmoidal$maximum_y
  decisionList$threshold_dsm_tmax_IntensityRatio <- threshold_dsm_tmax_IntensityRatio
  decisionList$test.dsm_tmax_IntensityRatio <- decisionList$dsm_tmax_IntensityRatio <
    threshold_dsm_tmax_IntensityRatio
  decisionList$sm_startPoint_x <- parameterVectorSigmoidal$startPoint_x
  decisionList$test.sm_startPoint_x <- parameterVectorSigmoidal$startPoint_x >
    0
  decisionList$dsm_startPoint_x <- parameterVectorDoubleSigmoidal$startPoint_x
  decisionList$test.dsm_startPoint_x <- parameterVectorDoubleSigmoidal$startPoint_x >
    0
  decisionList$sm_startIntensity <- sigmoidalFitFormula_h0(x = 0,
                                                           maximum = parameterVectorSigmoidal$maximum_y, slopeParam = parameterVectorSigmoidal$slopeParam_Estimate,
                                                           midPoint = parameterVectorSigmoidal$midPoint_Estimate, h0 = parameterVectorSigmoidal$h0_Estimate)    # changed but don't know if it works yet
  decisionList$threshold_t0_max_int <- threshold_t0_max_int
  decisionList$test.sm_startIntensity <- decisionList$sm_startIntensity <
    threshold_t0_max_int
  decisionList$dsm_startIntensity <- doubleSigmoidalFitFormula_h0(x = 0,
                                                                  finalAsymptoteIntensityRatio = parameterVectorDoubleSigmoidal$finalAsymptoteIntensityRatio_Estimate,
                                                                  maximum = parameterVectorDoubleSigmoidal$maximum_y,
                                                                  slope1Param = parameterVectorDoubleSigmoidal$slope1Param_Estimate,
                                                                  midPoint1Param = parameterVectorDoubleSigmoidal$midPoint1Param_Estimate,
                                                                  slope2Param = parameterVectorDoubleSigmoidal$slope2Param_Estimate,
                                                                  midPointDistanceParam = parameterVectorDoubleSigmoidal$midPointDistanceParam_Estimate,
                                                                  h0 = parameterVectorDoubleSigmoidal$h0_Estimate) # changed but dont know if it works yet
  decisionList$test.dsm_startIntensity <- decisionList$dsm_startIntensity <
    threshold_t0_max_int
  if (!is.na(decisionList$test.name)) {
    if (!decisionList$test.name) {
      stop("dataNames of sigmoidal and double-sigmoidal fits must be same")
    }
  }
  if (!decisionList$test.sm_modelCheck) {
    stop("provided sigmoidal model must be fitted by sigFit()")
  }
  if (!decisionList$test.dsm_modelCheck) {
    stop("provided double_sigmoidal model must be fitted by sigFit2()")
  }
  if (!decisionList$test.dataScalingParameters) {
    stop("data scaling parameters of provided sigmoidal and double_sigmoidal parameterVectors must be same")
  }
  if (!decisionList$test.sigmoidalAdditionalParameters) {
    stop("additional parameters for sigmoidal fit must be calculated")
  }
  if (!decisionList$test.doublesigmoidalAdditionalParameters) {
    stop("additional parameters for double_sigmoidal fit must be calculated")
  }
  choices <- c("no_signal", "sigmoidal", "double_sigmoidal",
               "ambiguous")
  decisonSteps <- c()
  if (!decisionList$test.minimum_for_intensity_maximum) {
    decisonSteps <- c(decisonSteps, "1a")
    choices <- setdiff(choices, c("sigmoidal", "double_sigmoidal",
                                  "ambiguous"))
  }
  if (!decisionList$test.minimum_for_intensity_maximum) {
    decisonSteps <- c(decisonSteps, "1b")
    choices <- setdiff(choices, c("sigmoidal", "double_sigmoidal",
                                  "ambiguous"))
  }
  if (!setequal(choices, c("no_signal"))) {
    decisonSteps <- c(decisonSteps, "1c")
    choices <- setdiff(choices, c("no_signal"))
  }
  if (!decisionList$test.sigmoidalFit) {
    decisonSteps <- c(decisonSteps, "2a")
    choices <- setdiff(choices, c("sigmoidal"))
  }
  if (!decisionList$test.sigmoidalFit) {
    decisonSteps <- c(decisonSteps, "2b")
    choices <- setdiff(choices, c("double_sigmoidal"))
  }
  if (!decisionList$test.sigmoidalAIC) {
    decisonSteps <- c(decisonSteps, "3a")
    choices <- setdiff(choices, c("sigmoidal"))
  }
  if (!decisionList$test.doublesigmoidalAIC) {
    decisonSteps <- c(decisonSteps, "3b")
    choices <- setdiff(choices, c("double_sigmoidal"))
  }
  if (!decisionList$test.sm_startPoint_x) {
    decisonSteps <- c(decisonSteps, "4a")
    choices <- setdiff(choices, c("sigmoidal"))
  }
  if (!decisionList$test.dsm_startPoint_x) {
    decisonSteps <- c(decisonSteps, "4b")
    choices <- setdiff(choices, c("double_sigmoidal"))
  }
  #if (!decisionList$test.sm_startIntensity) {
  #  decisonSteps <- c(decisonSteps, "5a")
  #  choices <- setdiff(choices, c("sigmoidal"))
  #}
  #if (!decisionList$test.dsm_startIntensity) {
  #  decisonSteps <- c(decisonSteps, "5b")
  #  choices <- setdiff(choices, c("double_sigmoidal"))
  #}
  if (!decisionList$test.dsm_tmax_IntensityRatio) {
    decisonSteps <- c(decisonSteps, "6")
    choices <- setdiff(choices, c("double_sigmoidal"))
  }
  if (!decisionList$test.sm_tmax_IntensityRatio) {
    decisonSteps <- c(decisonSteps, "7")
    choices <- setdiff(choices, c("sigmoidal"))
  }
  if (length(intersect(choices, c("sigmoidal", "double_sigmoidal"))) >
      0) {
    decisonSteps <- c(decisonSteps, "8")
    choices <- setdiff(choices, c("ambiguous"))
  }
  if (length(intersect(choices, c("sigmoidal", "double_sigmoidal"))) ==
      2) {
    decisonSteps <- c(decisonSteps, "9")
    if (decisionList$sigmoidalAIC + threshold_bonus_sigmoidal_AIC <
        decisionList$doublesigmoidalAIC) {
      choices <- setdiff(choices, c("double_sigmoidal"))
    }
    if (decisionList$sigmoidalAIC + threshold_bonus_sigmoidal_AIC >
        decisionList$doublesigmoidalAIC) {
      choices <- setdiff(choices, c("sigmoidal"))
    }
  }
  if (!length(choices) == 1) {
    stop("At this point length of choice must be 1")
  }
  decisionList$decisonSteps <- paste0(decisonSteps, collapse = "_")
  decisionList$decision <- paste0(choices, collapse = "_")
  if (showDetails) {
    utils::str(decisionList)
  }
  return(decisionList)
}

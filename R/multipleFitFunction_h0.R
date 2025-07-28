#' Run multiple fitting attempts and return the best sigmoidal or double‐sigmoidal fit
#'
#' Performs repeated fits of either the sigmoidal or double‐sigmoidal model to the
#' same dataset until a minimum number of successful fits is reached or a maximum
#' number of attempts is exhausted.  Returns the parameter vector from the best
#' fit (lowest residual sum of squares), together with counts of total and successful runs.
#'
#' @param dataInput A data frame (with columns `time` and `intensity`) or a list
#'   containing `timeIntensityData` and optional `dataScalingParameters`.  See
#'   \code{\link{dataCheck_h0}} for details.
#' @param dataInputName Optional name to assign to the input dataset; used for
#'   tracking/identification when fitting multiple curves.  Defaults to NA.
#' @param model Character; either `"sigmoidal"` or `"doublesigmoidal"`, indicating
#'   which underlying model to fit.
#' @param n_runs_min Integer; the minimum number of *successful* fits to perform
#'   before stopping.  Defaults to 20.
#' @param n_runs_max Integer; the maximum total number of fitting attempts allowed.
#'   Defaults to 500.
#' @param showDetails Logical; if TRUE, prints progress counters (best residual,
#'   number of successful fits, etc.) after each successful fit.  Defaults to FALSE.
#' @param ... Additional arguments passed to the underlying fitting function:
#'   \code{\link{sigFit}} when `model = "sigmoidal"`, or
#'   \code{\link{sigFit2}} when `model = "doublesigmoidal"`.
#'
#' @return A named list (invisibly) with components:
#'   * All fit parameters and goodness‐of‐fit metrics from the best run,
#'   * \code{betterFit}: count of times a new best residual was found,
#'   * \code{correctFit}: count of successful fits,
#'   * \code{totalFit}: total number of attempts made.
#'
#' @examples
#' \dontrun{
#' # single‐sigmoidal over noisy data, require at least 10 good fits
#' res <- multipleFitFunction_h0(
#'   dataInput    = my_data,
#'   model        = "sigmoidal",
#'   n_runs_min   = 10,
#'   n_runs_max   = 200
#' )
#' plot(res$time, res$intensity, type = "l")
#' }
#'
#' @export
multipleFitFunction_h0 <- function (dataInput, dataInputName = NA, model, n_runs_min = 20,
                                    n_runs_max = 500, showDetails = FALSE, ...)
{
  dataInputCheck <- dataCheck_h0(dataInput)
  if (!(model %in% c("sigmoidal", "doublesigmoidal"))) {
    stop("model should be one of sigmoidal, doublesigmoidal")
  }
  counterBetterFit <- 0
  counterCorrectFit <- 0
  counterTotalFit <- 0
  residual_Sum_of_Squares_min <- Inf
  storedModelOutput <- c()
  storedModelOutput$residual_Sum_of_Squares <- Inf
  while (counterCorrectFit < n_runs_min & counterTotalFit <
         n_runs_max) {
    counterTotalFit <- counterTotalFit + 1
    if (model == "sigmoidal") {
      modelOutput <- sigFit(dataInput = dataInput,
                            tryCounter = counterTotalFit, ...)
    }
    if (model == "doublesigmoidal") {
      modelOutput <- sigFit2(dataInput = dataInput,
                             tryCounter = counterTotalFit, ...)
    }
    if (is.na(dataInputName)) {
      isalist <- (is.list(dataInput) & !is.data.frame(dataInput))
      if (isalist) {
        modelOutput$dataInputName <- as.character(dataInput$dataInputName)
      }
      if (!isalist) {
        modelOutput$dataInputName <- NA
      }
    }
    if (!is.na(dataInputName)) {
      isalist <- (is.list(dataInput) & !is.data.frame(dataInput))
      if (isalist) {
        if (is.na(dataInput$dataInputName)) {
          modelOutput$dataInputName <- as.character(dataInputName)
        }
        if (!is.na(dataInput$dataInputName)) {
          if (dataInput$dataInputName != dataInputName) {
            stop("the input data has already have a name")
          }
          if (dataInput$dataInputName == dataInputName) {
            modelOutput$dataInputName <- as.character(dataInputName)
          }
        }
      }
      if (!isalist) {
        modelOutput$dataInputName <- as.character(dataInputName)
      }
    }
    if (modelOutput[["isThisaFit"]]) {
      counterCorrectFit <- counterCorrectFit + 1
      if (residual_Sum_of_Squares_min > modelOutput$residual_Sum_of_Squares) {
        counterBetterFit <- counterBetterFit + 1
        residual_Sum_of_Squares_min <- modelOutput$residual_Sum_of_Squares
        storedModelOutput <- modelOutput
      }
    }
    if (showDetails) {
      print(c(betterFit = counterBetterFit, correctFit = counterCorrectFit,
              totalFit = counterTotalFit, currentOutput = modelOutput$residual_Sum_of_Squares,
              bestOutput = storedModelOutput$residual_Sum_of_Squares))
    }
  }
  storedModelOutput$betterFit <- counterBetterFit
  storedModelOutput$correctFit <- counterCorrectFit
  storedModelOutput$totalFit <- counterTotalFit
  return(storedModelOutput)
}

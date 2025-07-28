#' Normalize Time–Intensity Data to \code{[0,1]} Range
#'
#' Scales a two-column data frame of raw time and intensity values so that
#' time runs from 0 to 1 and intensity runs from 0 to 1.  Returns both the
#' normalized data and the scaling parameters needed to reverse the transform.
#'
#' @param dataInput A data.frame with numeric columns \code{time} and
#'   \code{intensity}.
#' @param dataInputName Optional character string giving the name or ID of
#'   the input dataset.
#' @return A list with components:
#'   \item{timeIntensityData}{data.frame of normalized \code{time} and
#'     \code{intensity} in \code{[0,1]}.}
#'   \item{dataScalingParameters}{numeric vector with
#'     \code{timeRange}, \code{intensityMin}, \code{intensityMax}, and
#'     \code{intensityRange} used for normalization.}
#'   \item{dataInputName}{the original \code{dataInputName} value.}
#' @examples
#' df <- data.frame(time = 0:10, intensity = rnorm(11, 5, 2))
#' normalizeData_h0(df, dataInputName = "example")
#' @export
normalizeData_h0 <- function (dataInput, dataInputName = NA)
{
  dataInputCheckVariable <- dataCheck_h0(dataInput)
  timeData <- dataInput$time
  timeRange <- max(timeData, na.rm = T)
  timeData <- timeData/timeRange
  intensityMin <- min(dataInput$intensity, na.rm = T)
  intensityMax <- max(dataInput$intensity, na.rm = T)
  intensityData <- dataInput$intensity - intensityMin
  intensityRange <- max(intensityData, na.rm = T)
  intensityData <- intensityData/intensityRange
  dataOutput <- data.frame(time = timeData, intensity = intensityData)
  return(list(timeIntensityData = dataOutput, dataScalingParameters = c(timeRange = timeRange,
                                                                        intensityMin = intensityMin, intensityMax = intensityMax,
                                                                        intensityRange = intensityRange), dataInputName = dataInputName))
}


#' Unnormalize Time–Intensity Data
#'
#' Reverses the transformation applied by \code{\link{normalizeData_h0}},
#' restoring raw time and intensity values from scaled data.
#'
#' @param dataInput A list as returned by \code{\link{normalizeData_h0}},
#'   containing \code{timeIntensityData} and \code{dataScalingParameters}.
#' @return A list with component:
#'   \item{timeIntensityData}{data.frame of original \code{time} and
#'     \code{intensity} values.}
#' @examples
#' norm <- normalizeData_h0(data.frame(time=0:10, intensity=rnorm(11)))
#' unnormalizeData_h0(norm)
#' @export
unnormalizeData_h0 <- function (dataInput)
{
  dataInputCheckVariable <- dataCheck_h0(dataInput)
  time <- dataInput$dataScalingParameters[["timeRange"]] *
    dataInput$timeIntensityData[["time"]]
  intensity <- dataInput$dataScalingParameters[["intensityMin"]] +
    dataInput$dataScalingParameters[["intensityRange"]] *
    dataInput$timeIntensityData[["intensity"]]
  dataOutput <- list(timeIntensityData = data.frame(time,
                                                    intensity))
  return(dataOutput)
}

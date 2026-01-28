#' Calculate saturation vapour pressure using the Buck equation
#'
#' @param temp Temperature in degrees Celsius.
#' @return Saturation vapour pressure in kPa.
#' @keywords internal
calc_sat_vp <- function(temp) {
  0.61121 * exp((18.678 - temp / 234.5) * (temp / (257.14 + temp)))
}

#' Calculate dew point temperature
#'
#' Uses the inverse of the Buck equation to derive dew point from
#' temperature and relative humidity.
#'
#' @param temp Air temperature in degrees Celsius.
#' @param RH Relative humidity in percent.
#' @return Dew point temperature in degrees Celsius.
#' @keywords internal
calc_dew_point <- function(temp, RH) {
  es <- calc_sat_vp(temp)
  e <- (RH / 100) * es
  # Inverse of Buck equation
  log_term <- log(e / 0.61121)
  (257.14 * log_term) / (18.678 - log_term)
}

#' Calculate psychrometric wet bulb temperature
#'
#' Uses the Stull (2011) approximation for aspirated (psychrometric)
#' wet bulb temperature.
#'
#' @param temp Air temperature in degrees Celsius.
#' @param RH Relative humidity in percent.
#' @return Psychrometric (aspirated) wet bulb temperature in degrees Celsius.
#'
#' @references
#' Stull, R. (2011) Wet-Bulb Temperature from Relative Humidity and Air
#' Temperature. Journal of Applied Meteorology and Climatology, 50(11),
#' 2267-2269.
#' @keywords internal
calc_aspirated_wb <- function(temp, RH) {
  temp * atan(0.151977 * sqrt(RH + 8.313659)) +
    atan(temp + RH) - atan(RH - 1.676331) +
    0.00391838 * RH^1.5 * atan(0.023101 * RH) - 4.686035
}

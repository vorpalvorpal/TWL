#' Physical and physiological constants for TWL calculation
#'
#' A list of constants used throughout the TWL calculations:
#' \describe{
#'   \item{STEFAN_BOLTZMANN}{Stefan-Boltzmann constant, 5.67e-8 W/(m^2.K^4)}
#'   \item{ABS_ZERO}{Absolute zero offset, 273.15 K}
#'   \item{LATENT_HEAT_EVAP}{Latent heat of evaporation at skin temp, 2430 kJ/kg}
#'   \item{VIEW_EMISSIVITY_FACTOR}{Combined view (0.8) and emissivity (0.95) factor, 0.76}
#'   \item{BULB_DIA}{Standard wick diameter, 0.004 m}
#'   \item{AIR_THERMAL_CONDUCTIVITY}{Thermal conductivity of air, 0.028 W/(m.K)}
#'   \item{GLOBE_DIAMETER}{Black globe diameter, 0.15 m}
#'   \item{GLOBE_EMISSIVITY}{Globe emissivity, 0.95}
#' }
#' @keywords internal
TWL_CONSTANTS <- list(
  STEFAN_BOLTZMANN = 5.67e-8,
  ABS_ZERO = 273.15,
  LATENT_HEAT_EVAP = 2430,
  VIEW_EMISSIVITY_FACTOR = 0.76,
  BULB_DIA = 0.004,
  AIR_THERMAL_CONDUCTIVITY = 0.028,
  GLOBE_DIAMETER = 0.15,
  GLOBE_EMISSIVITY = 0.95
)

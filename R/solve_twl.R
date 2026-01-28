#' Solve TWL for a single observation
#'
#' Iteratively solves the heat balance equations to find the maximum
#' sustainable metabolic rate using a three-zone evaporation model and
#' Newton's method.
#'
#' @param temp Air temperature in degrees Celsius.
#' @param wind_speed Wind speed in m/s.
#' @param RH Relative humidity in percent.
#' @param direct_solar Direct solar radiation in W/m^2.
#' @param diffuse_solar Diffuse solar radiation in W/m^2.
#' @param pressure Atmospheric pressure in kPa.
#' @param solar_radiation Total solar radiation in W/m^2.
#' @param pa Actual vapour pressure in kPa.
#' @param temp_dewpoint Dew point temperature in degrees Celsius.
#' @param wet_bulb Natural wet bulb temperature in degrees Celsius.
#' @param globe_temp Globe temperature in degrees Celsius.
#' @param trad Mean radiant temperature in degrees Celsius.
#' @param hc Convective heat transfer coefficient in W/(m^2.K).
#' @param max_core_temp Maximum core temperature in degrees Celsius.
#' @param max_sweat_rate Maximum sweat rate in kg/(m^2.hr).
#' @param Icl Clothing insulation in clo.
#' @param icl Clothing vapour permeability (0--1).
#' @param LR Lewis relation in K/kPa.
#' @param lambda Latent heat of evaporation in kJ/kg.
#' @param fr Posture factor.
#' @param index Observation index for warnings (optional).
#'
#' @return TWL in W/m^2.
#' @keywords internal
solve_twl_single <- function(temp, wind_speed, RH, direct_solar, diffuse_solar,
                             pressure, solar_radiation, pa, temp_dewpoint,
                             wet_bulb, globe_temp, trad, hc,
                             max_core_temp, max_sweat_rate,
                             Icl, icl, LR, lambda, fr,
                             index = NULL) {

  # Constants
  ABS_ZERO <- 273.15
  STEFAN_BOLTZMANN <- 5.67e-8

  # Clothing parameters
  Rcl <- Icl * 0.155
  fcl <- 1 + 0.31 * Icl
  Recl <- Rcl / (icl * LR)

  # Radiative heat transfer coefficient with posture factor
  hr <- 4.61 * (1 + (trad + 35) / 546)^3 * fr

  # Combined heat transfer coefficient
  h <- hr + hc

  # Operative temperature
  toper <- (hr * trad + hc * temp) / h

  # Clothing corrections
  Fcle <- fcl / (1 + fcl * h * Rcl)
  he <- LR * hc
  Fpcl <- 1 / (1 + fcl * he * Recl)

  # Maximum evaporative capacity at full skin wettedness
  # Using skin temp of 35 degrees C as initial estimate
  ps_35 <- calc_sat_vp(35)
  E_max <- Fpcl * fcl * he * (ps_35 - pa)

  # Physiological conductance from Wyndham's data
  K_cs_base <- 14.3

  # Iterative solution for TWL
  t_skin <- 35
  tolerance <- 0.5
  max_iterations <- 100

  best_M <- NA
  best_error <- Inf

  for (iter in 1:max_iterations) {
    # Saturated vapour pressure at skin
    ps <- calc_sat_vp(t_skin)

    # Maximum evaporative heat loss
    E_max_skin <- Fpcl * fcl * he * (ps - pa)
    E_max_skin <- max(0, E_max_skin)

    # Limit by maximum sweat rate
    max_E_sweat <- max_sweat_rate * lambda / 3.6

    # Three-zone evaporation model based on skin wettedness
    # Zone 1: w < 0.4 (efficient)
    # Zone 2: 0.4 <= w < 1.0 (reduced efficiency)
    # Zone 3: w = 1.0 (dripping, constant)
    if (E_max_skin <= 0) {
      E <- 0
    } else if (E_max_skin <= 0.4 * max_E_sweat) {
      E <- E_max_skin
    } else if (E_max_skin <= max_E_sweat) {
      w <- E_max_skin / max_E_sweat
      efficiency <- 1 - 0.5 * (w - 0.4) / 0.6
      E <- E_max_skin * efficiency
    } else {
      E <- max_E_sweat * 0.85
    }

    # Sensible heat loss (convection + radiation)
    CR <- Fcle * h * (t_skin - toper)

    # Heat transfer from core to skin
    K_cs <- K_cs_base * (1 + 0.001 * E)
    H <- K_cs * (max_core_temp - t_skin)

    # Heat balance: H should equal CR + E
    heat_balance_error <- H - (CR + E)

    # Track best solution
    if (abs(heat_balance_error) < best_error) {
      best_error <- abs(heat_balance_error)

      # Calculate metabolic rate from core-skin heat transfer
      # M = H + B (respiratory losses)
      # B = 0.0014*M*(34-temp) + 0.0173*M*(5.87-pa)
      # Solving algebraically: M = H / (1 - resp_coeff)
      respiratory_coeff <- 0.0014 * (34 - temp) + 0.0173 * (5.87 - pa)

      if (respiratory_coeff < 0.95 && respiratory_coeff >= 0) {
        best_M <- H / (1 - respiratory_coeff)
      } else {
        best_M <- H * 1.05
      }
    }

    # Check convergence
    if (abs(heat_balance_error) < tolerance) {
      break
    }

    # Adjust skin temperature
    adjustment <- -heat_balance_error / (K_cs + Fcle * h + 50)
    adjustment <- max(-0.5, min(0.5, adjustment))
    t_skin <- t_skin + adjustment

    # Keep skin temperature in valid range
    t_skin_min <- temp_dewpoint + 1
    t_skin_max <- max_core_temp - 0.5
    t_skin <- max(t_skin_min, min(t_skin_max, t_skin))
  }

  # Return best metabolic rate found
  if (!is.na(best_M) && best_M > 0) {
    return(best_M)
  }

  # Fallback calculation if iteration failed
  ps_35 <- calc_sat_vp(35)
  E_max_est <- Fpcl * fcl * he * (ps_35 - pa)
  E_max_est <- max(0, min(max_E_sweat, E_max_est))
  CR_est <- Fcle * h * (35 - toper)
  H_est <- E_max_est + CR_est

  respiratory_coeff <- 0.0014 * (34 - temp) + 0.0173 * (5.87 - pa)
  if (respiratory_coeff < 0.95 && respiratory_coeff >= 0) {
    M_est <- H_est / (1 - respiratory_coeff)
  } else {
    M_est <- H_est * 1.05
  }

  return(max(60, M_est))
}

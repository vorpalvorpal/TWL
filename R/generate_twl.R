#' Calculate Thermal Work Limit (TWL)
#'
#' Calculates the maximum sustainable metabolic rate (W/m^2) that
#' well-hydrated, acclimatised individuals can maintain in a specific
#' thermal environment. Based on the Brake & Bates (2002) methodology.
#'
#' Weather parameters (`temp`, `wind_speed`, `RH`, `direct_solar`,
#' `diffuse_solar`, `cloud_cover`, `pressure`) that are not supplied are
#' automatically retrieved from the
#' [Open-Meteo API](https://open-meteo.com/) using the provided `datetime`,
#' `latitude` and `longitude`. An internet connection is required for
#' auto-fetching.
#'
#' @param datetime POSIXct datetime vector (required).
#' @param latitude Latitude in decimal degrees (required).
#' @param longitude Longitude in decimal degrees (required).
#' @param temp Dry bulb air temperature in degrees Celsius, or `NULL` to
#'   fetch from Open-Meteo.
#' @param wind_speed Wind speed in m/s, or `NULL` to fetch from Open-Meteo.
#' @param RH Relative humidity in percent, or `NULL` to fetch from
#'   Open-Meteo.
#' @param direct_solar Direct beam solar radiation in W/m^2, or `NULL` to
#'   fetch from Open-Meteo.
#' @param diffuse_solar Diffuse sky solar radiation in W/m^2, or `NULL` to
#'   fetch from Open-Meteo.
#' @param cloud_cover Cloud cover fraction (0--1), or `NULL` to fetch from
#'   Open-Meteo.
#' @param pressure Barometric pressure in hPa (or kPa if `convert_pressure =
#'   FALSE`), or `NULL` to fetch from Open-Meteo.
#' @param albedo Ground albedo (0--1). Default 0.12 for asphalt.
#' @param Icl Intrinsic clothing thermal resistance in clo. Default 0.6.
#' @param icl Clothing vapour permeation efficiency (0--1). Default 0.45.
#' @param max_core_temp Maximum acceptable core temperature in degrees Celsius.
#'   Default 38.2.
#' @param max_sweat_rate Maximum acceptable sweat rate in kg/(m^2.hr).
#'   Default 0.67.
#' @param convert_pressure If `TRUE` (the default), converts pressure from hPa
#'   to kPa.
#' @param verbose If `TRUE` (the default), show progress and diagnostic
#'   messages.
#'
#' @return Numeric vector of TWL values in W/m^2.
#'
#' @details
#' TWL categories (from Brake & Bates 2002):
#' \itemize{
#'   \item > 220 W/m^2: Unrestricted work
#'   \item 140--220 W/m^2: Acclimatisation zone
#'   \item 115--140 W/m^2: Buffer zone
#'   \item < 115 W/m^2: Withdrawal required
#' }
#'
#' @references
#' Brake, D.J. and Bates, G.P. (2002) Limiting Metabolic Rate (Thermal Work
#' Limit) as an Index of Thermal Stress. Applied Occupational and Environmental
#' Hygiene, 17:3, 176-186.
#'
#' @export
generate_twl <- function(datetime,
                         latitude,
                         longitude,
                         temp = NULL,
                         wind_speed = NULL,
                         RH = NULL,
                         direct_solar = NULL,
                         diffuse_solar = NULL,
                         cloud_cover = NULL,
                         pressure = NULL,
                         albedo = 0.12,
                         Icl = 0.6,
                         icl = 0.45,
                         max_core_temp = 38.2,
                         max_sweat_rate = 0.67,
                         convert_pressure = TRUE,
                         verbose = TRUE) {

  # --- Determine which weather fields need fetching ---
  field_map <- c(
    temp          = "temperature_2m",
    wind_speed    = "wind_speed_10m",
    RH            = "relative_humidity_2m",
    direct_solar  = "direct_radiation",
    diffuse_solar = "diffuse_radiation",
    cloud_cover   = "cloud_cover",
    pressure      = "surface_pressure"
  )

  missing <- vapply(
    list(temp, wind_speed, RH, direct_solar, diffuse_solar,
         cloud_cover, pressure),
    is.null, logical(1)
  )
  names(missing) <- names(field_map)
  fields_needed <- field_map[missing]

  if (length(fields_needed) > 0L) {
    if (verbose) {
      cli_alert_info(
        "Fetching missing weather data: {paste(names(fields_needed), collapse = ', ')}"
      )
    }
    api_data <- fetch_openmeteo(
      datetime, latitude, longitude, unname(fields_needed), verbose = verbose
    )
    if (is.null(temp))          temp          <- api_data[["temperature_2m"]]
    if (is.null(wind_speed))    wind_speed    <- api_data[["wind_speed_10m"]]
    if (is.null(RH))            RH            <- api_data[["relative_humidity_2m"]]
    if (is.null(direct_solar))  direct_solar  <- api_data[["direct_radiation"]]
    if (is.null(diffuse_solar)) diffuse_solar <- api_data[["diffuse_radiation"]]
    if (is.null(cloud_cover))   cloud_cover   <- api_data[["cloud_cover"]]
    if (is.null(pressure))      pressure      <- api_data[["surface_pressure"]]
  }

  n_obs <- length(temp)

  if (verbose) {
    cli_h1("TWL Calculation")
    cli_alert_info("Processing {n_obs} observation{?s}")
  }

  # Constants
  LR <- 16.5
  lambda <- 2430
  fr <- 0.72

  # Unit conversions
  if (convert_pressure) {
    pressure <- pressure / 10
  }

  # Constrain wind speed to valid range (Brake & Bates recommend 0.2-4.0 m/s)
  wind_speed_orig <- wind_speed
  wind_speed <- pmax(0.2, pmin(4.0, wind_speed))
  n_capped <- sum(wind_speed_orig != wind_speed, na.rm = TRUE)
  if (verbose && n_capped > 0) {
    cli_alert_warning(
      "{n_capped} wind speed value{?s} constrained to [0.2, 4.0] m/s range"
    )
  }

  # Convert cloud cover to fraction if needed (assume % if max > 1)
  if (max(cloud_cover, na.rm = TRUE) > 1) {
    cloud_cover <- cloud_cover / 100
  }

  # Calculate solar position
  if (verbose) cli_alert("Calculating solar position...")
  solar_pos <- calculate_solar_position(datetime, latitude, longitude)

  # Total solar radiation
  solar_radiation <- direct_solar + diffuse_solar

  # Calculate globe temperature
  if (verbose) cli_alert("Calculating globe temperature...")
  globe_temp <- calculate_globe_temp(
    temp, wind_speed, direct_solar, diffuse_solar,
    solar_pos$zenith, albedo, verbose = verbose
  )

  # Mean radiant temperature (approximation from globe temp)
  trad <- globe_temp + 1.5

  # Psychrometric calculations
  if (verbose) cli_alert("Calculating psychrometric variables...")
  es <- calc_sat_vp(temp)
  pa <- (RH / 100) * es
  temp_dewpoint <- calc_dew_point(temp, RH)

  # Convective heat transfer coefficient
  hc <- 8.3 * sqrt(wind_speed)

  # Natural wet bulb temperature
  if (verbose) cli_alert("Calculating natural wet bulb temperature...")
  wet_bulb <- calculate_natural_wet_bulb(
    temp, RH, pressure, wind_speed, globe_temp,
    verbose = verbose, show_progress = FALSE
  )

  # Solve for TWL
  if (verbose) cli_alert("Computing TWL...")

  pb_id <- NULL
  if (verbose && n_obs > 100) {
    pb_id <- cli_progress_bar("Computing TWL", total = n_obs)
  }

  TWL <- pmap_dbl(
    list(
      temp, wind_speed, RH, direct_solar, diffuse_solar,
      pressure, solar_radiation, pa, temp_dewpoint,
      wet_bulb, globe_temp, trad, hc, seq_len(n_obs)
    ),
    function(temp, wind_speed, RH, direct_solar, diffuse_solar,
             pressure, solar_radiation, pa, temp_dewpoint,
             wet_bulb, globe_temp, trad, hc, idx) {

      if (verbose && n_obs > 100 && idx %% 10 == 0) {
        cli_progress_update(id = pb_id)
      }

      # Return NA if inputs are NA
      if (is.na(temp) || is.na(wind_speed) || is.na(RH) ||
          is.na(direct_solar) || is.na(diffuse_solar) || is.na(pressure)) {
        return(NA_real_)
      }

      # Solve for this point
      result <- solve_twl_single(
        temp, wind_speed, RH, direct_solar, diffuse_solar,
        pressure, solar_radiation, pa, temp_dewpoint,
        wet_bulb, globe_temp, trad, hc,
        max_core_temp, max_sweat_rate,
        Icl, icl, LR, lambda, fr,
        index = if (verbose) idx else NULL
      )

      # Apply withdrawal limits
      if (!is.na(result)) {
        # DB > 44 degrees C withdrawal limit
        if (temp > 44) {
          result <- min(result, 115)
        }

        # WB > 32 degrees C withdrawal limit
        if (!is.na(wet_bulb) && wet_bulb > 32) {
          result <- min(result, 115)
        }

        # Constrain to valid range (60-380 W/m^2)
        if (result < 60 || result > 380) {
          if (verbose && (result < 50 || result > 400)) {
            cli_alert_warning(
              "Observation {idx}: TWL ({round(result, 1)} W/m\\u00b2) outside valid range, constrained to [60, 380]"
            )
          }
          result <- max(60, min(380, result))
        }
      }

      result
    }
  )

  if (verbose && n_obs > 100) {
    cli_progress_done(id = pb_id)
  }

  # Summary statistics
  if (verbose) {
    n_valid <- sum(!is.na(TWL))
    n_invalid <- sum(is.na(TWL))

    if (n_valid > 0) {
      twl_range <- range(TWL, na.rm = TRUE)
      twl_mean <- mean(TWL, na.rm = TRUE)

      cli_alert_success(
        "TWL calculation complete: {n_valid} valid result{?s}, {n_invalid} invalid"
      )
      cli_alert_info(
        "TWL range: [{round(twl_range[1], 1)}, {round(twl_range[2], 1)}] W/m\\u00b2, mean: {round(twl_mean, 1)} W/m\\u00b2"
      )

      # Categorise results
      n_withdrawal <- sum(TWL < 115, na.rm = TRUE)
      n_buffer <- sum(TWL >= 115 & TWL < 140, na.rm = TRUE)
      n_acclimatisation <- sum(TWL >= 140 & TWL < 220, na.rm = TRUE)
      n_unrestricted <- sum(TWL >= 220, na.rm = TRUE)

      cli_h2("TWL Categories")
      cli_ul(c(
        "Withdrawal (<115 W/m\\u00b2): {n_withdrawal} ({round(100*n_withdrawal/n_valid, 1)}%)",
        "Buffer (115-140 W/m\\u00b2): {n_buffer} ({round(100*n_buffer/n_valid, 1)}%)",
        "Acclimatisation (140-220 W/m\\u00b2): {n_acclimatisation} ({round(100*n_acclimatisation/n_valid, 1)}%)",
        "Unrestricted (>=220 W/m\\u00b2): {n_unrestricted} ({round(100*n_unrestricted/n_valid, 1)}%)"
      ))

      if (n_withdrawal > 0) {
        cli_alert_warning(
          "{n_withdrawal} observation{?s} in withdrawal zone"
        )
      }
    } else {
      cli_alert_danger("No valid TWL results calculated")
    }
  }

  return(TWL)
}

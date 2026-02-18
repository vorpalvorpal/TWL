#' Calculate solar position (zenith and azimuth angles)
#'
#' Calculates solar zenith and azimuth angles from date, time and location
#' using astronomical equations based on the Julian Day (low-precision
#' algorithm accurate to within ~0.01° for dates within a few centuries of
#' J2000.0).
#'
#' Currently only `zenith` is consumed by [calculate_globe_temp()]; `azimuth`
#' is returned for future use (e.g. anisotropic reflected radiation).
#'
#' @param datetime POSIXct datetime vector. Sub-daily precision is used; times
#'   should be in UTC (or have a UTC offset applied) for correct results.
#' @param latitude Latitude in decimal degrees (-90 to 90).
#' @param longitude Longitude in decimal degrees (-180 to 180).
#'
#' @return A named list with two numeric vectors:
#' \describe{
#'   \item{zenith}{Solar zenith angle in degrees (0 = overhead, 90 = horizon).}
#'   \item{azimuth}{Solar azimuth angle in degrees clockwise from north
#'     (0/360 = north, 90 = east, 180 = south, 270 = west). Currently unused
#'     by callers but available for future directional radiation models.}
#' }
#' @keywords internal
calculate_solar_position <- function(datetime, latitude, longitude) {
  # Convert to Julian Day
  JD <- as.numeric(julian(datetime, origin = as.POSIXct("2000-01-01 12:00:00", tz = "UTC"))) + 2451545.0

  # Days since J2000.0
  n <- JD - 2451545.0

  # Mean longitude (degrees)
  L <- (280.460 + 0.9856474 * n) %% 360

  # Mean anomaly (degrees)
  g <- (357.528 + 0.9856003 * n) %% 360
  g_rad <- g * pi / 180

  # Ecliptic longitude (degrees)
  lambda <- (L + 1.915 * sin(g_rad) + 0.020 * sin(2 * g_rad)) %% 360
  lambda_rad <- lambda * pi / 180

  # Obliquity of ecliptic (degrees)
  epsilon <- 23.439 - 0.0000004 * n
  epsilon_rad <- epsilon * pi / 180

  # Right ascension (degrees)
  RA <- atan2(cos(epsilon_rad) * sin(lambda_rad), cos(lambda_rad)) * 180 / pi
  RA <- (RA + 360) %% 360

  # Declination (degrees)
  delta <- asin(sin(epsilon_rad) * sin(lambda_rad)) * 180 / pi
  delta_rad <- delta * pi / 180

  # Greenwich Mean Sidereal Time (degrees) — already encodes the full UT
  # from n (Julian days since J2000.0, including fractional day).
  # Local Sidereal Time (LST) = GMST + east longitude.
  # Note: adding hour * 15 here would double-count the time-of-day, since n
  # already carries the sub-day fraction via the Julian Day.
  GMST <- (280.460 + 360.9856474 * n) %% 360
  LST  <- (GMST + longitude) %% 360

  # Local Hour Angle: angle of the sun west of the observer's meridian
  H <- (LST - RA) %% 360
  H <- ifelse(H > 180, H - 360, H)
  H_rad <- H * pi / 180

  # Latitude in radians
  lat_rad <- latitude * pi / 180

  # Solar zenith angle
  cos_zenith <- sin(lat_rad) * sin(delta_rad) + cos(lat_rad) * cos(delta_rad) * cos(H_rad)
  cos_zenith <- pmax(-1, pmin(1, cos_zenith))
  zenith <- acos(cos_zenith) * 180 / pi

  # Solar azimuth angle
  sin_azimuth <- -cos(delta_rad) * sin(H_rad) / sin(zenith * pi / 180)
  cos_azimuth <- (sin(delta_rad) - sin(lat_rad) * cos(zenith * pi / 180)) /
    (cos(lat_rad) * sin(zenith * pi / 180))

  # Handle edge cases
  sin_azimuth <- pmax(-1, pmin(1, sin_azimuth))
  cos_azimuth <- pmax(-1, pmin(1, cos_azimuth))

  azimuth <- atan2(sin_azimuth, cos_azimuth) * 180 / pi
  azimuth <- (azimuth + 360) %% 360

  list(zenith = zenith, azimuth = azimuth)
}

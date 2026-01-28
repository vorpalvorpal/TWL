#' Fetch weather data from the Open-Meteo API
#'
#' Queries the Open-Meteo forecast or historical archive API to retrieve
#' hourly weather observations for the given location and times. Only the
#' variables listed in `fields` are requested.
#'
#' @param datetime POSIXct datetime vector.
#' @param latitude Latitude in decimal degrees.
#' @param longitude Longitude in decimal degrees.
#' @param fields Character vector of Open-Meteo hourly variable names to
#'   request (e.g. `"temperature_2m"`, `"wind_speed_10m"`).
#' @param verbose If `TRUE`, report progress.
#'
#' @return A named list of numeric vectors, one per requested field, aligned
#'   to the input `datetime` values. Values that cannot be matched are `NA`.
#' @keywords internal
fetch_openmeteo <- function(datetime, latitude, longitude, fields,
                            verbose = FALSE) {
  if (length(fields) == 0L) return(list())

  dates <- as.Date(datetime)
  date_min <- min(dates, na.rm = TRUE)
  date_max <- max(dates, na.rm = TRUE)
  today <- Sys.Date()

  fields_csv <- paste(fields, collapse = ",")

  # Choose endpoint based on date range
  if (date_max <= today + 16L && date_min >= today - 92L) {
    # Forecast API can serve past_days up to 92 and forecast_days up to 16
    past <- as.integer(today - date_min)
    fore <- as.integer(date_max - today) + 1L
    past <- max(0L, min(92L, past))
    fore <- max(1L, min(16L, fore))
    url <- sprintf(
      "https://api.open-meteo.com/v1/forecast?latitude=%.6f&longitude=%.6f&hourly=%s&wind_speed_unit=ms&past_days=%d&forecast_days=%d&timeformat=iso8601",
      latitude, longitude, fields_csv, past, fore
    )
  } else {
    # Historical archive API
    url <- sprintf(
      "https://archive-api.open-meteo.com/v1/archive?latitude=%.6f&longitude=%.6f&hourly=%s&wind_speed_unit=ms&start_date=%s&end_date=%s&timeformat=iso8601",
      latitude, longitude, fields_csv,
      format(date_min, "%Y-%m-%d"), format(date_max, "%Y-%m-%d")
    )
  }

  if (verbose) {
    cli_alert("Fetching weather data from Open-Meteo...")
  }

  resp <- tryCatch(
    jsonlite::fromJSON(url),
    error = function(e) {
      stop(
        "Failed to fetch data from Open-Meteo: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (isTRUE(resp$error)) {
    stop("Open-Meteo API error: ", resp$reason, call. = FALSE)
  }

  # Parse API timestamps and match to input datetimes
  api_times <- as.POSIXct(resp$hourly$time, format = "%Y-%m-%dT%H:%M", tz = "UTC")
  input_hours <- as.POSIXct(
    format(datetime, "%Y-%m-%d %H:00:00"), tz = "UTC"
  )
  idx <- match(input_hours, api_times)

  result <- lapply(fields, function(f) {
    vals <- resp$hourly[[f]]
    if (is.null(vals)) {
      rep(NA_real_, length(datetime))
    } else {
      as.numeric(vals[idx])
    }
  })
  names(result) <- fields

  n_matched <- sum(!is.na(idx))
  if (verbose) {
    cli_alert_success(
      "Matched {n_matched}/{length(datetime)} observation{?s} to API data"
    )
  }

  result
}

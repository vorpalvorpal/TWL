#' Fetch weather data from the Open-Meteo API
#'
#' Queries the Open-Meteo forecast or historical archive API to retrieve
#' hourly weather observations for the given location and times. Only the
#' variables listed in `fields` are requested.
#'
#' Input `datetime` values must be in UTC (or a UTC-equivalent timezone).
#' Sub-hourly timestamps are truncated to the hour before matching against
#' API data.
#'
#' @param datetime POSIXct datetime vector (UTC).
#' @param latitude Latitude in decimal degrees.
#' @param longitude Longitude in decimal degrees.
#' @param fields Character vector of Open-Meteo hourly variable names to
#'   request (e.g. `"temperature_2m"`, `"wind_speed_10m"`).
#' @param verbose If `TRUE`, report progress.
#'
#' @return A named list of numeric vectors, one per requested field, aligned
#'   to the input `datetime` values. Values that cannot be matched are `NA`.
#' @keywords internal
fetch_openmeteo <- function(
  datetime,
  latitude,
  longitude,
  fields,
  verbose = FALSE
) {
  if (length(fields) == 0L) {
    return(list())
  }

  # Always work in UTC to avoid session-timezone mismatches
  datetime_utc <- as.POSIXct(format(datetime, tz = "UTC"), tz = "UTC")
  dates <- as.Date(datetime_utc, tz = "UTC")
  date_min <- min(dates, na.rm = TRUE)
  date_max <- max(dates, na.rm = TRUE)
  today <- Sys.Date()

  fields_csv <- paste(fields, collapse = ",")

  # Choose endpoint based on date range.
  # Forecast API: past 92 days + up to 16 days ahead.
  # Archive API: everything else.
  if (date_max <= today + 16L && date_min >= today - 92L) {
    past <- max(0L, min(92L, as.integer(today - date_min)))
    fore <- max(1L, min(16L, as.integer(date_max - today) + 1L))
    base_url <- sprintf(
      "https://api.open-meteo.com/v1/forecast?latitude=%.6f&longitude=%.6f&hourly=%s&wind_speed_unit=ms&past_days=%d&forecast_days=%d&timeformat=iso8601&timezone=UTC",
      latitude, longitude, fields_csv, past, fore
    )
  } else {
    base_url <- sprintf(
      "https://archive-api.open-meteo.com/v1/archive?latitude=%.6f&longitude=%.6f&hourly=%s&wind_speed_unit=ms&start_date=%s&end_date=%s&timeformat=iso8601&timezone=UTC",
      latitude, longitude, fields_csv,
      format(date_min, "%Y-%m-%d"),
      format(date_max, "%Y-%m-%d")
    )
  }

  if (verbose) {
    cli_alert("Fetching weather data from Open-Meteo...")
  }

  # Use httr2 for proper HTTP handling (timeout, retry, user-agent)
  resp <- tryCatch({
    req <- httr2::request(base_url) |>
      httr2::req_user_agent("TWL R package") |>
      httr2::req_timeout(30) |>
      httr2::req_retry(max_tries = 3, backoff = ~ 2)

    raw <- httr2::req_perform(req)
    httr2::resp_body_json(raw, simplifyVector = TRUE)
  }, error = function(e) {
    stop(
      "Failed to fetch data from Open-Meteo: ",
      conditionMessage(e),
      call. = FALSE
    )
  })

  if (isTRUE(resp$error)) {
    stop("Open-Meteo API error: ", resp$reason, call. = FALSE)
  }

  # Parse API timestamps (UTC) and match to input datetimes truncated to hour
  api_times <- as.POSIXct(
    resp$hourly$time,
    format = "%Y-%m-%dT%H:%M",
    tz = "UTC"
  )
  input_hours <- as.POSIXct(
    format(datetime_utc, "%Y-%m-%d %H:00:00", tz = "UTC"),
    tz = "UTC"
  )
  idx <- match(input_hours, api_times)

  n_matched <- sum(!is.na(idx))
  n_total   <- length(datetime)

  if (verbose) {
    cli_alert_success(
      "Matched {n_matched}/{n_total} observation{?s} to API data"
    )
  }

  n_unmatched <- sum(is.na(idx))
  if (n_unmatched > 0L) {
    cli_alert_warning(
      "{n_unmatched} observation{?s} could not be matched to API data and will be NA"
    )
  }

  result <- lapply(fields, function(f) {
    vals <- resp$hourly[[f]]
    if (is.null(vals)) {
      rep(NA_real_, n_total)
    } else {
      raw_vals <- as.numeric(vals[idx])
      n_na <- sum(is.na(raw_vals))
      if (n_na > 0L && verbose) {
        cli_alert_warning(
          "Field '{f}': {n_na} matched value{?s} are NA (may be at archive boundary)"
        )
      }
      raw_vals
    }
  })
  names(result) <- fields

  result
}

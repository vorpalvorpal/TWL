#' Categorise TWL values into work restriction zones
#'
#' Maps TWL values to the four standard zones defined by Brake & Bates (2002).
#'
#' @param twl Numeric vector of TWL values in W/m^2.
#' @return Character vector of categories: `"Withdrawal"`, `"Buffer"`,
#'   `"Acclimatisation"`, or `"Unrestricted"`.
#'
#' @examples
#' categorise_twl(c(100, 130, 180, 250))
#'
#' @export
categorise_twl <- function(twl) {
  case_when(
    is.na(twl) ~ NA_character_,
    twl < 115 ~ "Withdrawal",
    twl < 140 ~ "Buffer",
    twl < 220 ~ "Acclimatisation",
    TRUE ~ "Unrestricted"
  )
}

#' Get TWL category colour
#'
#' Returns hex colour codes corresponding to each TWL zone for use in
#' visualisations.
#'
#' @param twl Numeric vector of TWL values in W/m^2.
#' @return Character vector of hex colour codes:
#' \describe{
#'   \item{`"#D32F2F"`}{Red -- Withdrawal}
#'   \item{`"#FF9800"`}{Orange -- Buffer}
#'   \item{`"#FFC107"`}{Amber -- Acclimatisation}
#'   \item{`"#4CAF50"`}{Green -- Unrestricted}
#'   \item{`"#CCCCCC"`}{Grey -- NA}
#' }
#'
#' @examples
#' twl_colour(c(100, 130, 180, 250, NA))
#'
#' @export
twl_colour <- function(twl) {
  case_when(
    is.na(twl) ~ "#CCCCCC",
    twl < 115 ~ "#D32F2F",
    twl < 140 ~ "#FF9800",
    twl < 220 ~ "#FFC107",
    TRUE ~ "#4CAF50"
  )
}

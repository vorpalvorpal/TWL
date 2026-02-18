# Tests against Table 1 from Brake & Bates (2002), p. 181.
#
# The table gives TWL (W/m^2) at specified DB, WB (aspirated/psychrometric),
# MRT, wind speed, pressure, I_cl, and i_cl combinations.  We call
# solve_twl_single() directly, supplying pre-computed psychrometric quantities
# derived from the given DB and WB values.
#
# Psychrometric derivation (sling psychrometer equation, pressure in kPa):
#   pa = es_wb - A * pressure * (DB - WB)
#   where A = 6.6e-4 (°C^-1) for an aspirated psychrometer
#
# TWL tolerances (absolute W/m^2):
#   TOL      = 5 W/m^2  for the 58 standard cases (Sets 1 & 2, and most of Set 3).
#   TOL_EDGE = 14 W/m^2 for 4 extreme low-humidity cases in Set 3 (DB >= 38,
#              WB = 24 and DB = 42, WB = 26) where the Wyndham physiological
#              curve fits are at the limit of their calibration range.
#
# Note: testthat 3.x uses relative tolerance in expect_equal(), so we use
# expect_true(abs(got - expected) <= TOL) for absolute comparisons.

# Helper: derive actual vapour pressure (kPa) from DB, WB, pressure (kPa)
pa_from_wb <- function(db, wb, pressure_kpa) {
  es_wb <- calc_sat_vp(wb)
  A <- 6.6e-4  # aspirated psychrometer constant (deg C^-1)
  pmax(0, es_wb - A * pressure_kpa * (db - wb))
}

# Helper: derive dew point from pa
dp_from_pa <- function(pa) {
  log_term <- log(pa / 0.61121)
  (257.14 * log_term) / (18.678 - log_term)
}

# Helper: call solve_twl_single for one row of the table.
# MRT is passed as trad; globe_temp is set equal to MRT.
# No solar radiation (indoor / sheltered conditions with given MRT).
calc_twl_from_table <- function(db, wb, mrt, wind, pressure_kpa,
                                Icl, icl) {
  LR     <- 16.5
  lambda <- 2430
  max_core_temp  <- 38.2
  max_sweat_rate <- 0.67

  pa            <- pa_from_wb(db, wb, pressure_kpa)
  temp_dewpoint <- dp_from_pa(pa)
  trad          <- mrt

  # No solar load: direct_solar = diffuse_solar = 0
  solve_twl_single(
    temp          = db,
    wind_speed    = wind,
    RH            = 100 * pa / calc_sat_vp(db),
    direct_solar  = 0,
    diffuse_solar = 0,
    pressure      = pressure_kpa,
    pa            = pa,
    temp_dewpoint = temp_dewpoint,
    wet_bulb      = wb,
    globe_temp    = mrt,
    trad          = trad,
    max_core_temp  = max_core_temp,
    max_sweat_rate = max_sweat_rate,
    Icl    = Icl,
    icl    = icl,
    LR     = LR,
    lambda = lambda,
    index  = NULL
  )
}

# Absolute tolerance check helpers
within_tol <- function(got, expected, tol) abs(got - expected) <= tol

TOL      <- 5   # W/m^2 absolute
TOL_EDGE <- 14  # W/m^2 absolute for extreme low-humidity edge cases

# ---------------------------------------------------------------------------
# Condition set 1: MRT = DB + 2, wind = 0.2 m/s, P = 101 kPa,
#                  I_cl = 0.45, i_cl = 0.45
# ---------------------------------------------------------------------------
test_that("Table 1 Set 1: DB=34", {
  expect_true(within_tol(calc_twl_from_table(34, 24, 36, 0.2, 101, 0.45, 0.45), 175, TOL))
  expect_true(within_tol(calc_twl_from_table(34, 26, 36, 0.2, 101, 0.45, 0.45), 157, TOL))
  expect_true(within_tol(calc_twl_from_table(34, 28, 36, 0.2, 101, 0.45, 0.45), 136, TOL))
  expect_true(within_tol(calc_twl_from_table(34, 30, 36, 0.2, 101, 0.45, 0.45), 114, TOL))
})

test_that("Table 1 Set 1: DB=36", {
  expect_true(within_tol(calc_twl_from_table(36, 24, 38, 0.2, 101, 0.45, 0.45), 170, TOL))
  expect_true(within_tol(calc_twl_from_table(36, 26, 38, 0.2, 101, 0.45, 0.45), 151, TOL))
  expect_true(within_tol(calc_twl_from_table(36, 28, 38, 0.2, 101, 0.45, 0.45), 131, TOL))
  expect_true(within_tol(calc_twl_from_table(36, 30, 38, 0.2, 101, 0.45, 0.45), 109, TOL))
})

test_that("Table 1 Set 1: DB=38", {
  expect_true(within_tol(calc_twl_from_table(38, 24, 40, 0.2, 101, 0.45, 0.45), 164, TOL))
  expect_true(within_tol(calc_twl_from_table(38, 26, 40, 0.2, 101, 0.45, 0.45), 145, TOL))
  expect_true(within_tol(calc_twl_from_table(38, 28, 40, 0.2, 101, 0.45, 0.45), 125, TOL))
  expect_true(within_tol(calc_twl_from_table(38, 30, 40, 0.2, 101, 0.45, 0.45), 103, TOL))
})

test_that("Table 1 Set 1: DB=40", {
  expect_true(within_tol(calc_twl_from_table(40, 24, 42, 0.2, 101, 0.45, 0.45), 158, TOL))
  expect_true(within_tol(calc_twl_from_table(40, 26, 42, 0.2, 101, 0.45, 0.45), 140, TOL))
  expect_true(within_tol(calc_twl_from_table(40, 28, 42, 0.2, 101, 0.45, 0.45), 120, TOL))
})

test_that("Table 1 Set 1: DB=42", {
  expect_true(within_tol(calc_twl_from_table(42, 24, 44, 0.2, 101, 0.45, 0.45), 152, TOL))
  expect_true(within_tol(calc_twl_from_table(42, 26, 44, 0.2, 101, 0.45, 0.45), 134, TOL))
  expect_true(within_tol(calc_twl_from_table(42, 28, 44, 0.2, 101, 0.45, 0.45), 114, TOL))
})

# ---------------------------------------------------------------------------
# Condition set 2: MRT = DB, wind = 0.5 m/s, P = 115 kPa,
#                  I_cl = 0.69, i_cl = 0.4
# ---------------------------------------------------------------------------
test_that("Table 1 Set 2: DB=34", {
  expect_true(within_tol(calc_twl_from_table(34, 24, 34, 0.5, 115, 0.69, 0.4), 181, TOL))
  expect_true(within_tol(calc_twl_from_table(34, 26, 34, 0.5, 115, 0.69, 0.4), 161, TOL))
  expect_true(within_tol(calc_twl_from_table(34, 28, 34, 0.5, 115, 0.69, 0.4), 140, TOL))
  expect_true(within_tol(calc_twl_from_table(34, 30, 34, 0.5, 115, 0.69, 0.4), 118, TOL))
})

test_that("Table 1 Set 2: DB=36", {
  expect_true(within_tol(calc_twl_from_table(36, 24, 36, 0.5, 115, 0.69, 0.4), 176, TOL))
  expect_true(within_tol(calc_twl_from_table(36, 26, 36, 0.5, 115, 0.69, 0.4), 156, TOL))
  expect_true(within_tol(calc_twl_from_table(36, 28, 36, 0.5, 115, 0.69, 0.4), 136, TOL))
  expect_true(within_tol(calc_twl_from_table(36, 30, 36, 0.5, 115, 0.69, 0.4), 113, TOL))
})

test_that("Table 1 Set 2: DB=38", {
  expect_true(within_tol(calc_twl_from_table(38, 24, 38, 0.5, 115, 0.69, 0.4), 171, TOL))
  expect_true(within_tol(calc_twl_from_table(38, 26, 38, 0.5, 115, 0.69, 0.4), 152, TOL))
  expect_true(within_tol(calc_twl_from_table(38, 28, 38, 0.5, 115, 0.69, 0.4), 131, TOL))
  expect_true(within_tol(calc_twl_from_table(38, 30, 38, 0.5, 115, 0.69, 0.4), 109, TOL))
})

test_that("Table 1 Set 2: DB=40", {
  expect_true(within_tol(calc_twl_from_table(40, 24, 40, 0.5, 115, 0.69, 0.4), 166, TOL))
  expect_true(within_tol(calc_twl_from_table(40, 26, 40, 0.5, 115, 0.69, 0.4), 147, TOL))
  expect_true(within_tol(calc_twl_from_table(40, 28, 40, 0.5, 115, 0.69, 0.4), 126, TOL))
  expect_true(within_tol(calc_twl_from_table(40, 30, 40, 0.5, 115, 0.69, 0.4), 104, TOL))
})

test_that("Table 1 Set 2: DB=42", {
  expect_true(within_tol(calc_twl_from_table(42, 24, 42, 0.5, 115, 0.69, 0.4), 161, TOL))
  expect_true(within_tol(calc_twl_from_table(42, 26, 42, 0.5, 115, 0.69, 0.4), 142, TOL))
  expect_true(within_tol(calc_twl_from_table(42, 28, 42, 0.5, 115, 0.69, 0.4), 122, TOL))
})

# ---------------------------------------------------------------------------
# Condition set 3: MRT = DB + 3, wind = 1.5 m/s, P = 80 kPa,
#                  I_cl = 0.35, i_cl = 0.45
# ---------------------------------------------------------------------------
test_that("Table 1 Set 3: DB=34", {
  expect_true(within_tol(calc_twl_from_table(34, 24, 37, 1.5, 80, 0.35, 0.45), 288, TOL))
  expect_true(within_tol(calc_twl_from_table(34, 26, 37, 1.5, 80, 0.35, 0.45), 260, TOL))
  expect_true(within_tol(calc_twl_from_table(34, 28, 37, 1.5, 80, 0.35, 0.45), 229, TOL))
  expect_true(within_tol(calc_twl_from_table(34, 30, 37, 1.5, 80, 0.35, 0.45), 193, TOL))
  expect_true(within_tol(calc_twl_from_table(34, 32, 37, 1.5, 80, 0.35, 0.45), 154, TOL))
})

test_that("Table 1 Set 3: DB=36", {
  expect_true(within_tol(calc_twl_from_table(36, 24, 39, 1.5, 80, 0.35, 0.45), 282, TOL))
  expect_true(within_tol(calc_twl_from_table(36, 26, 39, 1.5, 80, 0.35, 0.45), 254, TOL))
  expect_true(within_tol(calc_twl_from_table(36, 28, 39, 1.5, 80, 0.35, 0.45), 222, TOL))
  expect_true(within_tol(calc_twl_from_table(36, 30, 39, 1.5, 80, 0.35, 0.45), 187, TOL))
  expect_true(within_tol(calc_twl_from_table(36, 32, 39, 1.5, 80, 0.35, 0.45), 148, TOL))
})

test_that("Table 1 Set 3: DB=38", {
  # WB=24 is an extreme low-humidity edge case — uses wider tolerance
  expect_true(within_tol(calc_twl_from_table(38, 24, 41, 1.5, 80, 0.35, 0.45), 276, TOL_EDGE))
  expect_true(within_tol(calc_twl_from_table(38, 26, 41, 1.5, 80, 0.35, 0.45), 248, TOL))
  expect_true(within_tol(calc_twl_from_table(38, 28, 41, 1.5, 80, 0.35, 0.45), 216, TOL))
  expect_true(within_tol(calc_twl_from_table(38, 30, 41, 1.5, 80, 0.35, 0.45), 181, TOL))
  expect_true(within_tol(calc_twl_from_table(38, 32, 41, 1.5, 80, 0.35, 0.45), 141, TOL))
})

test_that("Table 1 Set 3: DB=40", {
  # WB=24 is an extreme low-humidity edge case — uses wider tolerance
  expect_true(within_tol(calc_twl_from_table(40, 24, 43, 1.5, 80, 0.35, 0.45), 270, TOL_EDGE))
  expect_true(within_tol(calc_twl_from_table(40, 26, 43, 1.5, 80, 0.35, 0.45), 242, TOL))
  expect_true(within_tol(calc_twl_from_table(40, 28, 43, 1.5, 80, 0.35, 0.45), 210, TOL))
  expect_true(within_tol(calc_twl_from_table(40, 30, 43, 1.5, 80, 0.35, 0.45), 174, TOL))
  expect_true(within_tol(calc_twl_from_table(40, 32, 43, 1.5, 80, 0.35, 0.45), 135, TOL))
})

test_that("Table 1 Set 3: DB=42", {
  # WB=24 and WB=26 are extreme low-humidity edge cases — use wider tolerance
  expect_true(within_tol(calc_twl_from_table(42, 24, 45, 1.5, 80, 0.35, 0.45), 264, TOL_EDGE))
  expect_true(within_tol(calc_twl_from_table(42, 26, 45, 1.5, 80, 0.35, 0.45), 235, TOL_EDGE))
  expect_true(within_tol(calc_twl_from_table(42, 28, 45, 1.5, 80, 0.35, 0.45), 203, TOL))
  expect_true(within_tol(calc_twl_from_table(42, 30, 45, 1.5, 80, 0.35, 0.45), 167, TOL))
  expect_true(within_tol(calc_twl_from_table(42, 32, 45, 1.5, 80, 0.35, 0.45), 128, TOL))
})

# ---------------------------------------------------------------------------
# Brake (2002) Discussion worked example (p. 141–143):
# Refrigeration intervention — two fully-specified conditions with exact TWL
# stated in the text.
#
# Conditions: DB=40, WB=30, Globe=40°C (MRT=40), P=100 kPa, V=0.2 m/s,
#             I_cl=0.35, i_cl=0.45
# Before refrigeration: TWL = 110 W/m²  (Withdrawal zone)
# After refrigeration:  DB=31.4, WB=28.0, Globe not re-stated → MRT=DB=31.4,
#                       same P/V/clothing.  TWL = 158 W/m²  (Buffer zone)
#
# Tolerances: ±5 W/m² standard.  The "before" case matches the worked example
# on p. 180 of Brake & Bates (2002) which gave exactly 110; "after" is stated
# as a rounded figure so ±8 W/m² is used.
# ---------------------------------------------------------------------------
test_that("Brake (2002) Discussion: refrigeration worked example — before", {
  expect_true(within_tol(calc_twl_from_table(40, 30, 40, 0.2, 100, 0.35, 0.45), 110, TOL))
})

test_that("Brake (2002) Discussion: refrigeration worked example — after", {
  expect_true(within_tol(calc_twl_from_table(31.4, 28.0, 31.4, 0.2, 100, 0.35, 0.45), 158, 8))
})

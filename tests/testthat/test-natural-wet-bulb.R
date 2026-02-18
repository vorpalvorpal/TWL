# Tests for solve_natural_wb_single() and calculate_natural_wet_bulb()
#
# Physical properties of the natural (unventilated) wet bulb temperature:
#
#  1. At RH = 100 %, WBn == DB (no evaporative cooling possible).
#  2. WBn is always less than DB.
#  3. When globe_temp == DB (no net radiation), WBn is very close to the
#     aspirated (psychrometric) wet bulb -- within ~0.3 °C.
#  4. When globe_temp > DB (radiative loading), WBn > aspirated WB.
#  5. Higher globe temperature raises WBn (more radiant heat input).
#  6. Higher wind speed lowers the radiation-driven excess, bringing WBn
#     closer to the aspirated WB.
#  7. Higher pressure raises WBn (lower pressure accelerates evaporation).
#  8. Higher RH raises WBn (monotone relationship).
#  9. NA inputs propagate to NA output.
# 10. The vectorised wrapper calculate_natural_wet_bulb() agrees with the
#     scalar solver element-wise.

TOL_TIGHT <- 0.05   # °C -- near-exact physical limits
TOL_LOOSE <- 0.3    # °C -- allowable spread from aspirated WB at globe≈DB

# ── 1. RH = 100 %: natural WB equals DB ────────────────────────────────────
test_that("RH=100% gives natural WB equal to DB", {
  for (db in c(25, 30, 35, 40)) {
    wb <- solve_natural_wb_single(db, 100, 101.325, 0.5, db)
    expect_equal(wb, db, tolerance = TOL_TIGHT,
                 label = paste("DB =", db))
  }
})

# ── 2. Natural WB is always strictly less than DB ──────────────────────────
test_that("natural WB is always less than DB", {
  cases <- list(
    list(db = 30, rh = 60, globe = 35, wind = 0.5),
    list(db = 35, rh = 50, globe = 40, wind = 0.5),
    list(db = 40, rh = 40, globe = 45, wind = 1.0),
    list(db = 40, rh = 70, globe = 50, wind = 0.2),
    list(db = 25, rh = 30, globe = 30, wind = 2.0)
  )
  for (cx in cases) {
    wb <- solve_natural_wb_single(cx$db, cx$rh, 101.325, cx$wind, cx$globe)
    expect_lt(wb, cx$db,
              label = paste("DB =", cx$db, "RH =", cx$rh))
  }
})

# ── 3. Globe == DB: WBn close to aspirated WB ──────────────────────────────
test_that("when globe == DB, natural WB is within 0.3 °C of aspirated WB", {
  for (rh in c(40, 60, 80)) {
    db    <- 35
    wb_nat <- solve_natural_wb_single(db, rh, 101.325, 0.5, db)
    wb_asp <- calc_aspirated_wb(db, rh)
    expect_equal(wb_nat, wb_asp, tolerance = TOL_LOOSE,
                 label = paste("RH =", rh))
  }
})

# ── 4. Globe > DB: natural WB exceeds aspirated WB ─────────────────────────
test_that("globe_temp > DB causes natural WB to exceed aspirated WB", {
  cases <- list(
    list(db = 30, rh = 60, globe = 40, wind = 0.5),
    list(db = 35, rh = 50, globe = 50, wind = 0.5),
    list(db = 40, rh = 40, globe = 55, wind = 0.5),
    list(db = 38, rh = 70, globe = 48, wind = 0.3)
  )
  for (cx in cases) {
    wb_nat <- solve_natural_wb_single(cx$db, cx$rh, 101.325, cx$wind, cx$globe)
    wb_asp <- calc_aspirated_wb(cx$db, cx$rh)
    expect_gt(wb_nat, wb_asp,
              label = paste("DB =", cx$db, "RH =", cx$rh, "globe =", cx$globe))
  }
})

# ── 5. Higher globe temperature raises natural WB ──────────────────────────
test_that("increasing globe temperature raises natural WB monotonically", {
  db <- 35; rh <- 50; pressure <- 101.325; wind <- 0.5
  globe_temps <- c(35, 40, 45, 50, 55)
  wbs <- vapply(globe_temps, function(g)
    solve_natural_wb_single(db, rh, pressure, wind, g), numeric(1))
  # Each successive value should be higher
  expect_true(all(diff(wbs) > 0),
              info = paste("WBs:", paste(round(wbs, 3), collapse = ", ")))
})

# ── 6. Higher wind speed reduces the radiative excess ──────────────────────
# With globe > DB the natural WB excess over aspirated WB should
# shrink as wind increases (convection overwhelms radiation).
test_that("increasing wind speed reduces radiative excess in natural WB", {
  db <- 35; rh <- 50; pressure <- 101.325; globe <- 50
  wb_asp <- calc_aspirated_wb(db, rh)
  wind_speeds <- c(0.2, 0.5, 1.0, 2.0, 4.0)
  excesses <- vapply(wind_speeds, function(ws) {
    solve_natural_wb_single(db, rh, pressure, ws, globe) - wb_asp
  }, numeric(1))
  # Excess should decrease monotonically with wind
  expect_true(all(diff(excesses) < 0),
              info = paste("Excesses:", paste(round(excesses, 3), collapse = ", ")))
})

# ── 7. Higher pressure raises natural WB ───────────────────────────────────
test_that("higher pressure raises natural WB (less evaporation)", {
  db <- 35; rh <- 50; wind <- 0.5; globe <- 45
  pressures <- c(80, 90, 101.325, 115)
  wbs <- vapply(pressures, function(p)
    solve_natural_wb_single(db, rh, p, wind, globe), numeric(1))
  expect_true(all(diff(wbs) > 0),
              info = paste("WBs:", paste(round(wbs, 3), collapse = ", ")))
})

# ── 8. Higher RH raises natural WB ─────────────────────────────────────────
test_that("increasing RH raises natural WB monotonically", {
  db <- 35; pressure <- 101.325; wind <- 0.5; globe <- 40
  rhs <- c(30, 40, 50, 60, 70, 80, 90)
  wbs <- vapply(rhs, function(rh)
    solve_natural_wb_single(db, rh, pressure, wind, globe), numeric(1))
  expect_true(all(diff(wbs) > 0),
              info = paste("WBs:", paste(round(wbs, 3), collapse = ", ")))
})

# ── 9. NA inputs propagate ─────────────────────────────────────────────────
test_that("NA inputs return NA", {
  expect_true(is.na(solve_natural_wb_single(NA,  50, 101.325, 0.5, 35)))
  expect_true(is.na(solve_natural_wb_single(35,  NA, 101.325, 0.5, 35)))
  expect_true(is.na(solve_natural_wb_single(35,  50,      NA, 0.5, 35)))
  expect_true(is.na(solve_natural_wb_single(35,  50, 101.325,  NA, 35)))
  expect_true(is.na(solve_natural_wb_single(35,  50, 101.325, 0.5, NA)))
})

# ── 10. Vectorised wrapper matches scalar solver ────────────────────────────
test_that("calculate_natural_wet_bulb matches scalar solver element-wise", {
  temps  <- c(30, 35, 40)
  rhs    <- c(60, 50, 40)
  press  <- rep(101.325, 3)
  winds  <- c(0.5, 1.0, 0.2)
  globes <- c(35, 45, 50)

  scalar_results <- mapply(solve_natural_wb_single, temps, rhs, press, winds, globes)
  vector_results <- calculate_natural_wet_bulb(temps, rhs, press, winds, globes)

  expect_equal(vector_results, scalar_results, tolerance = 1e-10)
})

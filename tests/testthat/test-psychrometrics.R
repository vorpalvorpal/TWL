test_that("calc_sat_vp returns known values", {
  # At 20C, saturation VP is approximately 2.34 kPa

  expect_equal(calc_sat_vp(20), 2.338, tolerance = 0.01)
  # At 0C, approximately 0.611 kPa
  expect_equal(calc_sat_vp(0), 0.61121, tolerance = 0.001)
})

test_that("calc_dew_point is consistent with calc_sat_vp", {
  # At 100% RH, dew point should equal air temperature
  expect_equal(calc_dew_point(25, 100), 25, tolerance = 0.1)
})

test_that("calc_aspirated_wb is less than or equal to temp", {
  wb <- calc_aspirated_wb(35, 50)
  expect_true(wb <= 35)
  expect_true(wb > 0)
})

test_that("categorise_twl returns correct categories", {
  cats <- categorise_twl(c(100, 130, 180, 250, NA))
  expect_equal(cats, c("Withdrawal", "Buffer", "Acclimatisation", "Unrestricted", NA))
})

test_that("twl_colour returns correct hex colours", {
  cols <- twl_colour(c(100, 130, 180, 250, NA))
  expect_equal(cols, c("#D32F2F", "#FF9800", "#FFC107", "#4CAF50", "#CCCCCC"))
})

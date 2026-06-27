
# batch_fit tests  -------------------------------------------------------------



# .stan_data_prep tests  -------------------------------------------------------



# .stan_file tests  ------------------------------------------------------------

test_that("stan_file returns an vb existing linear file", {
  
  path <- .stan_file("vb", "linear")
  
  expect_true(file.exists(path))
  
})

test_that("stan_file returns an existing categorical file", {
  
  path <- .stan_file("vb", "categorical")
  
  expect_true(file.exists(path))
  
})

test_that("stan_file returns an existing random effect file", {
  
  path <- .stan_file("vb", "random")
  
  expect_true(file.exists(path))
  
})

test_that("invalid model throws an error", {
  
  expect_error(
    .stan_file("fake_model", "linear")
  )
  
})

test_that("invalid effect structure throws an error", {
  
  expect_error(
    .stan_file("vb", "fake_effect")
  )
  
})



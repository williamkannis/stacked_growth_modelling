
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


# .stan_diag tests  ------------------------------------------------------------

.diag_test_helper <- function(out.type) {
  if(out.type == "good") file <- "good_out.rds"
  if(out.type == "bad") file <- "bad_out.rds"
  
  out <- readRDS(
    test_path(
      "test_data",
      file
    )
  )
  
  diag <- .stan_diag(out)
  
}

test_that("returns list", {
  # RUN MODEL SO OUTPUT IS CEAN AND STORE IN TEST DATA
  
  # diagnostics
  diag <- .diag_test_helper("good")
  expect_true(is.list(diag))
  
})

test_that("returns list", {
  # RUN MODEL SO OUTPUT IS CEAN AND STORE IN TEST DATA
  
  # diagnostics
  diag <- .diag_test_helper("good")
  expect_true(is.list(diag))
  
})


test_that("model outputs with no issues are detected", {
  # RUN MODEL SO OUTPUT IS CEAN AND STORE IN TEST DATA
  
  # diagnostics
  diag <- .diag_test_helper("good")
  expect_true(diag$divergence > 0)
  
})

test_that("model outputs with divergent transitions are detected", {
  # RUN MODEL SO OUTPUT HAS DIVERGNET AND TREE TRANSIONS AND STORE IN TEST DATA
})

test_that("model outputs with large tree depths are detected", {
  # RUN MODEL SO OUTPUT HAS DIVERGNET AND TREE TRANSIONS AND STORE IN TEST DATA
})

test_that("model outputs with small effect sample sizes are detected", {
  # RUN MODEL SO OUTPUT HAS DIVERGNET AND TREE TRANSIONS AND STORE IN TEST DATA
})

test_that("model outputs with high rhat are detected", {
  # RUN MODEL SO OUTPUT HAS DIVERGNET AND TREE TRANSIONS AND STORE IN TEST DATA
})



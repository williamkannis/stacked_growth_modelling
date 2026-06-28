# stan_diag tests  ------------------------------------------------------------

# .diag_test_helper <- function(out.type) {
#   if(out.type == "good") file <- "good_out.rds"
#   if(out.type == "bad") file <- "bad_out.rds"
#   
#   out <- readRDS(
#     test_path(
#       "test_data",
#       file
#     )
#   )
#   
#   diag <- .stan_diag(out)
#   
# }
# 
# test_that("returns list", {
#   # RUN MODEL SO OUTPUT IS CEAN AND STORE IN TEST DATA
#   
#   # diagnostics
#   diag <- .diag_test_helper("good")
#   expect_true(is.list(diag))
#   
# })
# 
# test_that("returns list", {
#   # RUN MODEL SO OUTPUT IS CEAN AND STORE IN TEST DATA
#   
#   # diagnostics
#   diag <- .diag_test_helper("good")
#   expect_true(is.list(diag))
#   
# })
# 
# 
# test_that("model outputs with no issues are detected", {
#   # RUN MODEL SO OUTPUT IS CEAN AND STORE IN TEST DATA
#   
#   # diagnostics
#   diag <- .diag_test_helper("good")
#   expect_true(diag$divergence > 0)
#   
# })
# 
# test_that("model outputs with divergent transitions are detected", {
#   # RUN MODEL SO OUTPUT HAS DIVERGNET AND TREE TRANSIONS AND STORE IN TEST DATA
# })
# 
# test_that("model outputs with large tree depths are detected", {
#   # RUN MODEL SO OUTPUT HAS DIVERGNET AND TREE TRANSIONS AND STORE IN TEST DATA
# })
# 
# test_that("model outputs with small effect sample sizes are detected", {
#   # RUN MODEL SO OUTPUT HAS DIVERGNET AND TREE TRANSIONS AND STORE IN TEST DATA
# })
# 
# test_that("model outputs with high rhat are detected", {
#   # RUN MODEL SO OUTPUT HAS DIVERGNET AND TREE TRANSIONS AND STORE IN TEST DATA
# })
# 
# 

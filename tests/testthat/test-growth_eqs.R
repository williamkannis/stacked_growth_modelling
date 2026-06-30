
# .length2growth tests  --------------------------------------------------------

# output format
test_that("Returns numeric values length of input",{
  expect_true(is.numeric(.length2growth(20,"gz",34,.002,10)))
  input <- c(10,13,15)
  input_length <- length(input)
  g <- .length2growth(input,"gz",34,.002,10)
  expect_true(length(g) == input_length)
  expect_true(all(is.numeric(g)))
})

# coefficent values
test_that("throws error when impossible model coefficnets given",{
  # expect_error(.length2growth(-20,"gz",34,.002,10))
  expect_error(.length2growth(20,"gz",-34,.002,10))
  expect_error(.length2growth(20,"gz",34,-.002,10))
})

# model forms
test_that(
  "accepted model forms return results and inccorect model forms return error",{
    expect_no_error(
      sapply(c("gz","lg","vb"), function(mod) {
      .length2growth(20,mod,34,.002,10)
        }
      )
    )
    expect_error(.length2growth(20,"gf",34,.002,10))
})

# Values below asympotote
test_that("value below asympote return plausible",{
  mods <- c("vb","gz","lg")
  cond <- lapply(mods, function(mod){
    Linf <- 34
    input <- 1:(Linf-1)
    g <- .length2growth(input,mod,34,.002,10)
    list(
      postive = all(g >= 0),  # negative ages can idicate sizes expected before birth
      not_na = all(!is.na(g)),
      not_nan = all(!is.nan(g)),
      finite = all(is.finite(g))
    )
  })
  names(cond) <- mods
  cond_t <- purrr::transpose(cond)
  expect_all_true(unlist(cond_t$postive)) 
  expect_all_true(unlist(cond_t$not_na))
  expect_all_true(unlist(cond_t$not_nan))
  expect_all_true(unlist(cond_t$finite))
})

# Values above asynpote
test_that("Values above and at asymopote return zero",{
  g <-sapply(c("gz","lg","vb"),function(mod){
    linf <- 34
    input <- linf + 1
    .length2growth(input,"gz",linf,.002,10)
  }
  )
  expect_all_true(g == 0)
  
  g <-sapply(c("gz","lg","vb"),function(mod){
    linf <- 34
    input <- linf
    .length2growth(input,"gz",linf,.002,10)
  }
  )
  expect_all_true(g == 0)
})

# are correct values returned?
test_that("Model returns correct results",{
  
  # params
  Linf = 34
  g = 0.04
  inf = 10
  
  # input
  age <- seq(20, 100, by = 10)
  
  g <- sapply(c("vb","gz","lg"), function(mod){
    L <- .age2length(age,mod,Linf,g,inf)
    
    # Expected growth from the derivative of the age-length equation
    if(mod == "vb") growth_expected <- Linf * g * exp(-g * (age - inf))
    if(mod == "gz") growth_expected <- L * g * exp(-g * (age - inf))
    if(mod == "lg") growth_expected <- g * L * (1 - L / Linf)
    
    # Growth from your function
    growth_function <- .length2growth(
      input = L,
      g.mod = mod,
      Linf = Linf,
      g = g,
      inf = inf
    )
    
    all(round(growth_expected,10) == round(growth_function,10))
    
  }
  )
  expect_all_true(g)
  
})


# .length2age tests  --------------------------------------------------------


# coefficient values
test_that("throws error when impossible model coefficnets given",{
  # expect_error(.length2age(-20,"gz",34,.002,10))
  expect_error(.length2age(20,"gz",-34,.002,10))
  expect_error(.length2age(20,"gz",34,-.002,10))
})

# model forms
test_that(
  "accepted model forms return results and inccorect model forms return error",{
    expect_no_error(
      sapply(c("gz","lg","vb"), function(mod) {
        .length2age(20,mod,34,.002,10)
      }
      )
    )
    expect_error(.length2age(20,"gf",34,.002,10))
  })

# output format
test_that("Returns numeric values length of input",{
  expect_true(is.numeric(.length2age(20,"gz",34,.002,10)))
  input <- c(10,13,15)
  input_length <- length(input)
  g <- .length2age(input,"gz",34,.002,10)
  expect_true(length(g) == input_length)
  expect_true(all(is.numeric(g)))
})

# Values below asympotote
test_that("value below asympote return plausible",{
  mods <- c("vb","gz","lg")
  cond <- lapply(mods, function(mod){
    Linf <- 34
    input <- 1:(Linf-1)
    g <- .length2age(input,mod,34,.002,10)
    list(
      # postive = all(g >= 0),  # negative ages can idicate sizes expected before birth
      not_na = all(!is.na(g)),
      not_nan = all(!is.nan(g)),
      finite = all(is.finite(g))
    )
  })
  names(cond) <- mods
  cond_t <- purrr::transpose(cond)
  # expect_all_true(unlist(cond_t$postive)) 
  expect_all_true(unlist(cond_t$not_na))
  expect_all_true(unlist(cond_t$not_nan))
  expect_all_true(unlist(cond_t$finite))
})

# Values above asynpote
test_that("Values above and at asymopote return Inf",{
  g <-sapply(c("gz","lg","vb"),function(mod){
    linf <- 34
    input <- sapply(0:10, function(x) (linf + 1)+x)
    all(.length2age(input,"gz",linf,.002,10) == Inf)
  }
  )
  expect_all_true(g)
  
  g <-sapply(c("gz","lg","vb"),function(mod){
    linf <- 34
    input <- linf
    .length2age(input,"gz",linf,.002,10)
  }
  )
  expect_all_true(g == Inf)
})

# Negative inputs
# test_that("negative ages return non NA values?",{
#   mods <- c("vb","gz","lg")
#   cond <- lapply(mods, function(mod){
#     input <- -34:0
#     g <- .length2age(input,mod,34,.002,10)
#     list(
#       # postive = all(g >= 0),  # negative ages can idicate sizes expected before birth
#       not_na = all(!is.na(g)),
#       not_nan = all(!is.nan(g)),
#       finite = all(is.finite(g))
#     )
#   })
#   names(cond) <- mods
#   cond_t <- purrr::transpose(cond)
#   # expect_all_true(unlist(cond_t$postive)) 
#   expect_all_true(unlist(cond_t$not_na))
#   expect_all_true(unlist(cond_t$not_nan))
#   expect_all_true(unlist(cond_t$finite))
# })
## MAYBE CHANGE CODE TO PRODUCE AGES OF ZERO???

# are correct values returned?
test_that("Model returns correct results",{
  
  # params
  Linf = 34
  g = 0.04
  inf = 10
  
  # input
  real_a <- 0:200
  
  cond <- sapply(c("vb","gz","lg"), function(mod){
    l <- .age2length(real_a,mod,Linf,g,inf)
    fun_age <- .length2age(l,mod,Linf,g,inf)
    all(round(fun_age,10) == round(real_a,10))
  }
  )
  expect_all_true(cond)
  
})


# .age2length tests  --------------------------------------------------------


# coefficient values
test_that("throws error when impossible model coefficnets given",{
  # expect_error(.age2length(-20,"gz",34,.002,10))
  expect_error(.age2length(20,"gz",-34,.002,10))
  expect_error(.age2length(20,"gz",34,-.002,10))
})

# model forms
test_that(
  "accepted model forms return results and inccorect model forms return error",{
    expect_no_error(
      sapply(c("gz","lg","vb"), function(mod) {
        .age2length(20,mod,34,.002,10)
      }
      )
    )
    expect_error(.age2length(20,"gf",34,.002,10))
  })

# output format
test_that("Returns numeric values length of input",{
  expect_true(is.numeric(.age2length(20,"gz",34,.002,10)))
  input <- c(10,13,15)
  input_length <- length(input)
  g <- .age2length(input,"gz",34,.002,10)
  expect_true(length(g) == input_length)
  expect_true(all(is.numeric(g)))
})

# NA Values for postive inputs
test_that("return plausible values",{
  mods <- c("vb","gz","lg")
  cond <- lapply(mods, function(mod){
    input <- 0:360
    g <- .age2length(input,mod,34,.002,10)
    list(
      # postive = all(g >= 0),  # negative ages can idicate sizes expected before birth
      not_na = all(!is.na(g)),
      not_nan = all(!is.nan(g)),
      finite = all(is.finite(g))
    )
  })
  names(cond) <- mods
  cond_t <- purrr::transpose(cond)
  # expect_all_true(unlist(cond_t$postive)) 
  expect_all_true(unlist(cond_t$not_na))
  expect_all_true(unlist(cond_t$not_nan))
  expect_all_true(unlist(cond_t$finite))
})

# Negative inputs
test_that("negative ages return non NA values?",{
  mods <- c("vb","gz","lg")
  cond <- lapply(mods, function(mod){
    input <- -360:0
    g <- .age2length(input,mod,34,.002,10)
    list(
      # postive = all(g >= 0),  # negative ages can idicate sizes expected before birth
      not_na = all(!is.na(g)),
      not_nan = all(!is.nan(g)),
      finite = all(is.finite(g))
    )
  })
  names(cond) <- mods
  cond_t <- purrr::transpose(cond)
  # expect_all_true(unlist(cond_t$postive)) 
  expect_all_true(unlist(cond_t$not_na))
  expect_all_true(unlist(cond_t$not_nan))
  expect_all_true(unlist(cond_t$finite))
})


# are correct values returned?
test_that("Model returns correct results",{
  
  # params
  Linf = 34
  g = 0.04
  inf = 10
  
  # input
  l <- seq(1, Linf-1)

  cond <- sapply(c("vb","gz","lg"), function(mod){
    age <- .length2age(l,mod,Linf,g,inf)
    returned_l <- .age2length(age,mod,Linf,g,inf)
    all(round(l,10) == round(returned_l,10))
  }
  )
  expect_all_true(cond)
  
})


# .age2growth  -----------------------------------------------------------------

# output format
test_that("Returns numeric values length of input",{
  expect_true(is.numeric(.age2growth(20,"gz",34,.002,10)))
  input <- c(10,13,15)
  input_length <- length(input)
  g <- .age2growth(input,"gz",34,.002,10)
  expect_true(length(g) == input_length)
  expect_true(all(is.numeric(g)))
})


# .age2interval_growth  -------------------------------------------------------
wt_df_helper <- function(){
  data.frame(
    a =-4.782,
    b=3.042,
    c= 1
      )
}
# output format
test_that("Returns numeric values length of input",{
  wt.df <- wt_df_helper()
  interval = 30
  expect_true(
    is.numeric(
      .age2interval_growth(
        20,
        "gz",
        34,
        .002,
        10,
        interval=interval,
        wt.df=wt.df
        )
      )
    )
  input <- c(10,13,15)
  input_length <- length(input)
  g <- .age2interval_growth(
    input,
    "gz",
    34,
    .002,
    10,
    interval=interval,
    wt.df=wt.df
    )
  expect_true(length(g) == input_length)
  expect_true(all(is.numeric(g)))
})


# .length2interval_growth  -----------------------------------------------------

# output format
test_that("Returns numeric values length of input",{
  wt.df <- wt_df_helper()
  interval = 30
  expect_true(
    is.numeric(
      .length2interval_growth(
        20,
        "gz",
        34,
        .002,
        10,
        interval=interval,
        wt.df=wt.df
        )
    )
  )
  input <- c(10,13,15)
  input_length <- length(input)
  g <- .length2interval_growth(
    input,
    "gz",
    34,
    .002,
    10,
    interval=interval,
    wt.df=wt.df
    )
  expect_true(length(g) == input_length)
  expect_true(all(is.numeric(g)))
})


# .length_forcast  -------------------------------------------------------------

# interval tests
## THIS CAN POSSIBLY BE DISBALED, BUT BLOCK THIS FOR NOW
test_that("negative interval throws error", {
  expect_error(.length_forcast(20,interval = -3, "gz",34,.002,10))
})

# output format
test_that("Returns numeric values length of input",{
  expect_true(is.numeric(.length_forcast(20,30,"gz",34,.002,10)))
  input <- c(10,13,15)
  input_length <- length(input)
  g <- .length_forcast(input,30,"gz",34,.002,10)
  expect_true(length(g) == input_length)
  expect_true(all(is.numeric(g)))
})

# values below asympote
test_that("all growth forms return values for inputs below asympote",{
  
  # params
  Linf = 34
  g = 0.04
  inf = 10
  
  #
  input = 1:(Linf-1)
  
  # input
  
  cond <- lapply(c("vb","gz","lg"), function(mod){
    .length_forcast(input,20,mod,Linf,g,inf)
    
    list(
      postive = all(g >= 0),  # negative ages can idicate sizes expected before birth
      not_na = all(!is.na(g)),
      not_nan = all(!is.nan(g)),
      finite = all(is.finite(g))
    )
  }
  )
  
  names(cond) <- mods
  cond_t <- purrr::transpose(cond)
  expect_all_true(unlist(cond_t$postive)) 
  expect_all_true(unlist(cond_t$not_na))
  expect_all_true(unlist(cond_t$not_nan))
  expect_all_true(unlist(cond_t$finite))
  })

# values above aympote
test_that("inputs at or above asympotote returns input",{
  cond <- sapply(c("vb","gz","vb"),function (mod){
    Linf = 34
    input <- sapply(0:10,function(x) Linf+x)
  
    output <- .length_forcast(input,30,mod,34,.002,10)
    all(input == output)
  }
  )
  expect_all_true(cond)
})








# fit_growth tests  -------------------------------------------------------------



# .prep_stan_data tests  -------------------------------------------------------



## all data structure

## throws errors ###
# Missing length data from age.df (cant be any)
# Missing age data from age.df (cant be any)
# missing length data from len.df (cant be all)


## Linear predictors  ##

## Cat predictors

## linear predictions  ##

## Mean length ##

## Bridge dataframe

test_that(".prep_stan_data returns a list",{
  df <- .data_prep_test_helper()
  prep <- .prep_stan_data(
    sp = "a", 
    age.df = df, 
    sample.groups = c("date","site"),   
    fixed.effect = "random"
  )
  expect_true(is.list(prep))
})

## incorrect input errors - age.df  ##

test_that(
  ".prep_stan_data returns an error with missing species column from age.df",{
    df <- .data_prep_test_helper() %>% dplyr::select(-species)
    expect_error(.prep_stan_data(
      sp = "a", 
      age.df = df, 
      sample.groups = c("date","site")  
    ))
  })
test_that(
  ".prep_stan_data returns an error with missing species value from age.df",{
    df <- .data_prep_test_helper()
    expect_error(.prep_stan_data(
      sp = "b", 
      len.df = len_df,
      age.df = df, 
      sample.groups = c("date","site")  
    ))
  })
test_that(
  ".prep_stan_data returns an error when age.df is missing length column",{
    df <- .data_prep_test_helper() %>% dplyr::select(-length)
    expect_error(.prep_stan_data(
      sp = "a", 
      age.df = df, 
      sample.groups = c("date","site")  
    ))
    
  })
test_that(
  ".prep_stan_data returns an error when age.df is missing age column",{
    df <- .data_prep_test_helper() %>% dplyr::select(-age)
    expect_error(.prep_stan_data(
      sp = "a", 
      age.df = df, 
      sample.groups = c("date","site")  
    ))
    
  })


## ##

test_that(
  ".prep_stan_data returns an error with missing species column from len.df",{
    df <- .data_prep_test_helper() 
    len_df <- 
      data.frame(
        speciexs = rep("b",1000),
        length = rnorm(1000))
    expect_error(.prep_stan_data(
      sp = "a", 
      age.df = df, 
      len.df = len_df,
      sample.groups = c("date","site")  
    ))
  })

test_that(
  ".prep_stan_data returns an error with missing species value from len.df",{
    df <- .data_prep_test_helper()
    len_df <- 
      data.frame(
        species = rep("b",1000),
        length = rnorm(1000))
    expect_error(.prep_stan_data(
      sp = "a", 
      age.df = df, 
      len.df = len_df,
      sample.groups = c("date","site")  
    ))
    
  })
test_that(
  ".prep_stan_data returns an error with missing length column  from len.df",{
    df <- .data_prep_test_helper()
    len_df <- 
      data.frame(
        species = rep("a",1000),
        lengths = rnorm(1000))
    expect_error(.prep_stan_data(
      sp = "a", 
      len.df = len_df,
      age.df = df, 
      sample.groups = c("date","site")  
    ))
  })

## Errors - groupings  ###

test_that(
  ".prep_stan_data returns an error when no groupings are supplied",{
    df <- .data_prep_test_helper()
    expect_error(.prep_stan_data(
      sp = "a", 
      age.df = df, 
      fixed.effect = "random"
    ))
  }
)
test_that(
  ".prep_stan_data returns an error when incorrect groupings are supplied",{
  df <- .data_prep_test_helper()
  expect_error(.prep_stan_data(
    sp = "a", 
    age.df = df, 
    sample.groups = c("date","site","plot"),   
    fixed.effect = "random"
  ))
})

## Errors - effect sturucture  ##
test_that(
  ".prep_stan_data returns an error with invalid effect strucuture",{
    df <- .data_prep_test_helper()
    expect_error(.prep_stan_data(
      sp = "a", 
      age.df = df, 
      sample.groups = c("date","site"),   
      fixed.effect = "no effects"
    ))
    
  })

test_that(
  ".prep_stan_data returns an error when incorrect predictors are supplied",{
    df <- .data_prep_test_helper()
    expect_error(.prep_stan_data(
      sp = "a", 
      age.df = df, 
      sample.groups = c("date","site"),   
      fixed.effect = "linear",
      predictors = "depth"
    ))
  })
test_that(
  ".prep_stan_data returns an error when predictors are missing but expected",{
    df <- .data_prep_test_helper()
    expect_error(.prep_stan_data(
      sp = "a", 
      age.df = df, 
      sample.groups = c("date","site"),   
      fixed.effect = "linear"
    ))
  })

test_that(
  ".prep_stan_data returns an error when incorrect categories are supplied",{
    df <- .data_prep_test_helper()
    expect_error(.prep_stan_data(
      sp = "a", 
      age.df = df, 
      sample.groups = c("date","site"),   
      fixed.effect = "categorical",
      category  = "state"
    ))
  })
test_that(
  ".prep_stan_data returns an error when categories are missing but expected",{
    df <- .data_prep_test_helper()
    expect_error(.prep_stan_data(
      sp = "a", 
      age.df = df, 
      sample.groups = c("date","site"),   
      fixed.effect = "categorical"
    ))
  })

test_that(
  paste0(
    ".prep_stan_data returns an error when categories are provided but not ",
    "needed"
  ),{
    df <- .data_prep_test_helper()
    expect_error(.prep_stan_data(
      sp = "a", 
      age.df = df, 
      sample.groups = c("date","site"),   
      fixed.effect = "linear",
      category  = "cat_name"
    ))
  })

test_that(
  paste0(
    ".prep_stan_data returns an error when predictors are provided but not ",
    "needed"
  )
  ,{
    df <- .data_prep_test_helper()
    expect_error(.prep_stan_data(
      sp = "a", 
      age.df = df, 
      sample.groups = c("date","site"),   
      fixed.effect = "linear",
      category  = "cat_name"
    ))
  })

test_that(
  paste0(
    ".prep_stan_data returns an error when linear predictions are specified ",
    "but not possible due to random effect strucutre"
    ),{
    df <- .data_prep_test_helper()
    expect_error(.prep_stan_data(
      sp = "a", 
      age.df = df, 
      sample.groups = c("date","site"),   
      linear.predictions = T
    ))
  })

test_that(
  paste0(
    ".prep_stan_data returns an error when linear predictions are specified ",
    "but not possible due to cateorical effect strucutre"
  ),{
    df <- .data_prep_test_helper()
    expect_error(.prep_stan_data(
      sp = "a", 
      age.df = df, 
      sample.groups = c("date","site"),   
      fixed.effect = "categorical",
      linear.predictions = T
    ))
  })

## Random effect test  ##


test_that("number of groups in .prep_stan_data output equal to input data",{
  df <- .data_prep_test_helper()
  prep <- .prep_stan_data(
    sp = "a", 
    age.df = df, 
    sample.groups = c("date","site"),   
    fixed.effect = "random"
  )
  in_n_group <- df %>% dplyr::distinct(date,site) %>% nrow()
  out_n_group <- length(unique(prep$stan_data$ID))
  expect_true(in_n_group == out_n_group)
})

test_that(".prep_stan_data outputs have correct number of items",{
  df <- .data_prep_test_helper()
  prep <- .prep_stan_data(
    sp = "a", 
    age.df = df, 
    sample.groups = c("date","site"),   
    fixed.effect = "random"
  )
  prep_length <- lapply(prep$stan_data,length)
  in_length <- nrow(df)
  expect_true(all(
    in_length == prep_length$LENGTH,
    in_length == prep_length$AGE,
    in_length == prep_length$ID,
    in_length == prep$stan_data$N))
})

# linear effect tests
# prep <- .prep_stan_data(
#   sp, 
#   age.df, 
#   len.df = NULL,
#   NU = 0,
#   sample.groups,   
#   fixed.effect = NULL,
#   category = NULL, 
#   predictors = NULL,  
#   scale = T,
#   linear.predictions = F,  
#   pred.len = 100)

test_that(".prep_stan_data returns predicotr data when specified",{
  df <- .data_prep_test_helper()
  prep <- .prep_stan_data(
    sp = "a", 
    age.df = df, 
    sample.groups = c("date","site"),   
    fixed.effect = "random"
  )
  prep_length <- lapply(prep$stan_data,length)
  in_length <- nrow(df)
  expect_true(all(
    in_length == prep_length$LENGTH,
    in_length == prep_length$AGE,
    in_length == prep_length$ID,
    in_length == prep$stan_data$N))
})

# .extract_stan_file tests  ----------------------------------------------------

test_that("stan_file returns an vb existing linear file", {
  
  path <- .extract_stan_file("vb", "linear")
  
  expect_true(file.exists(path))
  
})

test_that("stan_file returns an existing categorical file", {
  
  path <- .extract_stan_file("vb", "categorical")
  
  expect_true(file.exists(path))
  
})

test_that("stan_file returns an existing random effect file", {
  
  path <- .extract_stan_file("vb", "random")
  
  expect_true(file.exists(path))
  
})

test_that("invalid model throws an error", {
  
  expect_error(
    .extract_stan_file("fake_model", "linear")
  )
  
})

test_that("invalid effect structure throws an error", {
  
  expect_error(
    .extract_stan_file("vb", "fake_effect")
  )
  
})



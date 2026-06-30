
# Create length and age input data frame ---------------------------------------
.data_prep_test_helper <- function(){
  
  data.frame(
    species = rep("a",300),
    length = runif(300,5,40),
    age = runif(300,0,360),
    date = rep(rep(1:3,10),10),
    site = rep(rep(1:10, each = 3),10),
    cat_name = sample(1:3,300,T),
    pred1 = rnorm(300,0,1),
    pred2 = rnorm(300,0,1),
    pred3 = rnorm(300,0,1)
  )
}





# Create testing wt conversion data.frames  ------------------------------------
wt_df_helper <- function(x){
  
  if(x == "good"){
    df <- data.frame(
      a =-4.782,
      b=3.042,
      c= 1
    )
  }
  
  if(x == "bad"){
    df <- data.frame(
      x =-4.782,
      y=3.042,
      z= 1
    )
  }
  
  if(x == "missing_a"){
    df <- data.frame(
      x =-4.782,
      b=3.042,
      c= 1
    )
  }
  
  if(x == "missing_b"){
    df <- data.frame(
      a =-4.782,
      x=3.042,
      c= 1
    )
  }
  
  if(x == "missing_c"){
    df <- data.frame(
      a =-4.782,
      b=3.042,
      x= 1
    )
  }
  if(x == "extra_rows"){
    df <- data.frame(
      a =c(-4.782,3,5),
      b=c(3.042,3,5),
      c= c(1,1,1)
    )
  }
  df
}

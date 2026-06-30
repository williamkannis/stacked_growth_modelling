# Create testing wt conversion data.frames
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

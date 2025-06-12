
# Function that performs Random Walk on a normalized contact matrix # 

RW <- function(normmatrix, step=3){
  Q <- normmatrix
  if(step<=0 | !(step%%1==0)){
    stop("Step must be a positive integer!")
  }else if(step==1){
    return(Q)
  }else{
    for(k in 1:(step-1)){
      Q <- Q%*%normmatrix
    }
    return(Q)
  }
}



# Function that performs RWR on a normalized contact matrix # 

RWR <- function(normmatrix, alpha=0.1){
  if(alpha<=0 | alpha>=1){
    stop("Restart probability must be between 0 and 1!")
  }else{
    n <- dim(normmatrix)[1]
    Q <- diag(n)
    for(i in 1:300){
      Q <- ((1-alpha)*Q%*%normmatrix+alpha*diag(n))
    }
    return(Q) 
  }
}


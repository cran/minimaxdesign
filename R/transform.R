is.whole <- function(a) {
  (is.numeric(a) && floor(a)==a) ||
    (is.complex(a) && floor(Re(a)) == Re(a) && floor(Im(a)) == Im(a))
}

CtoA <- function(D, by=ifelse(ncol(D)>2,1e-3,-1), num_proc=parallel:::detectCores()){
  ret <- CtoAA(D,by,num_proc)
  return(ret)
}

CtoB <- function(D, by=ifelse(ncol(D)>2,1e-3,-1), num_proc=parallel:::detectCores()){
  pp <- ncol(D)
  if (pp > 2){
    ret <- CtoBp(D,by,num_proc)
  }else{
    ret <- CtoB2(D,by,num_proc)
  }
  return(ret)
}

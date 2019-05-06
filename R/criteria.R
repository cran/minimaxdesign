########################################################################
# Minimax distance
########################################################################

mMdist=function(D,neval=1e5,method="lattice",region="hypercube",const=NA){

  bd <- c(0,1) #Default limits are (0,1)
  p <- ncol(D)

  # Define transformation function
  switch(region,
         hypercube = {
           tf <- function(D,by,num_proc){return(D)}
           checkBounds <- function(D){ #function for checking bound
             return(D)
           }
         },
         simplex = {
           tf <- CtoA;
           checkBounds <- function(D){ #function for checking bound
             return(D)
           }
         },
         ball = {
           tf <- CtoB;
           bd <- c(-1,1);
           checkBounds <- function(D){ #function for checking bound
             good.idx <- which(rowSums(D^2)<=1)
             return(D[good.idx,])
           }
         },
         custom = {
           tf <- function(D,by,num_proc){ #randomly sample until enough points
             num_pts <- nrow(D)
             samp <- matrix(NA,nrow=num_pts,ncol=p)
             cur_pts <- 0
             while (cur_pts < num_pts){
               xx <- stats::runif(p)
               if (const(xx)){
                 samp[cur_pts+1,] <- xx
                 cur_pts <- cur_pts + 1
               }
             }
             return(samp)
           }
           checkBounds <- function(D){ #function for checking bound
             good.idx <- apply(D,1,const)
             return(D[good.idx,])
           }
         }
         )
  regionby=ifelse(ncol(D)>2,1e-3,-1)
  by <- regionby

  # Generating evaluation points
  switch(method,
         lattice = {
           m <- nrow(D)
           d <- ncol(D)
           M <- tf(as.matrix(Lattice(max(conf.design::primes(neval)),d)),by,parallel::detectCores())
           M <- rbind(M,tf(gtools::permutations(3,p,c(0.0,0.5,1.0),repeats.allowed=TRUE))) #add 3^p design
           M <- checkBounds(M)
         },
         sobol = {
           m <- nrow(D)
           d <- ncol(D)
           M <- tf(randtoolbox::sobol(neval,d),by,parallel::detectCores())
           M <- rbind(M,tf(gtools::permutations(3,p,c(0.0,0.5,1.0),repeats.allowed=TRUE))) #add 3^p design
           M <- checkBounds(M)
         }
  )

  # Compute minimax distance
  mM=as.matrix(pdist::pdist(M,D))
  ind=apply(mM,1,which.min)
  distD=numeric(m)
  for(i in 1:m){
    distD[i]=max(mM[ind==i,i])
  }

  maxD <- max(distD)
  maxI <- which.max(distD)
  ptI <- which(ind==maxI)
  ptI <- ptI[which.max(mM[ind==maxI,maxI])]
  return(list(dist=max(distD),dist.vec=distD,far.pt=M[ptI,],far.despt=D[maxI,]))
}

# ########################################################################
# # Maximin distance
# ########################################################################
#
# Mmdist=function(D){
#   distD <- max(stats::dist(D))
#   return(distD)
# }

########################################################################
# Computing shifted lattices
########################################################################
Lattice=function(n,p,shift=FALSE){
  #n: number of points in the lattice rule, scalar, must be prime
  #n = 101
  #s_max: number of dimensions, scalar
  s_max = p
  #omega: function for the varying kernel part of the shift-invariant kernel function (assumed to be symmetric)
  omega=function(x)
    2*pi^2*(x^2-x+1/6)
  #gamma:gamma parameters for AC weighting per dimension, vector of length s_max
  gamma = rep(1,s_max) / s_max
  #beta: beta parameters for DC weighting per dimension, vector of length s_max
  beta = rep(1,s_max)
  ############################################

  #Define fuction "powmod"
  #Calculate x^a (mod n) using the Russian Peasant method
  powmod=function(x,a,n){
    y=1
    u=x
    while(a>0){
      if(a%%2==1)
        y=(y*u)%%n
      u=(u*u)%%n
      a = floor(a / 2)
    }
    return(y)
  }

  #Define fuction "generatorp"
  #Return a generator for the cyclic group multiplication modulo p
  #Also called a primitive root of p
  #input:
  #p:a prime number
  #output:
  #g:the generator

  #Get all irreducible factors of an integer
  generatorp <- function(p){
    primef=as.numeric(unique(conf.design::factorize(p-1)))
    g=2
    i=1
    while(i<=length(primef)){
      if(powmod(g,(p-1)/primef[i],p)==1){
        g=g+1
        i=0
      }
      i=i+1
    }
    return(g)
  }

  ############################################
  ###Generate the lattice rule###

  #s_max=p
  m = (n-1)/2          #assume the omega function symmetric around 1/2
  E2 = rep(0,m)       # the vector $\tilde{\vec{E}}^2$ in the text
  cumbeta = cumprod(beta)

  g = generatorp(n) #generator $g$ for $\{1, 2, \ldots, n-1\}$
  perm=rep(0,m) #permutation formed by positive powers of $g$
  perm[1]=1
  for(j in 1:(m-1))
    perm[j+1]=(perm[j]*g)%%n
  perm = apply(cbind(n - perm, perm),1,min) #map everything back to $[1, n/2)$
  psi= apply(matrix(perm/n,ncol=1),1,omega) #the vector $\vec{\psi}'$
  psi0 = omega(0)       # zero index: $\psi(0)$
  fft_psi = stats::fft(psi)

  #z:generating vector of the lattice rule, vector of length s_max
  #e2:optimal square error per dimension (= one for each iteration), vector of length s_max
  z=rep(NA,s_max)
  e2=rep(NA,s_max)
  q=rep(1,m) #permuted product vector $\vec{q}'$ (without zero index)
  q0 = 1     #zero index of permuted product vector: $q(0)$

  for(s in 1:s_max){
    #step 2a: circulant matrix-vector multiplication
    E2=stats::fft(fft_psi*stats::fft(q),inverse=T)/m
    E2 = Re(E2)#remove imaginary rounding errors
    #step 2b: choose $w_s$ and $z_s$ which give minimal value
    min_E2=min(E2)
    w=which.min(E2)#pick index of minimal value
    if(s==1){
      w=1
      noise = abs(E2[1] - min_E2)
    }
    z[s]=perm[w]
    #extra: we want to know the exact value of the worst-case error
    e2[s] = -cumbeta[s] + ( beta[s] * (q0 + 2*sum(q)) +
                              gamma[s] * (psi0*q0 + 2*min_E2) ) / n
    #step 2c: update $\vec{q}$
    q = (beta[s] + gamma[s] * psi[c(seq(w,1,-1),seq(m,(w+1),-1))])* q
    q0 = (beta[s] + gamma[s] * psi0) * q0
  }

  #Output
  #z:generating vector of the lattice rule, vector of length s_max
  #e2:optimal square error per dimension (= one for each iteration), vector of length s_max
  #output=data.frame(z,e2,sqrt(e2))



  tseq=function(gamma)
  {
    x=((((1:n))*gamma)%%n)/n
    return(x)
  }
  D=matrix(0,nrow=n,ncol=p)
  for(i in 1:p)
    D[,i]=cbind(tseq(z[i]))
  l=apply(D,2,min)
  u=apply(D,2,max)
  D=.5/n+(1-1/n)*(D-rep(1,n)%*%t(l))/(rep(1,n)%*%t(u-l))

  if (shift){ # random shift
    D <- apply(D,2,function(xx){ ((xx+stats::runif(1))%%1) })
  }

  return(D)
}

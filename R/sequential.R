minimax.seq <- function(D,Nseq,q=10,region="hypercube",const=NA,subsp=NA,wt.subsp=FALSE,
                    params_pso=list(w=0.72,c1=1.49,c2=1.49),
                    npart=5,nclust=1e5,neval=nclust,
                    itmax_pso=50,itmax_pp=100,itmax_inn=1e4,jit=0.1/sqrt(Nseq+nrow(D))){

  if (any(is.na(subsp))){
    p <- ncol(D)
  }else{ #need to define projected domain on active subspace
    p <- ncol(D)
    ini.samp <- matrix(runif(nclust*nrow(subsp)),nrow=nclust,ncol=nrow(subsp))
    ini.samp <- t(t(subsp)%*%t(ini.samp))
    idx.hull <- geometry::convhulln(ini.samp) #indices of convex hull
    dim.min <- apply(ini.samp[c(idx.hull),],2,min)
    dim.max <- apply(ini.samp[c(idx.hull),],2,max)
  }
  N <- nrow(D)

  #Setting transformation type
  regionby=ifelse(p>2,1e-3,-1) #for simplex transformation
  bd <- c(0,1) #Default limits are (0,1)

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
         },
         subspace = {
           tf <- function(D,by,num_proc){ #randomly sample until enough points
             num_pts <- nrow(D)
             if (!wt.subsp){
               samp <- matrix(NA,nrow=num_pts,ncol=p)
               cur_pts <- 0
               while (cur_pts < num_pts){ #rejection sampling
                 xx <- stats::runif(p) * (dim.max-dim.min) + dim.min
                 if ( geometry::inhulln(idx.hull,matrix(xx,nrow=1)) ){
                   samp[cur_pts+1,] <- xx
                   cur_pts <- cur_pts + 1
                 }
               }
               return(samp)
             }else{
               plot.pts <- matrix(runif(num_pts*nrow(subsp)),ncol=nrow(subsp))
               ret <- t(t(subsp)%*%t(plot.pts))
               return(ret)
             }
           }
           checkBounds <- function(D){ #function for checking bound
             return(D)
           }
         }
         )
  by <- regionby

  print("Generating clustering points ...")
  #Generate point
  clust_pts = tf(as.matrix(Lattice(max(conf.design::primes(nclust)),p))) #clustering points
  clust_pts <- checkBounds(clust_pts)

  #Generate eval_pts
  eval_pts = tf(as.matrix(Lattice(max(conf.design::primes(neval)),p,shift=TRUE))) #minimax approximation points for post-processing
  eval_pts = rbind(eval_pts,tf(gtools::permutations(3,p,c(0.0,0.5,1.0),repeats.allowed=TRUE))) #add 3^p design
  eval_pts <- checkBounds(eval_pts)

  # Generate initial center particles

  #automatically set initialization flag
  by.clust <- 1e-4
  cluster_center = rbind(D,tf(randtoolbox::sobol(Nseq,p),by.clust,parallel::detectCores())) #initial cluster centers
  for (i in 1:(npart-1)){ #add scrambled randtoolbox::sobol sets
    cluster_center = cbind(cluster_center,
                           rbind(D,tf(randtoolbox::sobol(Nseq,p,init=FALSE),by.clust,parallel::detectCores())) )
  }
  fix_ind <- c(rep(1,nrow(D)),rep(0,Nseq))

  # mMc-PSO: function coded in C++
  t1 = Sys.time()
  D = kmeanspso(clust_pts, eval_pts, cluster_center, q, 2.0,
                params_pso$w,params_pso$c1,params_pso$c2,
                2*npart, itmax_pso, itmax_pp, 1000, 1000, 1e-4, 1e-4,
                1e-4, itmax_inn,
                parallel::detectCores(), jit, bd[1], bd[2], fix_ind)
  fitpre <- D$gbes_centers
  fitpost <- fitpre;

  t2 = Sys.time()
  tm <- difftime(t2,t1)
  # print(5)

  # D$gbes_centers <- fitpost
  # D$time <- tm
  # return(D$gbes_centers)
  return(D$gbes_centers[(N+1):(N+Nseq),])
}

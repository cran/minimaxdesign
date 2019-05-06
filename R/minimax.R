# #set compiler flags for openMP (parallel processing in C++)
# Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
# Sys.setenv("PKG_LIBS"="-fopenmp")

########################################################################
# Main function for mMcPSO
########################################################################
minimax <- function(N,p,q=10,region="hypercube",ini=NA,const=NA,
                  params_pso=list(w=0.72,c1=1.49,c2=1.49),
                  npart=5,nclust=1e5,neval=nclust,
                  itmax_pso=50,itmax_pp=100,itmax_inn=1e4,jit=0.1/sqrt(N)){

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
         })
  by <- regionby

  print("Generating clustering points ...")
  # if (is.na(clust_pts)){
    #Generate point
    # point = tf(randtoolbox::sobol(nclust,p),by,parallel::detectCores()) #clustering points
    clust_pts = tf(as.matrix(Lattice(max(conf.design::primes(nclust)),p))) #clustering points
    clust_pts <- checkBounds(clust_pts)
  # }

  # if (is.na(eval_pts)){
    #Generate eval_pts
    # eval_pts = tf(randtoolbox::sobol(neval,p),by,parallel::detectCores()) #minimax approximation points for post-processing
    eval_pts = tf(as.matrix(Lattice(max(conf.design::primes(neval)),p,shift=TRUE))) #minimax approximation points for post-processing
    eval_pts = rbind(eval_pts,tf(gtools::permutations(3,p,c(0.0,0.5,1.0),repeats.allowed=TRUE))) #add 3^p design
    eval_pts <- checkBounds(eval_pts)
  # }

  #   offset = seq(0,1-1/part_num,1/part_num)

  # Generate initial center particles

  #automatically set initialization flag
  by.clust <- 1e-4
  if (is.na(ini)){
    if (region=="hypercube"){
      if ( (is.whole(log(N)/log(2))) && (N <= 2^p) ){ # ff initialization
        ini <- "ff"
      }else if( (N>2^p) && (N<=1.25*2^p) ){
        ini <- "ff"
      }else{
        ini <- "sobol"
      }
    }else{
      ini <- "sobol"
    }
  }
  # Initialize
  if (ini=="sobol"){ # initialize via randomized sobol'
    cluster_center = tf(randtoolbox::sobol(N,p),by.clust,parallel::detectCores()) #initial cluster centers
    for (i in 1:(npart-1)){ #add scrambled randtoolbox::sobol sets
      cluster_center = cbind(cluster_center,tf(randtoolbox::sobol(N,p,init=FALSE),by.clust,parallel::detectCores()))
    }
  }else if (ini=="ff"){ # initialize via 2^(k-p) fractional factorial
    if ( N >= 2*(2^p) ){
      stop("Design size too large for FF initialization!")
    }else{
      pwf = floor(log(N)/log(2));
      cur_center = tf(0.25*DoE.base::desnum(FrF2::FrF2(2^pwf,p))+0.5,by.clust,parallel::detectCores()) #initial cluster centers
      cur_center = rbind(cur_center,tf(jitter(randtoolbox::sobol(N-2^pwf,p)),by.clust,parallel::detectCores()))
      cluster_center = cur_center
      for (i in 1:(npart-1)){
        cur_center = tf(0.25*DoE.base::desnum(FrF2::FrF2(2^pwf,p))+0.5,by.clust,parallel::detectCores()) #initial cluster centers
        cur_center = rbind(cur_center,tf(jitter(randtoolbox::sobol(N-2^pwf,p,init=FALSE)),by.clust,parallel::detectCores()))
        cluster_center = cbind(cluster_center,cur_center)
      }
    }
  }

  # print("4")
  # mMc-PSO: function coded in C++
  t1 = Sys.time()
  D = kmeanspso(clust_pts, eval_pts, cluster_center, q, 2.0,
                params_pso$w,params_pso$c1,params_pso$c2,
                2*npart, itmax_pso, itmax_pp, 1000, 1000, 1e-4, 1e-4,
                1e-4, itmax_inn,
                parallel::detectCores(), jit, bd[1], bd[2])
  fitpre <- D$gbes_centers
  fitpost <- fitpre;

  t2 = Sys.time()
  tm <- difftime(t2,t1)
  # print(5)

  # D$gbes_centers <- fitpost
  # D$time <- tm
  return(D$gbes_centers)
  # return(eval_pts)
}

########################################################################
# Functions for minimax projection designs
########################################################################

#Main function for generating minimax projection designs
miniMaxPro <- function(N,p,mMdes=NA, mMtol=1e-3*p,
                       neval=1e5, itmax_refine=100, ...){
  #Generate minimax design from mMc-PSO
  if (is.na(mMdes)){
    mMdes <- minimax(N,p,...)
  }
  # if (is.na(refine_pts)){
    # refine_pts <- randtoolbox::sobol(refine_num,p)
  # }
  eval_pts = as.matrix(Lattice(max(conf.design::primes(neval)),p,shift=TRUE)) #minimax approximation points for post-processing
  eval_pts = rbind(eval_pts,gtools::permutations(3,p,c(0.0,0.5,1.0),repeats.allowed=TRUE)) #add 3^p design

  #Refinement
  ret <- refine(mMdes,eval_pts=eval_pts,pp_itmax=itmax_refine,mM_tol=mMtol)

  return( list(minimax=mMdes,miniMaxPro=ret) )

}

#Refinement step for miniMaxPro designs
refine <- function(mMdes,eval_pts=NA,
                   pp_itmax=50,mM_tol=1e-3*(ncol(mMdes)),jitter.scl=5){

  mp=function(xx)
  {
    #Blockwise MaxPro objective
    dif <- DD - rep(1,dim(DD)[1])%*%t(xx)
    d=apply(dif, 1, prod)
    val=sum(1/(d^2+10^(-10))) #For numerical stability
    lfn <- log(val)

    #... gradient
    B <- (t(1/dif)%*%matrix(d,ncol=1))
    G <- 2*B/val

    return(list("objective"=lfn,"gradient"=G))
  }

  con=function(xx)
  {
    #Constraint function
    return(list("constraints"= sum((xx-cur_pt)^2) - toler^2,
                "jacobian"= 2*(xx-cur_pt) ))
  }

  #Compute minimax criterion
  N <- nrow(mMdes)
  p <- ncol(mMdes)
  orig_vec <- mMcrit_allpts( mMdes, eval_pts )
  orig_mM <- max(mMcrit_allpts( mMdes, eval_pts )) #Current minimax distance

  #Perform MaxPro adjustment
  numCores <- parallel::detectCores()
  cluster = parallel::makeCluster(numCores,type="SOCK")
  doParallel::registerDoParallel(cluster)

  cur_des <- mMdes
  itcur <- 0

  print(paste0("MaxPro refinement ..."))
  while (itcur < pp_itmax){
    if(itcur>1) printBar(itcur/pp_itmax)

    #Blockwise MaxPro for each design point
    dst_vec <- mMcrit_allpts( cur_des, eval_pts ) #First compute the minimum distance for each design point
    dst_vec <- pmax(orig_mM + mM_tol - dst_vec,0)
    mM_ind <- which(dst_vec==0)
    prev_des <- cur_des

    indices <- gtools::permute(setdiff(1:N,mM_ind))
    flg_cont <- F

    while (!flg_cont){
      biglist <- foreach::'%dopar%'(foreach::foreach (i = indices, .packages="nloptr"),
      # biglist <- foreach::'%do%'(foreach::foreach (i = indices, .packages="nloptr"),
      {
        toler <- dst_vec[i]/jitter.scl # Ensures the jittered initial points do not exceed minimax
        # print(toler)
        cur_pt <- cur_des[i,] #Current point
        DD <- cur_des[-i,] #Design without current point

        #         ini_pt <- pmin(pmax(cur_des[i,]+toler/sqrt(p)*stats::runif(p,min=-1),0),1)
        ini_pt <- pmin(pmax(cur_pt+toler*stats::runif(p,min=-1),0),1)
        # ini_pt <- cur_pt
        tryCatch({a=nloptr::nloptr( x0=ini_pt,
                            mp,lb=rep(0,p),ub=rep(1,p),eval_g_ineq = con,
                            opts = list("algorithm"="NLOPT_LD_MMA","maxeval"=100))

                  list("obj"=a$objective,"cur_sol"=a$sol,"prev_obj"=mp(cur_des[i,])$objective)}
                 ,error = function(err){
                   list("obj"=mp(cur_des[i,])$objective,"cur_sol"=cur_des[i,],"prev_obj"=mp(cur_des[i,])$objective)
                 })

      })

      prev_obj <- rep(0,length(indices))
      cur_obj <- rep(0,length(indices))
      solvec <- vector("list",length(indices))
      for (i in 1:length(indices)){
        #If there is an improvement, then change. Otherwise, stick with original point
        if (biglist[[i]]$obj < biglist[[i]]$prev_obj){
          cur_des[indices[i],] <- biglist[[i]]$cur_sol
        }
      }
      flg_cont=T

#       #Check if minimax distance is violated
#       cur_dst_vec <- mMcrit_allpts( cur_des, eval_pts )
#       if (max(cur_dst_vec)-orig_mM <= mM_tol){
#         #... terminate loop
#         flg_cont = T
#         scl <- jitter.scl
#       }else{
#         cur_des <- prev_des
#         scl <- scl + 5
# #         print(paste0("Decreasing jitter: ", scl))
#
#         #If scl hits maxscl, terminate algorithm
#         if (scl >= maxscl){
#           cur_des <- prev_des
#           flg_cont = T
#           itcur <- pp_itmax
#         }
#
#       }

    }

#     plot(mMdes,xlim=c(0,1),ylim=c(0,1))
#     points(cur_des,col="red",pch=4)

    #increment
    itcur <- itcur + 1
  }

  parallel::stopCluster(cluster)

  return(cur_des)
}

########################################################################
# Functions for generating minimax designs from images
########################################################################

# (N,p,q=10,region="uh",
#  pso=list(w=0.72,c1=1.49,c2=1.49),
#  npart=10,2*npart=5,
#  nclust=1e5,neval=10*nclust,point=NA,eval_pts=NA,
#  itmax_pso=200,itmax_pp=50,itmax_inn=1e4,
#  it_lim_pso=25,it_lim_pp=10,
#  tol_pso=1e-4,tol_pp=1e-4,tol_inn=1e-4,
#  regionby=ifelse(p>2,1e-3,-1),
#  jit=ifelse(region=="simp",0,0.1/sqrt(N)),
#  pp_flag=F)

minimax.map = function(N,img,p=2,q=10,
                      params_pso=list(w=0.72,c1=1.49,c2=1.49),
                      npart=5,nclust=1e5,neval=nclust,
                      itmax_pso=50,itmax_pp=100,itmax_inn=1e4,jit=0.1/sqrt(N)){

  #Read in image
  mapfile <- round(img)

  #Get evaluation points for minimax PSO
  eval_pts <- which(mapfile==0,arr.ind=T)
  eval_pts[,1] <- eval_pts[,1]/nrow(mapfile)
  eval_pts[,2] <- eval_pts[,2]/ncol(mapfile)

  # Get clustering data
  clust_pts <- eval_pts[sample.int(nrow(eval_pts),nclust),]
  # inv_pt <- cbind(point[,2],point[,1])
  # inv_pt[,2] <- 1 - inv_pt[,2]
  # plot(inv_pt)

  #Generate initial centers
  cluster_center <- eval_pts[sample.int(nrow(eval_pts),N),]
  for (i in 1:(npart[1]-1)){
    cluster_center = cbind(cluster_center, eval_pts[sample.int(nrow(eval_pts),N),])
  }

  # Do mMc-PSO
  bd <- c(0,1)
  D <- kmeanspso(clust_pts, eval_pts, cluster_center, q, 2.0,
                params_pso$w,params_pso$c1,params_pso$c2,
                2*npart, itmax_pso, itmax_pp, 1000, 1000, 1e-4, 1e-4,
                1e-4, itmax_inn,
                parallel::detectCores(), jit, bd[1], bd[2])
  # D$gbes_centers <- cbind(D$gbes_centers[,2],1-D$gbes_centers[,1])
  return(D$gbes_centers)
}

# if (ini=="kmeans"){ # initialize via k-means clustering
#   cluster_center = stats::kmeans(clust_pts, jitter(tf(matrix(stats::runif(N*p),ncol=p),by.clust,parallel::detectCores())) )$centers
#   for (i in 1:(npart-1)){ #add scrambled randtoolbox::sobol sets
#     cluster_center = cbind(cluster_center,
#                            stats::kmeans(clust_pts, jitter(tf(randtoolbox::sobol(N,p,init=FALSE),by.clust,parallel::detectCores())) )$centers
#     )
#   }
# }else

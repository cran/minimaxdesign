# #set compiler flags for openMP (parallel processing in C++)
# Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
# Sys.setenv("PKG_LIBS"="-fopenmp")

######################################################################
# Loading required libraries
######################################################################
library(randtoolbox)
library(DiceDesign)
library(gtools)
library(MaxPro)
library(parallel)
library(doParallel)
library(nloptr)
library(foreach)
library(Rcpp)

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

########################################################################
# Main function for mMcPSO
########################################################################
mMcPSO = function(N,p,q=10,
                  pso=list(w=0.72,c1=1.49,c2=1.49),
                  point_num=1e5,eval_num=10*point_num,point=NA,eval_pts=NA,
                  bd = c(0,1),
                  part_num=c(pso=10,pp=10),
                  it_max=c(pso=200,pp=150,inn=1e4),
                  it_lim=c(pso=25,pp=50),
                  it_tol=c(pso=1e-4,pp=1e-4,inn=sqrt(p)*1e-4),
                  region="uh",
                  regionby=ifelse(p>2,1e-3,-1),
                  jit=0.1/sqrt(N),
                  pp_flag=F){
  # Description of Inputs:
  # N                   - Number of design points desired
  # p                   - Dimension of design space
  # k                   - q in paper; approximation coefficient for minimax criterion
  # point_num           - Number of clustering data points to use
  # eval_num            - Number of evaluation points to use in minimax post-processing
  # part_num            - Number of PSO particles to use in minimax clustering
  # mM_part_num         - Number of PSO particles to use in minimax post-processing
  # it_max              - Maximum iterations for minimax clustering
  # mM_it_max           - Maximum iterations for minimax post-processing
  # inn_itmax           - Inner maximum iterations for AGD
  # trans_f             - Transformation function for randtoolbox::sobol' sequence
  # by                  - Approximation step-size for transformation
  # tol                 - Jitter tolerance
  # lb                  - Lower bound for design
  # ub                  - Upper bound for design
  # eval_pts            - Pre-computed evaluation points (NA - compute within function)
  # point               - Pre-computed clustering data (NA - compute within function)

  #Setting transformation type
  # print(1)
  switch(region,
         uh = {
           tf <- function(D,by,num_proc){return(D)}
           checkBounds <- function(D){ #function for checking bound
             return(D)
           }
         },
         simp = {
           tf <- CtoA;
           checkBounds <- function(D){ #function for checking bound
            return(D)
           }
         },
         ball = {
           # if (p==2){
           #   tf <- CtoB2;
           # }
           # else{
           #   tf <- CtoBp;
           # }
           tf <- CtoB;
           bd <- c(-1,1);
           checkBounds <- function(D){ #function for checking bound
              good.idx <- which(rowSums(D^2)<=1)
              return(D[good.idx,])
           }
         })
  by <- regionby

  # print(2.1)
  # print(by)
  # print(parallel::detectCores())
  print("Generating approximating point sets ...")
  if (is.na(point)){
    #Generate point
    point = tf(randtoolbox::sobol(point_num,p),by,parallel::detectCores()) #clustering points
    point <- checkBounds(point)
  }

  # print(2.2)
  if (is.na(eval_pts)){
    #Generate eval_pts
    eval_pts = tf(randtoolbox::sobol(eval_num,p,init=FALSE),by,parallel::detectCores()) #minimax approximation points for post-processing
    eval_pts <- checkBounds(eval_pts)
  }

  #   offset = seq(0,1-1/part_num,1/part_num)

  # print(3)
  # Generate initial center particles
  by.clust <- 1e-4
  if (pp_flag){
    cluster_center = stats::kmeans(point, jitter(tf(matrix(stats::runif(N*p),ncol=p),by.clust,parallel::detectCores())) )$centers
    for (i in 1:(part_num[1]-1)){ #add scrambled randtoolbox::sobol sets
      cluster_center = cbind(cluster_center,
                             stats::kmeans(point, jitter(tf(matrix(stats::runif(N*p),ncol=p),by.clust,parallel::detectCores())) )$centers
                             )
    }
  }else{
    cluster_center = tf(randtoolbox::sobol(N,p),by.clust,parallel::detectCores()) #initial cluster centers
    for (i in 1:(part_num[1]-1)){ #add scrambled randtoolbox::sobol sets
      cluster_center = cbind(cluster_center,tf(randtoolbox::sobol(N,p,init=FALSE),by.clust,parallel::detectCores()))
    }
  }

  # print(4)
  # mMc-PSO: function coded in C++
  t1 = Sys.time()
  D = kmeanspso(point, eval_pts, cluster_center, q, 2.0,
                pso$w,pso$c1,pso$c2,
                part_num[2], it_max[1], it_max[2], it_lim[1], it_lim[2], it_tol[1], it_tol[2],
                it_tol[3], it_max[3],
                parallel::detectCores(), jit, bd[1], bd[2])
  fitpre <- D$gbes_centers
  fitpost <- fitpre;

  t2 = Sys.time()
  tm <- difftime(t2,t1)
  # print(5)

  # D$gbes_centers <- fitpost
  # D$time <- tm
  return(D$gbes_centers)
}

########################################################################
# Functions for minimax projection designs
########################################################################

#Main function for generating minimax projection designs
miniMaxPro <- function(N,p,mMdes=NA,it_max=c(pso=100,pp=100,inn=1e4),
                       refine_num=1e5, refine_pts=NA, refine_itmax=100, mM_tol=1e-3*p, ...){
  #Generate minimax design from mMc-PSO
  if (is.na(mMdes)){
    mMdes <- mMcPSO(N,p,it_max=it_max, ...)
  }
  if (is.na(refine_pts)){
    refine_pts <- randtoolbox::sobol(refine_num,p)
  }

  #Refinement
  ret <- refine(mMdes,eval_pts=refine_pts,pp_itmax=refine_itmax,mM_tol=mM_tol)

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

  while (itcur < pp_itmax){
    print(paste0("MaxPro adjustment iteration: ", itcur))

    #Blockwise MaxPro for each design point
    dst_vec <- mMcrit_allpts( cur_des, eval_pts ) #First compute the minimum distance for each design point
    dst_vec <- pmax(orig_mM + mM_tol - dst_vec,0)
    mM_ind <- which(dst_vec==0)
    prev_des <- cur_des

    indices <- gtools::permute(setdiff(1:N,mM_ind))
    flg_cont <- F

    while (!flg_cont){
      biglist <- foreach::'%dopar%'(foreach::foreach (i = indices, .packages="nloptr"),
      {
        toler <- dst_vec[i]/jitter.scl # Ensures the jittered initial points do not exceed minimax
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

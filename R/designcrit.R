# ########################################################################
# # Computation of minimax criterion
# ########################################################################
# mMcrit <- function(D,eval_pts){
#   return( mMcritPt(D,eval_pts) );
# }
#
# Mmcrit_idx <- function(D){
#   dim_num <- ncol(D)
#   point_num <- nrow(D)
#   ret <- ( 1.0/choose(point_num,2)*sum(1/(stats::dist(D)^(2*dim_num))) ) ^ (-1/(2*dim_num))
#   return(ret)
# }
#
# mMcrit_pr <- function(D,k,eval_pts){
#   #Projected minimax criterion
#   p <- ncol(D)
#   N <- nrow(D)
#   indices <- t(utils::combn(1:p, k))
#   # if (k == 1){
#   #   eval_pts <- matrix(randtoolbox::sobol(eval_num,k),ncol=1)
#   # }else{
#   #   eval_pts <- randtoolbox::sobol(eval_num,k)
#   # }
#   # eval_pts <- randtoolbox::sobol(eval_num,p)
#
#   return(mMcrit_proj(D,eval_pts,indices))
# }
#
# avgcrit_pr <- function(D,k,eval_pts){
#   #Projected minimax criterion
#   p <- ncol(D)
#   N <- nrow(D)
#   indices <- t(utils::combn(1:p, k))
#   # if (k == 1){
#   #   eval_pts <- matrix(randtoolbox::sobol(eval_num,k),ncol=1)
#   # }else{
#   #   eval_pts <- randtoolbox::sobol(eval_num,k)
#   # }
#   # eval_pts <- randtoolbox::sobol(eval_num,p)
#
#   return(avgcrit_proj(D,eval_pts,indices))
# }
#
# Mmcrit_pr <- function(D,k){
#   #Projected minimax criterion
#   p <- ncol(D)
#   N <- nrow(D)
#   indices <- t(utils::combn(1:p, k))
#   dst <- 1e10 #arbitrarily large
#
#   for (i in 1:nrow(indices)){
#     dst <- min(dst,Mmcrit_idx( matrix(D[,indices[i,]],ncol=k) ) )
#   }
#
#   return(dst);
# }

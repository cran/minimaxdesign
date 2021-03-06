# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

CtoBp <- function(D, by, num_proc) {
    .Call('_minimaxdesign_CtoBp', PACKAGE = 'minimaxdesign', D, by, num_proc)
}

printBar <- function(prop) {
    invisible(.Call('_minimaxdesign_printBar', PACKAGE = 'minimaxdesign', prop))
}

kmeansreg <- function(Rcpp_point, Rcpp_cluster_center, p, pw, it_max, inn_tol, num_proc, fix_ind) {
    .Call('_minimaxdesign_kmeansreg', PACKAGE = 'minimaxdesign', Rcpp_point, Rcpp_cluster_center, p, pw, it_max, inn_tol, num_proc, fix_ind)
}

kmeanspso <- function(Rcpp_point, Rcpp_evalpts, Rcpp_cluster_center, p, pw, w, c1, c2, mM_part_num, it_max, mM_it_max, it_lim, mM_it_lim, it_tol, mM_it_tol, inn_tol, inn_itmax, num_proc, tol, lb, ub, Rcpp_fix_ind) {
    .Call('_minimaxdesign_kmeanspso', PACKAGE = 'minimaxdesign', Rcpp_point, Rcpp_evalpts, Rcpp_cluster_center, p, pw, w, c1, c2, mM_part_num, it_max, mM_it_max, it_lim, mM_it_lim, it_tol, mM_it_tol, inn_tol, inn_itmax, num_proc, tol, lb, ub, Rcpp_fix_ind)
}

mMcritPt <- function(Rcpp_point, Rcpp_evalpts) {
    .Call('_minimaxdesign_mMcritPt', PACKAGE = 'minimaxdesign', Rcpp_point, Rcpp_evalpts)
}

mMcrit_allpts <- function(Rcpp_point, Rcpp_evalpts) {
    .Call('_minimaxdesign_mMcrit_allpts', PACKAGE = 'minimaxdesign', Rcpp_point, Rcpp_evalpts)
}

mMcrit <- function(point, Rcpp_evalpts) {
    .Call('_minimaxdesign_mMcrit', PACKAGE = 'minimaxdesign', point, Rcpp_evalpts)
}

mMcrit_proj <- function(Rcpp_pts, Rcpp_evalpts, indices) {
    .Call('_minimaxdesign_mMcrit_proj', PACKAGE = 'minimaxdesign', Rcpp_pts, Rcpp_evalpts, indices)
}

avgcrit_proj <- function(Rcpp_pts, Rcpp_evalpts, indices) {
    .Call('_minimaxdesign_avgcrit_proj', PACKAGE = 'minimaxdesign', Rcpp_pts, Rcpp_evalpts, indices)
}

kmeansobj <- function(Rcpp_point, Rcpp_evalpts, p) {
    .Call('_minimaxdesign_kmeansobj', PACKAGE = 'minimaxdesign', Rcpp_point, Rcpp_evalpts, p)
}

CtoAA <- function(D, by, num_proc) {
    .Call('_minimaxdesign_CtoAA', PACKAGE = 'minimaxdesign', D, by, num_proc)
}

CtoB2 <- function(D, by, num_proc) {
    .Call('_minimaxdesign_CtoB2', PACKAGE = 'minimaxdesign', D, by, num_proc)
}

closestPt <- function(points, grd) {
    .Call('_minimaxdesign_closestPt', PACKAGE = 'minimaxdesign', points, grd)
}


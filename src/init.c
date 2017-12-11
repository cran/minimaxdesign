#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _minimaxdesign_avgcrit_proj(SEXP, SEXP, SEXP);
extern SEXP _minimaxdesign_closestPt(SEXP, SEXP);
extern SEXP _minimaxdesign_CtoAA(SEXP, SEXP, SEXP);
extern SEXP _minimaxdesign_CtoB2(SEXP, SEXP, SEXP);
extern SEXP _minimaxdesign_CtoBp(SEXP, SEXP, SEXP);
extern SEXP _minimaxdesign_kmeansobj(SEXP, SEXP, SEXP);
extern SEXP _minimaxdesign_kmeanspso(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _minimaxdesign_kmeansreg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _minimaxdesign_mMcrit_allpts(SEXP, SEXP);
extern SEXP _minimaxdesign_mMcrit_proj(SEXP, SEXP, SEXP);
extern SEXP _minimaxdesign_mMcritPt(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_minimaxdesign_avgcrit_proj",  (DL_FUNC) &_minimaxdesign_avgcrit_proj,   3},
    {"_minimaxdesign_closestPt",     (DL_FUNC) &_minimaxdesign_closestPt,      2},
    {"_minimaxdesign_CtoAA",         (DL_FUNC) &_minimaxdesign_CtoAA,          3},
    {"_minimaxdesign_CtoB2",         (DL_FUNC) &_minimaxdesign_CtoB2,          3},
    {"_minimaxdesign_CtoBp",         (DL_FUNC) &_minimaxdesign_CtoBp,          3},
    {"_minimaxdesign_kmeansobj",     (DL_FUNC) &_minimaxdesign_kmeansobj,      3},
    {"_minimaxdesign_kmeanspso",     (DL_FUNC) &_minimaxdesign_kmeanspso,     21},
    {"_minimaxdesign_kmeansreg",     (DL_FUNC) &_minimaxdesign_kmeansreg,      7},
    {"_minimaxdesign_mMcrit_allpts", (DL_FUNC) &_minimaxdesign_mMcrit_allpts,  2},
    {"_minimaxdesign_mMcrit_proj",   (DL_FUNC) &_minimaxdesign_mMcrit_proj,    3},
    {"_minimaxdesign_mMcritPt",      (DL_FUNC) &_minimaxdesign_mMcritPt,       2},
    {NULL, NULL, 0}
};

void R_init_minimaxdesign(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

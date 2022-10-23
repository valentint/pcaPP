//  VT::27.06.2017 - this file was added to fix the warning "Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols'"
//
//  About registration of native symbols see for example: https://www.r-bloggers.com/1-easy-package-registration/
//      also here http://r.789695.n4.nabble.com/Registration-of-native-routines-td4728874.html
//      - about Windows - take the 64 bit version of mingw!
//

#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include <R.h>              
#include <Rinternals.h>

/*
	EXPORT void C_kendallNlogN (double* arr1, double* arr2, int *pnPar, double *dRet) ;
	EXPORT void C_l1median_BFGS (int *pnParam_In, int *pnParam_Out, double *pdParam_In, double *pdParam_Out, double *pdData, double *pdMRet) ;
	EXPORT void C_l1median_CG(int *pnParam, int *pnParam_Out, double *pdParam, double *pdParam_Out, double *pdData, double *pdMRet) ;
	EXPORT void C_l1median_HoCr	(int *pnParIn, int *pnParOut, double *pdParIn, double *pdX, double *pdMed) ;
	EXPORT void C_l1median_NLM (int *pnParam, double *pdParam, double *pdData, double *pdMRet) ;
	EXPORT void C_l1median_NLM_Hess (int *pnParam, double *pdParam, double *pdData, double *pdMRet) ;
	EXPORT void C_l1median_NM(int *pnParam, double *pdParam, double *pdData, double *pdMRet) ;
	EXPORT void C_l1Median_VZ		(int *pnParIn, int *pnParOut, double *pdParIn, double *pdX, double *pdMed) ;
	EXPORT void C_PCAgrid (int *pnParamIn, int *pnParamOut, double *pdParamIn, double *pdData, double *pdLoadings, double *pdSDev, double *pdObj) ;
	EXPORT void C_sPCAgrid (int *pnParamIn, int *pnParamOut, double *pdParamIn, double *pdData, double *pdLoadings, double *pdSDev, double *pdObj, double *pdLambda, double *pdBackTransHD) ;
	EXPORT void C_pcaProj (int *pnParIn, double *pdParIn, double *pdX, double *pdZ, double *pdL, double *pdSDev) ;
	EXPORT void C_pcaProj_up (int *pnParIn, double *pdParIn, double *pdX, double *pdZ, double *pdL, double *pdSDev) ;
	EXPORT void C_qn (int *pnParIn, double *pdParIn, double *pdParOut, double *pdX) ;

*/
 
/* .C calls */
extern void C_kendallNlogN(void *, void *, void *, void *);
extern void C_l1median_BFGS(void *, void *, void *, void *, void *, void *);
extern void C_l1median_CG(void *, void *, void *, void *, void *, void *);
extern void C_l1median_HoCr(void *, void *, void *, void *, void *);
extern void C_l1median_NLM(void *, void *, void *, void *);
extern void C_l1median_NLM_Hess(void *, void *, void *, void *);
extern void C_l1median_NM(void *, void *, void *, void *);
extern void C_l1Median_VZ(void *, void *, void *, void *, void *);
extern void C_PCAgrid(void *, void *, void *, void *, void *, void *, void *);
extern void C_pcaProj(void *, void *, void *, void *, void *, void *);
extern void C_pcaProj_up(void *, void *, void *, void *, void *, void *);
extern void C_qn(void *, void *, void *, void *);
extern void C_sPCAgrid(void *, void *, void *, void *, void *, void *, void *, void *, void *);

static R_NativePrimitiveArgType C_kendallNlogN_t[] = {
    REALSXP, REALSXP, INTSXP, REALSXP
};

static R_NativePrimitiveArgType C_l1median_BFGS_t[] = {
    INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType C_l1median_CG_t[] = {
    INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType C_l1median_HoCr_t[] = {
    INTSXP, INTSXP, REALSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType C_l1median_NLM_t[] = {
    INTSXP, REALSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType C_l1median_NLM_Hess_t[] = {
    INTSXP, REALSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType C_l1median_NM_t[] = {
    INTSXP, REALSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType C_l1median_VZ_t[] = {
    INTSXP, INTSXP, REALSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType C_PCAgrid_t[] = {
    INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType C_sPCAgrid_t[] = {
    INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType C_pcaProj_t[] = {
    INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType C_pcaProj_up_t[] = {
    INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType C_qn_t[] = {
    INTSXP, REALSXP, REALSXP, REALSXP
};

static const R_CMethodDef CEntries[] = {
    {"C_kendallNlogN",      (DL_FUNC) &C_kendallNlogN,      4, C_kendallNlogN_t},
    {"C_l1median_BFGS",     (DL_FUNC) &C_l1median_BFGS,     6, C_l1median_BFGS_t},
    {"C_l1median_CG",       (DL_FUNC) &C_l1median_CG,       6, C_l1median_CG_t},
    {"C_l1median_HoCr",     (DL_FUNC) &C_l1median_HoCr,     5, C_l1median_HoCr_t},
    {"C_l1median_NLM",      (DL_FUNC) &C_l1median_NLM,      4, C_l1median_NLM_t},
    {"C_l1median_NLM_Hess", (DL_FUNC) &C_l1median_NLM_Hess, 4, C_l1median_NLM_Hess_t},
    {"C_l1median_NM",       (DL_FUNC) &C_l1median_NM,       4, C_l1median_NM_t},
    {"C_l1Median_VZ",       (DL_FUNC) &C_l1Median_VZ,       5, C_l1median_VZ_t},
    {"C_PCAgrid",           (DL_FUNC) &C_PCAgrid,           7, C_PCAgrid_t},
    {"C_pcaProj",           (DL_FUNC) &C_pcaProj,           6, C_pcaProj_t},
    {"C_pcaProj_up",        (DL_FUNC) &C_pcaProj_up,        6, C_pcaProj_up_t},
    {"C_qn",                (DL_FUNC) &C_qn,                4, C_qn_t},
    {"C_sPCAgrid",          (DL_FUNC) &C_sPCAgrid,          9, C_sPCAgrid_t},
    {NULL, NULL, 0}
};

void R_init_pcaPP(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}


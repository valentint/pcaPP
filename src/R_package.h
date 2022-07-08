#ifdef R_PACKAGE_FILE

#ifdef ES_DEV_ENV
	#include "../../../RDev/R.Inc.h"
	#include "../../../SMat/smat.def.h"
#else
	#include "R.Inc.h"
	#include "smat.def.h"
#endif



/////////////
//	pcaPP  //
/////////////

	EXPORT void C_pcaProj_up (int *pnParIn, double *pdParIn, double *pdX, double *pdZ, double *pdL, double *pdSDev) ;
	EXPORT void C_pcaProj (int *pnParIn, double *pdParIn, double *pdX, double *pdZ, double *pdL, double *pdSDev) ;

	EXPORT void C_sPCAgrid (int *pnParamIn, int *pnParamOut, double *pdParamIn, double *pdData, double *pdLoadings, double *pdSDev, double *pdObj/*, double *pdMaxMaha*/, double *pdLambda, double *pdBackTransHD) ;
	EXPORT void C_PCAgrid (int *pnParamIn, int *pnParamOut, double *pdParamIn, double *pdData, double *pdLoadings, double *pdSDev, double *pdObj/*, double *pdMaxMaha*/) ;

/////////////////
//	L1 Median  //
/////////////////

	EXPORT void Hess_Sub_R (int *pnPar, double *pdX_i, double *pdMu, double *pdHess) ;
	EXPORT void Hess_R (int *pnPar, double *pdX, double *pdMu, double *pdHess) ;

	EXPORT void C_l1Median_VZ		(int *pnParIn, int *pnParOut, double *pdParIn, double *pdX, double *pdMed/*, double *pdWeights*/) ;
	EXPORT void C_l1median_HoCr	(int *pnParIn, int *pnParOut, double *pdParIn, double *pdX, double *pdMed) ;

	EXPORT void C_l1median_NM(int *pnParam, double *pdParam, double *pdData, /*double *pdParScale, */double *pdMRet) ;
	EXPORT void C_l1median_CG(int *pnParam, int *pnParam_Out, double *pdParam, double *pdParam_Out, double *pdData, /*double *pdParScale, */double *pdMRet) ;
	EXPORT void C_l1median_BFGS (int *pnParam_In, int *pnParam_Out, double *pdParam_In, double *pdParam_Out, double *pdData, /*double *pdParScale, */double *pdMRet) ;
	EXPORT void l1median_SA (int *pnParam_In, int *pnParam_Out, double *pdParam_In, double *pdParam_Out, double *pdData, /*double *pdParScale, */double *pdMRet) ;

	EXPORT void C_l1median_NLM (int *pnParam, double *pdParam, double *pdData, double *pdMRet/*, double *pdTypSize*/) ;
	EXPORT void C_l1median_NLM_Hess (int *pnParam, double *pdParam, double *pdData, double *pdMRet/*, double *pdTypSize*/) ;

//////////
//	Qn  //
//////////

	EXPORT void C_qn (int *pnParIn, double *pdParIn, double *pdParOut, double *pdX) ;

///////////////////
//	Cov Kendall  //
///////////////////

	EXPORT void C_kendallNlogN (double* arr1, double* arr2, int *pnPar, double *dRet) ;

//////////////
//	SDoOut  //
//////////////

	EXPORT void SDoOut (int *pnParIn, double *pdX, double *pdMaxMaha, int *pnNChanged) ;


#endif	//	#ifdef R_PACKAGE_FILE

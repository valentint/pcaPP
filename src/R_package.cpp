#ifdef R_PACKAGE_FILE

#include "R_package.h"
#include "cov.kendall.h"

#include "PCAgrid.h"
#include "PCAproj.h"
#include "outSDo.h"

#include "L1Median.h"
#include "qnn.h"

#ifdef ES_DEV_ENV
	#include "..\..\..\RDev\perftimer.h"
	#include "..\..\..\RDev\R_meal.h"
#else
	#include "perftimer.h"
	#include "R_meal.h"
#endif	//	#ifdef ES_DEV_ENV

//R_MEAL_SETTINGS ("P.Filzmoser@tuwien.ac.at") ;	//	settings for the R meal - implementation
R_MEAL_SETTINGS ("Heinrich_Fritz@hotmail.com") ;	//	settings for the R meal - implementation

////////////////////////////////
//	exporting functions to R  //
////////////////////////////////

#ifndef ES_DEV_ENV
	void C_kendallNlogN (double* arr1, double* arr2, int *pnPar, double *dRet)		//	2do: wrap other fct (N^2) too; make choice depending on n (pnPar[0])
	{
		TRY( 
			*dRet = kendallNlogN (arr1, arr2,		//	arr1, arr2
								(size_t) pnPar[0],	//	length
								pnPar[1]) ;			//	cor
			)
	}
#endif

	void C_PCAgrid (int *pnParamIn, int *pnParamOut, double *pdParamIn, double *pdData, double *pdLoadings, double *pdSDev, double *pdObj/*, double *pdMaxMaha*/)
	{
		TRY( 
			CPCAGrid (pnParamIn, pnParamOut, pdParamIn, pdData, pdLoadings, pdSDev, pdObj/*, pdMaxMaha*/).Calc () ;
			)
	}

	void C_sPCAgrid (int *pnParamIn, int *pnParamOut, double *pdParamIn, double *pdData, double *pdLoadings, double *pdSDev, double *pdObj/*, double *pdMaxMaha*/, double *pdLambda, double *pdBackTransHD)
	{
		TRY( 
			CsPCAGrid (pnParamIn, pnParamOut, pdParamIn, pdData, pdLoadings, pdSDev, pdObj/*, pdMaxMaha*/, pdLambda, pdBackTransHD).Calc () ;
			)
	}

	void C_pcaProj_up (int *pnParIn, double *pdParIn, double *pdX, double *pdZ, double *pdL, double *pdSDev)
	{
		TRY( 
			CPCAprojU (pnParIn, pdParIn, pdX, pdZ, pdL, pdSDev).Calc () ;
			)
	}

	void C_pcaProj (int *pnParIn, double *pdParIn, double *pdX, double *pdZ, double *pdL, double *pdSDev)
	{
		TRY( 
			CPCAproj (pnParIn, pdParIn, pdX, pdZ, pdL, pdSDev).Calc () ;
			)
	}

	void C_l1Median_VZ (int *pnParIn, int *pnParOut, double *pdParIn, double *pdX, double *pdMed/*, double *pdWeights*/)
	{
		TRY( 
			CPerfTimer tim ;
			CL1Median_VZ (pnParIn, pnParOut, pdParIn, pdX, pdMed, NULL) ; //, pdWeights) ;
			pnParOut[2] = tim.GetTimeMS () ;
			)
	}

	void C_l1median_HoCr (int *pnParIn, int *pnParOut, double *pdParIn, double *pdX, double *pdMed)
	{
		const int &n = pnParIn [0], &p = pnParIn [1], &dwMaxit = pnParIn[2], &dwTrace = pnParIn [3] ;
		int &nCode = pnParOut[0] = 0, &dwIterCount = pnParOut [1] ;

		const double &dTol = pdParIn[0], &dZeroTol = pdParIn[1] ;

		TRY( 
			CPerfTimer tim ;
			nCode = l1median_HoCr (SMatD (pdX, n, p), SVecD (pdMed, p), dZeroTol, dTol, dwMaxit, dwTrace, &dwIterCount) ;
			pnParOut[2] = tim.GetTimeMS () ;
			)

		return ;
	}

	void C_qn (int *pnParIn, double *pdParIn, double *pdParOut, double *pdX)
	{
		int &n = pnParIn[0] ;
		double &dCorrFact = pdParIn [0] ;
		double &dQn = pdParOut [0] ;

		TRY( 
			dQn = qn_raw (pdX, n) ;
			dQn *= qn_corrN (n, dCorrFact) ;
			)
	}

	void SDoOut (int *pnParIn, double *pdX, double *pdMaxMaha, int *pnNChanged)
	{
		TRY( 
			CSDoOut (pnParIn, pdX, pdMaxMaha, pnNChanged).Calc () ;
			)
	}

#endif	//	#ifdef R_PACKAGE_FILE

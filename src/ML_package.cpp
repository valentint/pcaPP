#ifdef MATLAB_MEX_FILE


#include "ML_package.h"
#include "L1Median.h"

#include "qnn.h"

#include "PCAgrid.h"
#include "PCAproj.h"


	void mexFunction(int nlhs, mxArray* plhs[], int nrhs, mxArray *prhs[])
	{
		if (nrhs < 1)
		{
			meal_error ("At least one argument needed.") ;
			return ;
		}

		mxArray *pSwitch = prhs [0] ;

		if (mxGetM (pSwitch) != 1 || 
			mxGetN (pSwitch) != 1)
		{
			meal_error ("First argument has to be of size 1x1.") ;
			return ;
		}

		int nSwitch = (int) *mxGetPr (pSwitch) ;

		switch (nSwitch)
		{
		case 0:
			ML_l1median_HoCr (nlhs, plhs, nrhs - 1, prhs + 1) ;
			break ;

		case 1:
			ML_l1median_VaZh (nlhs, plhs, nrhs - 1, prhs + 1) ;
			break ;

		case 2:
			ML_qn (nlhs, plhs, nrhs - 1, prhs + 1) ;
			break ;

		case 3:
			ML_PCAgrid (nlhs, plhs, nrhs - 1, prhs + 1) ;
			break ;

		case 4:
			ML_sPCAgrid (nlhs, plhs, nrhs - 1, prhs + 1) ;
			break ;

		case 5:
			ML_PCAproj (nlhs, plhs, nrhs - 1, prhs + 1) ;
			break ;

		case 6:
			ML_PCAprojU (nlhs, plhs, nrhs - 1, prhs + 1) ;
			break ;

		
		default:
			meal_error ("Unknown switch argument received.") ;

			break ;
				

		} ;
	}

	void ML_l1median_HoCr (int nlhs, mxArray* plhs[], int nrhs, mxArray *prhs[])
		//(int *pnParIn, int *pnParOut, double *pdParIn, double *pdX, double *pdMed)
	{
		if (nlhs != 1)
			meal_error ("1 output argument expected") ;

		if (nrhs != 3)
			meal_error ("3 input arguments expected") ;

			//	Data Matrix

		mxArray *mxX = prhs [0] ;
		const int n = mxGetM (mxX), p = mxGetN (mxX) ;
		double *pdX = mxGetPr (mxX) ;

			//	Median Vector

		mxArray *mxMed = prhs [1] ;

		if (mxGetM (mxMed) * mxGetN (mxMed) != (t_size) p)
			meal_error ("Length of the median vector is expected to be p.") ;
		double *pdMed = mxGetPr (mxMed) ;

			//	Parameter Vector
		mxArray *mxPar = prhs [2] ;

		if (mxGetM (mxPar) * mxGetN (mxPar)!= 4)
			meal_error ("Parameter array of dimension 5x1 expected") ;

		double *pdParIn = mxGetPr (mxPar) ;
		const int &dwMaxit = (int) pdParIn[0], dwTrace = (int) pdParIn[1] ;
		const double &dTol = pdParIn[2], &dZeroTol = pdParIn[3] ;

		plhs [0] = mxCreateDoubleMatrix (1, 2, mxREAL) ;
		double *pdParOut = mxGetPr (plhs [0]) ;
		double &dCode = pdParOut [0] ;
		

		int nIterCount ;

		TRY( 
			dCode = l1median_HoCr (SMatD (pdX, n, p), SVecD (pdMed, p), dZeroTol, dTol, dwMaxit, dwTrace, &nIterCount) ;
			)

		//	Out Parameter Vector
		pdParOut [1] = nIterCount ;

		return ;
	}


	void ML_l1median_VaZh (int nlhs, mxArray* plhs[], int nrhs, mxArray *prhs[])
		//(int *pnParIn, int *pnParOut, double *pdParIn, double *pdX, double *pdMed)
	{
		if (nlhs != 1)
			meal_error ("1 output argument expected") ;

		if (nrhs != 3)
			meal_error ("3 input arguments expected") ;

			//	Data Matrix

		mxArray *mxX = prhs [0] ;
		const int n = mxGetM (mxX), p = mxGetN (mxX) ;
		double *pdX = mxGetPr (mxX) ;

			//	Median Vector

		mxArray *mxMed = prhs [1] ;

		if (mxGetM (mxMed) * mxGetN (mxMed) != (t_size) p)
			meal_error ("Length of the median vector is expected to be p.") ;
		double *pdMed = mxGetPr (mxMed) ;

			//	Parameter Vector
		mxArray *mxPar = prhs [2] ;

		if (mxGetM (mxPar) * mxGetN (mxPar)!= 4)
			meal_error ("Parameter array of dimension 5x1 expected") ;

		double *pdParIn = mxGetPr (mxPar) ;
		const int &dwMaxit = (int) pdParIn[0], dwTrace = (int) pdParIn[1] ;
		const double &dTol = pdParIn[2], &dZeroTol = pdParIn[3] ;

		int nCode, nIterCount ;

		TRY( 
			CL1Median_VZ (n, p, nCode, nIterCount, pdParIn, pdX, pdMed, NULL) ; //, pdWeights) ;
			)

		//	Out Parameter Vector
		plhs [0] = mxCreateDoubleMatrix (1, 2, mxREAL) ;
		double *pdParOut = mxGetPr (plhs [0]) ;
		pdParOut [0] = nCode ;
		pdParOut [1] = nIterCount ;

		return ;
	}


	void ML_qn (int nlhs, mxArray* plhs[], int nrhs, mxArray *prhs[])
	{

		if (nlhs != 1)
			meal_error ("1 output argument expected") ;

		if (nrhs != 2)
			meal_error ("2 input arguments expected") ;

	//	Get data

		const int n = mxGetM (prhs[0]) * mxGetN (prhs[0]) ;
		double *pdX = mxGetPr (prhs[0]) ;

	//	Get Corrfact

		if (mxGetM (prhs[1]) != 1 ||
			mxGetN (prhs[1]) != 1)
			meal_error  ("Parameter 2 is supposed to be a scalar.") ;

		const double dCorrFact = *mxGetPr (prhs[1]) ;

	//	GetQn (output)

		plhs[0] = mxCreateDoubleMatrix (1, 1, mxREAL) ;
		double &dQn = *mxGetPr (plhs [0]) ;

		TRY( 
			qn_nc (dQn, pdX, n) ;
			dQn *= qn_corrN (n, dCorrFact) ;
			)
	}

	void ML_PCAgrid (int nlhs, mxArray* plhs[], int nrhs, mxArray *prhs[])
	{

		if (nlhs != 3)
			meal_error ("3 output argument expected") ;

		if (nrhs != 4)
			meal_error ("4 input arguments expected") ;
		

	//	Get pnParIn
		if (mxGetM (prhs[0]) * mxGetN (prhs[0]) != 7)
			meal_error ("pnParamIn: vector of length 7 expected.") ;

		double *pdnParamIn = mxGetPr (prhs[0]) ;
		int pnParamIn[9] ;
		Double2Int (pnParamIn + 2, pdnParamIn, 7) ;

		int nk = pnParamIn[2] ;

	//	Get pdParIn

		if (mxGetM (prhs[1]) * mxGetN (prhs[1]) != 1)
			meal_error ("pdParIn: vector of length 1 expected.") ;

		double *pdParamIn = mxGetPr (prhs[1]) ;

	//	Get pdData

		int &n = pnParamIn[0], &p = pnParamIn[1] ;

		n = mxGetM (prhs[2]) ;
		p = mxGetN (prhs[2]) ;

		double *pdData = mxGetPr (prhs[2]) ;

	//	Get pdLoadings

		if (mxGetM (prhs[3]) != (t_size) p ||
			mxGetN (prhs[3]) != (t_size)	p)
			meal_error ("pdLoadings: matrix of dimension pxp expected.") ;

		double *pdLoadings = mxGetPr (prhs[3]) ;


	//	Get pnParOut (output)
		plhs[0] = mxCreateDoubleMatrix (1, 1, mxREAL) ;
		double *pdnParamOut = mxGetPr (plhs[0]) ;
		int pnParamOut [1] ;

	//	Get pdSdev (output)

		plhs[1] = mxCreateDoubleMatrix (1, nk, mxREAL) ;
		double *pdSDev = mxGetPr (plhs[1]) ;

	//	Get	pdObj (output)

		plhs[2] = mxCreateDoubleMatrix (1, nk, mxREAL) ;
		double *pdObj = mxGetPr (plhs[2]) ;


		TRY( 
			CPCAGrid (pnParamIn, pnParamOut, pdParamIn, pdData, pdLoadings, pdSDev, pdObj/*, pdMaxMaha*/).Calc () ;
			)		

		Int2Double (pdnParamOut, pnParamOut, 1) ;

	}


	void ML_sPCAgrid (int nlhs, mxArray* plhs[], int nrhs, mxArray *prhs[])
	{

		if (nlhs != 3)
			meal_error ("expected 3 output argument") ;

		if (nrhs != 6)
			meal_error ("expected 6 input arguments") ;
		

	//	Get pnParIn
		if (mxGetM (prhs[0]) * mxGetN (prhs[0]) != 10)
			meal_error ("pnParamIn: expected vector of length 10.") ;

		double *pdnParamIn = mxGetPr (prhs[0]) ;
		int pnParamIn[12] ;
		Double2Int (pnParamIn + 2, pdnParamIn, 10) ;

		int nk = pnParamIn[2] ;
		int nPhd = pnParamIn[10] ;

	//	Get pdParIn

		if (mxGetM (prhs[1]) * mxGetN (prhs[1]) != 1)
			meal_error ("pdParIn: expected vector of length 1.") ;

		double *pdParamIn = mxGetPr (prhs[1]) ;

	//	Get pdData

		int &n = pnParamIn[0], &p = pnParamIn[1] ;

		n = mxGetM (prhs[2]) ;
		p = mxGetN (prhs[2]) ;

		double *pdData = mxGetPr (prhs[2]) ;

	//	Get pdLoadings

		if (mxGetM (prhs[3]) != (t_size) p ||
			mxGetN (prhs[3]) != (t_size)	p)
			meal_error ("pdLoadings: expected matrix of dimension pxp.") ;

		double *pdLoadings = mxGetPr (prhs[3]) ;

	//	Get pdLambda
		if (mxGetM (prhs[4]) * mxGetN (prhs[4]) != (t_size)	nk)
			meal_error ("pdLambda: expected vector of length k.") ;

		double *pdLambda = mxGetPr (prhs[4]) ;

	//	Get pdBackTransHD

		double *pdBackTransHD = NULL ;
		if( mxGetM (prhs[5]) * mxGetN (prhs[5]))		//	only if this matrix holds any values..
		{
			if (mxGetM (prhs[5]) != (t_size) nPhd ||
				mxGetN (prhs[5]) != (t_size) p)

				meal_error ("pdBackTransHD: expected matrix of dimension k.") ;

			pdBackTransHD = mxGetPr (prhs[5]) ;
		}


	//	Get pnParOut (output)
		plhs[0] = mxCreateDoubleMatrix (1, 1, mxREAL) ;
		double *pdnParamOut = mxGetPr (plhs[0]) ;
		int pnParamOut [1] ;

	//	Get pdSdev (output)

		plhs[1] = mxCreateDoubleMatrix (1, nk, mxREAL) ;
		double *pdSDev = mxGetPr (plhs[1]) ;

	//	Get	pdObj (output)

		plhs[2] = mxCreateDoubleMatrix (1, nk, mxREAL) ;
		double *pdObj = mxGetPr (plhs[2]) ;


		TRY( 
			CsPCAGrid (pnParamIn, pnParamOut, pdParamIn, pdData, pdLoadings, pdSDev, pdObj/*, pdMaxMaha*/, pdLambda, pdBackTransHD).Calc () ;
			)		

		Int2Double (pdnParamOut, pnParamOut, 1) ;
	}

	void ML_PCAproj (int nlhs, mxArray* plhs[], int nrhs, mxArray *prhs[])
	{

		if (nlhs != 3)
			meal_error ("3 output argument expected") ;

		if (nrhs != 3)
			meal_error ("3 input arguments expected") ;
		

	//	Get pnParIn
		if (mxGetM (prhs[0]) * mxGetN (prhs[0]) != 4)
			meal_error ("pnParamIn: vector of length 4 expected.") ;

		double *pdnParamIn = mxGetPr (prhs[0]) ;
		int pnParamIn[6] ;
		Double2Int (pnParamIn + 2, pdnParamIn, 4) ;

		int nk = pnParamIn[3] ;

	//	Get pdParIn

		if (mxGetM (prhs[1]) * mxGetN (prhs[1]) != 1)
			meal_error ("pdParIn: vector of length 1 expected.") ;

		double *pdParamIn = mxGetPr (prhs[1]) ;

	//	Get pdX

		int &n = pnParamIn[0], &p = pnParamIn[1] ;

		n = mxGetM (prhs[2]) ;
		p = mxGetN (prhs[2]) ;

		double *pdX = mxGetPr (prhs[2]) ;

	//	Get pdZ (Scores)

		plhs[0] = mxCreateDoubleMatrix (n, nk, mxREAL) ;
		double *pdZ = mxGetPr (plhs[0]) ;

	//	Get pdL (Loadings)

		plhs[1] = mxCreateDoubleMatrix (p, nk, mxREAL) ;
		double *pdL = mxGetPr (plhs[1]) ;

	//	Get pdSdev (output)

		plhs[2] = mxCreateDoubleMatrix (1, nk, mxREAL) ;
		double *pdSDev = mxGetPr (plhs[2]) ;

		TRY( 
			CPCAproj (pnParamIn, pdParamIn, pdX, pdZ, pdL, pdSDev).Calc () ;
			)		
	}

	void ML_PCAprojU (int nlhs, mxArray* plhs[], int nrhs, mxArray *prhs[])
	{

		if (nlhs != 3)
			meal_error ("3 output argument expected") ;

		if (nrhs != 3)
			meal_error ("3 input arguments expected") ;

	//	Get pnParIn
		if (mxGetM (prhs[0]) * mxGetN (prhs[0]) != 6)
			meal_error ("pnParamIn: vector of length 4 expected.") ;

		double *pdnParamIn = mxGetPr (prhs[0]) ;
		int pnParamIn[8] ;
		Double2Int (pnParamIn + 2, pdnParamIn, 6) ;

		int nk = pnParamIn[3] ;

	//	Get pdParIn

		if (mxGetM (prhs[1]) * mxGetN (prhs[1]) != 1)
			meal_error ("pdParIn: vector of length 1 expected.") ;

		double *pdParamIn = mxGetPr (prhs[1]) ;

	//	Get pdX

		int &n = pnParamIn[0], &p = pnParamIn[1] ;

		n = mxGetM (prhs[2]) ;
		p = mxGetN (prhs[2]) ;

		double *pdX = mxGetPr (prhs[2]) ;

	//	Get pdZ (Scores)

		plhs[0] = mxCreateDoubleMatrix (n, nk, mxREAL) ;
		double *pdZ = mxGetPr (plhs[0]) ;

	//	Get pdL (Loadings)

		plhs[1] = mxCreateDoubleMatrix (p, nk, mxREAL) ;
		double *pdL = mxGetPr (plhs[1]) ;

	//	Get pdSdev (output)

		plhs[2] = mxCreateDoubleMatrix (1, nk, mxREAL) ;
		double *pdSDev = mxGetPr (plhs[2]) ;

		TRY( 
			CPCAprojU (pnParamIn, pdParamIn, pdX, pdZ, pdL, pdSDev).Calc () ;
			)		

	}

#endif	//	#ifdef MATLAB_MEX_FILE

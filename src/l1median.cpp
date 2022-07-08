#ifdef R_PACKAGE_FILE	//	this is an R specific source file, using built-in optimizers..

#include "R_package.h"

#include "R.h"
#include "R_ext/Applic.h"


#ifdef ES_DEV_ENV
	#include "../../../RDev/perftimer.h"
#else
	#include "perftimer.h"
#endif

	void Hess (int p, int n, double *pdX, double *pdMu, double *pdHess, double *pdTempP1, double *pdTempP2) ;	

	void ResetVect (double *p, int n, double v)						//	2do: replace this function by SVec - operation
	{
		double * const pEnd = p + n ;
		for (; p < pEnd; ++p )
			*p = v ;
	}

	void VectorMultVector (double *pA, double const * pB, int n)	//	2do: replace this function by SVec - operation
	{
		double * const pEndA = pA + n ;
		while (pA < pEndA)
		{
			*pA *= *pB ;
			++pA ;
			++pB ;
		}
	}

	inline double AddSqr (double &d, double pAdd)
	{
		return d += pAdd * pAdd ; 
	}

	class L1MinStruct									//2do: move to header file!
	{
	public:
		L1MinStruct (int n, int p, double *pdX, double *pdParScale = 0) ;
		~L1MinStruct () ;

		int m_n, m_p, m_pn ;
		double *m_pdX, *m_pdX_, *m_pdDi, *m_pdM, *m_pdParScale ; //, *m_pdM ;

		double calcall (double *pdM, double *pdMRet) ;
		double calObj (double *pdM) ;
		//void sqrtrowsumssq () ;

		int m_nCFn, m_nCGr ;

	} ;
	
	class L1MinStruct_Hess : public L1MinStruct			//2do: move to header file!
	{
	public:
		L1MinStruct_Hess (int n, int p, double *pdX, double *pdParScale = 0) :
		L1MinStruct (n, p, pdX, pdParScale),
			m_pdTemp1 (new double [p]), m_pdTemp2 (new double [p]) { }

		void calcHess (double *pdM, double *pdHess) ;

		~L1MinStruct_Hess ()
		{
			delete [] m_pdTemp1 ;
			delete [] m_pdTemp2 ;
		}

		double *m_pdTemp1, *m_pdTemp2 ;

	} ;
	

	L1MinStruct::L1MinStruct (int n, int p, double *pdX, double *pdParScale) : m_n (n), m_p (p), m_pn (p * n), m_pdX (pdX), m_pdParScale (pdParScale)
	{
		m_pdX_ = new double [n * p] ;
		m_pdDi = new double [n] ;
		m_pdM = new double [p] ;	//	set m_pdM to some strange value!!
		m_nCFn = m_nCGr = 0 ;
	}

	L1MinStruct::~L1MinStruct ()
	{
		delete [] m_pdX_ ;
		delete [] m_pdDi ;
		delete [] m_pdM ;
	}

/*	void Rtprintf (double *pd, int n)
	{
		for (; n; n--)
			Rprintf ("%0.10f\t", *pd++) ;
		Rprintf ("\r\n") ;
	}
*/
	void L1MinStruct_Hess ::calcHess (double *pdM, double *pdHess)
	{
		Hess (m_p, m_n, m_pdX, pdM, pdHess, m_pdTemp1, m_pdTemp2) ;
	}
	
	double L1MinStruct::calcall (double *pdM, double *pdMRet)
	{
		m_nCGr++ ;

//Rprintf ("grad:\t") ;
//Rtprintf (pdM, m_p) ;
		int r, c ;
		double *pdX = m_pdX + m_pn, *pdX_ = m_pdX_ + m_pn ;

		for (r = m_n - 1; r != -1; r--)
			m_pdDi[r] = 0 ;

        //    X. <- centr(X,m)

        //   d.i <- rowSums(X.^2)

		for (c = m_p - 1; c != -1; c--)
		{
			double dXm = m_pdParScale ? pdM[c] * m_pdParScale[c] : pdM[c] ;
			//for (r = 0; r < m_n; r++)
			for (r = m_n - 1; r != -1; r--)
			{
				*--pdX_ = *--pdX - dXm ;
				//*pdX_++ = *pdX++ - dXm ;
				m_pdDi[r] += *pdX_ * *pdX_ ;
			}
		}

		pdX_ = m_pdX_ + m_pn ;

		for (r = m_n - 1; r != -1; --r)
		{
			m_pdDi[r] = ::sqrt ((double) m_pdDi[r]) ;

		}

		//	MED <- colSums(X. / sqrt (d.i))


		for (c = m_p - 1; c != -1; c--)
		{
			double &dXm = pdMRet[c] = 0 ;
			for (r = m_n - 1; r != -1; r--)
				dXm -= *--pdX_ / m_pdDi[r] ;
		}

		return 0 ;
	}

	double L1MinStruct::calObj (double *pdM)
	{
		m_nCFn++ ;
//Rprintf ("obj: \t") ;
//Rtprintf (pdM, m_p) ;

		memcpy (m_pdM, pdM, m_p * sizeof (double)) ;

		if (m_pdParScale)
			VectorMultVector (m_pdM, m_pdParScale, m_p) ;

		//double *pdCur = m_pdX +  m_n * m_p ;

		double dSum = 0, dRowSum ;
		int r, c ;
		for (r = m_n - 1; r != -1; r--)
		{
			dRowSum = 0 ;
			for (c = m_p - 1; c != -1; c--)
				AddSqr (dRowSum, m_pdX [r + c * m_n] - m_pdM [c]) ;

			dSum += ::sqrt (dRowSum) ;
		}
//Rprintf ("Obj. val: %f\r\n", dSum) ;
		return dSum ;
	}

	void l1obj_hess (int n, int p, double *pdMu, double *pdHess, void *pDat)
	{
		((L1MinStruct_Hess *)pDat)->calcHess (pdMu, pdHess) ;
	}

	void l1objg(int n, double *pdCurCenter, double *pdRetCenter, void *pDat)
	{
		((L1MinStruct *) pDat)->calcall (pdCurCenter, pdRetCenter) ;
	}

	void l1obj_(int n, double *pdCurCenter, double *ddObj, void *pDat)
	{
		*ddObj = ((L1MinStruct *) pDat)->calObj (pdCurCenter) ;
	}

	double l1obj (int n, double *pdCurCenter, void *pDat)
	{
		return ((L1MinStruct *) pDat)->calObj (pdCurCenter) ;
	}

	void C_l1median_NM(int *pnParam, double *pdParam, double *pdData/*, double *pdParScale*/, double *pdMRet)
	{
		int &n = pnParam [0], &p = pnParam [1], &nMaxStep = pnParam[2], &nFail = pnParam [3], &nTrace = pnParam [4], &nFnCount = pnParam [5], &nTime = pnParam[6] ;
		double &dAbsTol = pdParam[0], &dRelTol = pdParam[1], &dRet = pdParam[2], &dAlpha = pdParam[3], &dBeta = pdParam[4], &dGamma = pdParam[5] ;

		double *pdStart = new double [p] ;

		memcpy (pdStart, pdMRet, sizeof (double) * p) ;

//		VectorDivVector (pdStart, &p, pdParScale) ;

		L1MinStruct minstruc (n, p, pdData/*, pdParScale*/) ;

		CPerfTimer tim ;
		nmmin (p, pdStart, pdMRet, &dRet, l1obj, &nFail, dAbsTol, dRelTol, (void *) &minstruc, dAlpha, dBeta, dGamma, nTrace, &nFnCount, nMaxStep) ;
		nTime = tim.GetTimeMS () ;

//		VectorMultVector (pdMRet, &p, pdParScale) ;

		delete [] pdStart ;
	}

	void C_l1median_CG(int *pnParam, int *pnParam_Out, double *pdParam, double *pdParam_Out, double *pdData/*, double *pdParScale*/, double *pdMRet)
	{
		int &n = pnParam [0], &p = pnParam [1], &nMaxStep = pnParam[2], &nTrace = pnParam [3], nType = pnParam [4] ;
		int &nFail = pnParam_Out [0], &nFnCount = pnParam_Out [1], &nGrCount = pnParam_Out [2], &nTime = pnParam_Out[3] ;

		double &dAbsTol = pdParam[0], &dRelTol = pdParam[1] ;
		double &dRet = pdParam_Out [0] ;
		double *pdStart = new double [p] ;

		memcpy (pdStart, pdMRet, sizeof (double) * p) ;

//		VectorDivVector (pdStart, &p, pdParScale) ;

		L1MinStruct minstruc (n, p, pdData/*, pdParScale*/) ;

		CPerfTimer tim ;
		cgmin (p, pdStart, pdMRet, &dRet, l1obj, l1objg, &nFail, dAbsTol, dRelTol, (void *) &minstruc, nType, nTrace, &nFnCount, &nGrCount, nMaxStep) ;
		nTime = tim.GetTimeMS () ;

//		VectorMultVector (pdMRet, &p, pdParScale) ;

		delete [] pdStart ;
	}

	void C_l1median_BFGS (int *pnParam_In, int *pnParam_Out, double *pdParam_In, double *pdParam_Out, double *pdData/*, double *pdParScale*/, double *pdMRet)
	{
		int &n = pnParam_In [0], &p = pnParam_In [1], &nMaxStep = pnParam_In[2], &nTrace = pnParam_In [3], nReport = pnParam_In [4] ;
		int &nFail = pnParam_Out [0], &nFnCount = pnParam_Out [1], &nGrCount = pnParam_Out [2], &nTime = pnParam_Out[3] ;

		double &dAbsTol = pdParam_In[0], &dRelTol = pdParam_In[1] ;
		double &dRet = pdParam_Out [0] ;

		dAbsTol = dRelTol ;

//		VectorDivVector (pdMRet, &p, pdParScale) ;

		L1MinStruct minstruc (n, p, pdData/*, pdParScale*/) ;

		int *pnMask = new int [p] ;
		for (int i = p - 1; i != -1; i--)
			pnMask[i] = 1 ;

		CPerfTimer tim ;
		vmmin (p, pdMRet, &dRet, l1obj, l1objg, nMaxStep, nTrace, pnMask, dAbsTol, dRelTol, nReport, (void *) &minstruc, &nFnCount, &nGrCount, &nFail) ;
		nTime = tim.GetTimeMS () ;

//		VectorMultVector (pdMRet, &p, pdParScale) ;

		delete [] pnMask ;
	}

	void l1median_SA (int *pnParam_In, int *pnParam_Out, double *pdParam_In, double *pdParam_Out, double *pdData/*, double *pdParScale*/, double *pdMRet)
	{
		int &n = pnParam_In [0], &p = pnParam_In [1], &nMaxStep = pnParam_In[2], &nTrace = pnParam_In [3], nTMax = pnParam_In [4] ;
		int &nFnCount = pnParam_Out [0], &nTime = pnParam_Out[1] ;

		double &dTempInit = pdParam_In [0] ;
		double &dRet = pdParam_Out [0] ;

//		VectorDivVector (pdMRet, &p, pdParScale) ;

		L1MinStruct minstruc (n, p, pdData/*, pdParScale*/) ;

		CPerfTimer tim ;
		samin (p, pdMRet, &dRet, l1obj, nMaxStep, nTMax, dTempInit, nTrace, (void *) &minstruc) ;
		nTime = tim.GetTimeMS () ;

		nFnCount = minstruc.m_nCFn ;

//		VectorMultVector (pdMRet, &p, pdParScale) ;
	}
	
	void C_l1median_NLM (int *pnParam, double *pdParam, double *pdData, double *pdMRet/*, double *pdTypSize*/)
	{
		int &n = pnParam [0], &p = pnParam [1], &nMaxStep = pnParam[2], &nFail = pnParam [3], &nTime = pnParam[5], &nMsg = pnParam[6]/*, &nTrace = pnParam[7]*/ ;
		double &dTol = pdParam[0], &dRet = pdParam[1] ;

		double *pdStart = new double [p] ;

		memcpy (pdStart, pdMRet, sizeof (double) * p) ;

		L1MinStruct minstruc (n, p, pdData) ;

//		nMsg = 8 ;

		double	*pdRetGrad = new double [p],
				*pdWrk = new double [p * 8],
				*pdA = new double [p * p],
				*pdTypSize = new double [p] ;

		
		ResetVect (pdTypSize, p, 1) ;
		
		CPerfTimer tim ;
		optif9 (	p,
					p,
					pdStart,
					l1obj_,
					l1objg,
					NULL,
					(void *) &minstruc,
					pdTypSize,
					1,			/*	fscale	*/
					1,			/*	method	*/
					1,			/*	iexp	*/
					&nMsg,
					-1,			/*	ndigit	*/	//	// C: -1 / R: 12
					nMaxStep,
					1,			/*	iagflag	*/
					0,			/*	iahflag	*/
					1e-6,		/*	dlt		*/		// 1, //1e-6,
					dTol, 		/*	gradtl	*/		//pow(DBL_EPSILON, 0.25), //0.1,
					1000,		/*	stepmx	*/		//	C:1e-6 / R: 1000
					dTol,		/*	steptl	*/		//sqrt(DBL_EPSILON), //1e-6,			
					pdMRet,	
					&dRet,
					pdRetGrad,
					&nFail,
					pdA,
					pdWrk,
					&nMaxStep) ;

		nTime = tim.GetTimeMS () ;

		delete [] pdStart ;
		delete [] pdRetGrad ;
		delete [] pdWrk ;
		delete [] pdA ;
		delete [] pdTypSize ;
	}


	void C_l1median_NLM_Hess (int *pnParam, double *pdParam, double *pdData, double *pdMRet/*, double *pdTypSize*/)
	{
		int &n = pnParam [0], &p = pnParam [1], &nMaxStep = pnParam[2], &nFail = pnParam [3], nMethod = pnParam[4], &nTime = pnParam[5], &nMsg = pnParam[6]/*, &nTrace = pnParam[7]*/, &nGFlag = pnParam[8], &nHFlag = pnParam[9], &nExp = pnParam[10], &nDigits = pnParam[11] ;
		double &dTol = pdParam[0], &dRet = pdParam[1] ;

		double *pdStart = new double [p] ;

		memcpy (pdStart, pdMRet, sizeof (double) * p) ;

		L1MinStruct_Hess minstruc (n, p, pdData) ;

//		nMsg = 8 ;

		double	*pdRetGrad = new double [p],
				*pdWrk = new double [p * 8],
				*pdA = new double [p * p],
				*pdTypSize = new double [p] ;

		ResetVect (pdTypSize, p, 1) ;

		CPerfTimer tim ;
		optif9 (	p,
					p,
					pdStart,
					l1obj_,
					l1objg,
					l1obj_hess,
					(void *) &minstruc,
					pdTypSize,
					1,			/*	fscale	*/
					nMethod,	/*	method	*/
					nExp,		/*	iexp	*/
					&nMsg,
					nDigits,	/*	ndigit	*/		//	C: -1 / R: 12
					nMaxStep,
					nGFlag,		/*	iagflag	*/
					nHFlag,		/*	iahflag	*/
					-1,		/*	dlt		*/			//	1, //1e-6,
					dTol, 		/*	gradtl	*/		//	pow(DBL_EPSILON, 0.25), //0.1,
					1000,		/*	stepmx	*/		//	C:1e-6 / R: 1000
					dTol,		/*	steptl	*/		//	sqrt(DBL_EPSILON), //1e-6,			
					pdMRet,	
					&dRet,
					pdRetGrad,
					&nFail,
					pdA,
					pdWrk,
					&nMaxStep) ;

		nTime = tim.GetTimeMS () ;

		delete [] pdStart ;
		delete [] pdRetGrad ;
		delete [] pdWrk ;
		delete [] pdA ;
		delete [] pdTypSize ;
	}
	
#endif	//	#ifdef R_PACKAGE_FILE

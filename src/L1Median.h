
#ifdef ES_DEV_ENV
	#include "../../../SMat/smat.h"
#else
	#include "smat.h"
#endif

	int l1median_HoCr (const SCMatD &mX, const SVecD &vdMedian, double dZeroTol, double dTol, int dwMaxit, int nTrace, int *pdwIterCount = NULL) ;

	class CL1Median_VZ
	{
	public:
		CL1Median_VZ (int *pnParIn, int *pnParOut, double *pdParIn, double *pdDat, double *pdMed, double *pdWeights = NULL) ;

		CL1Median_VZ (int n, int p, int &nCode, int &nIter, double *pdParIn, double *pdX, double *pdMed, double *pdWeights = NULL) ;

		BOOL Iter () ;

		t_size CheckRowSums (const double &dThreshold) ;

	protected:

		void Calc (double *pdWeights) ;

		t_size m_dwN, m_dwP, m_dwMaxIt, m_dwUseWeights ;			//	t_size input parameters 
		int m_nTrace ;												//	int input parameters
		int &m_nRetCode, &m_nIter ;									//	int output parameters

		double &m_dTol, &m_dZeroTol ;								//	double input parameters

		const t_size m_dwNHalf ;
		int m_nEqs ;

		SMatD m_mX, m_mXc ;
		SVecD m_vMed, m_vRt, m_vTt, m_vOldMed ;
		SVecD m_vWeights, m_vRowSums, m_vTemp ;
		SVecN m_mIsZero ;

	//	User Operators
	public:
		class AaCmD_BpaAmA			{ CALC_4_2(void) { a = c - d; b += sm_sqr (a) ; } } ;
		class if_C_ApaBdD			{ CALC_4_1(void) { if (c) a += b / d ; } } ;
		class if_C_Apa_inv_b		{ CALC_3_1(void) { if (c) a += 1 / b ; } } ;
		class Apa_abs_c_Bpa_abs_DmC	{ CALC_4_2(void) { a += fabs (c) ; b += fabs (d - c) ; } };
	} ;

	double calObj (const double *pdData, const double *pdM, int n, int p) ;

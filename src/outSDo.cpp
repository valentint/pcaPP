#include "outSDo.h"
#include "L1Median.h"

	double ApplyCenterMethod (const SCVecD &v, t_size dwMethod)
	{
		switch (dwMethod)
		{
		case 0:
			return 0 ;
		case 1:
			return mean (v) ;
		case 2:
			return median (v) ;
		}

		ASSERT (0) ;
		return 0 ;
	}

#define IM_OBS		0
#define IM_DIFFOBS	1
#define IM_RANDOM	2
#define IM_RANDDIFFOBS 3

	CSDoOut::CSDoOut (int *pnParIn, double *pdX, double *pdMaxMaha, int *pnNChanged)
		: m_dwN (pnParIn[0]), m_dwP (pnParIn[1]), m_dwIterMethod (pnParIn[2]), m_dwIterParam (pnParIn[3]), m_dwCenterMethod (pnParIn[4]), m_dwScatterMethod (pnParIn[5]), m_dwReset (pnParIn[6])
		, m_mX (pdX, m_dwN, m_dwP)
		, m_vMaxMaha (pdMaxMaha, m_dwN)
		, m_dwNDir ((m_dwIterMethod == IM_OBS) ? m_dwN : m_dwIterParam)	//	number of checked directions
		, m_vXProj (m_dwN), m_vCurDir (m_dwP)
		, m_pnNChanged (pnNChanged), m_pdDiff (pdX)
		, m_pdXProj (m_vXProj), m_pdEndXProj (m_vXProj.GetDataEnd ()), m_pdMaxMaha (m_vMaxMaha)
	{
	}

	void CSDoOut::Calc ()
	{
		if (m_dwReset)
			m_vMaxMaha.Reset (0) ;

		switch (m_dwIterMethod)
		{
		case IM_OBS :
			IterObs () ;
			break ;
		case IM_DIFFOBS :
			IterDiffObs (m_dwNDir) ;
			break ;
		case IM_RANDOM :
			IterRand (m_dwNDir) ;
			break ;
		case IM_RANDDIFFOBS :
			IterRandDiffObs (m_dwNDir) ;
			break ;
		}
	}

	void CSDoOut::IterObs ()
	{
		t_size i ;
		double dNorm ;
		int nChanged ;
//		for (i = m_dwN - 1; i != NAI; i)	//	can't be used because of m_vChanged
		for (i = 0; i < m_dwN; ++i)
		{
			CopyRow (*m_vCurDir, m_mX, i) ;

			EO<SOP::Apa_sqr_B>::SVc (dNorm = 0, m_vCurDir) ;
			EO<SOP::a_divide>::VSc (*m_vCurDir, sqrt (dNorm)) ;

			nChanged = DoDir (m_vCurDir) ;
			if (m_pnNChanged)
				m_pnNChanged[i] = nChanged ;
//			if (m_bSaveChanged)
//				m_vChanged(i) = nChanged ;
		}
	}

	void CSDoOut::IterDiffObs (int n)
	{
		
	}

	void CSDoOut::IterRandDiffObs (int n)
	{
		
	}

	void CSDoOut::IterRand (int n)
	{
		
	}

	int CSDoOut::DoDir (const SCVecD &vLoad)
	{
		m_vXProj.Reset (0) ;
		EO<SOP::ApaBmC>::VMcVct (*m_vXProj, m_mX, vLoad) ;

		const double dCenter = ApplyCenterMethod (m_vXProj, m_dwCenterMethod) ;
		const double dScatter = ApplyMethod (m_vXProj, m_dwScatterMethod) ;

		double *pProj = m_pdXProj, *pdMaxMaha = m_pdMaxMaha ;

		int nChanged = 0 ;

		double dDiff = 0, dCurMaha ;

		while (pProj < m_pdEndXProj)
		{
			//nChanged += sm_setmax_b  (*pdMaxMaha, (*pProj - dCenter) / dScatter) ;

			dCurMaha = fabs(*pProj - dCenter) / dScatter ;

			if (dCurMaha > *pdMaxMaha)
			{
				*pdMaxMaha = dCurMaha ;
				++nChanged ;
				dDiff += dCurMaha - *pdMaxMaha ;
			}

			++ pProj ;
			++ pdMaxMaha ;
		}

		return nChanged ;
	}

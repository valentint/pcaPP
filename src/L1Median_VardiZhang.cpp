#include "L1Median.h"

	t_size CL1Median_VZ::CheckRowSums (const double &dThreshold)
	{		//	counts the elements in array m_vRowSums which are greater than dThreshold
		t_size dwRet = 0 ;

		ASSERT (m_mIsZero.size () == m_vRowSums.size ()) ;

		const double *pRS = m_vRowSums, * const pEndRS = m_vRowSums.GetDataEnd () ;
		int *pZero = m_mIsZero ;

		while (pRS < pEndRS)
		{
			*pZero = *pRS > dThreshold ;
			if (*pZero)
				dwRet ++ ;
			++pRS ;
			++pZero ;
		}
		
		return m_dwN - dwRet ;	
	}

	BOOL CL1Median_VZ::Iter ()
	{
		m_vRowSums.Reset (0) ;
		EO<AaCmD_BpaAmA>::MVMcVct (!m_mXc, *m_vRowSums, m_mX, m_vMed) ;
		EO<SOP::a_sqrt>::V (*m_vRowSums) ;

		double dMin = min (m_vRowSums) ;

		t_size dwGreater = 0 ;
		EO<SOP::inc_a_if_b_lesseq_c>::SScVc (dwGreater, dMin / m_dZeroTol, m_vRowSums) ;

		if (dwGreater * 2 > m_dwN)
		{	//	some of the elements of m_vRowSums are zero
			m_nEqs++ ;

			t_size dwZero = CheckRowSums (median (m_vRowSums) * m_dZeroTol) ;

			if (dwZero > m_dwNHalf)	//	there's more than half of the observations concentrated in one single point
			{
				if (m_nTrace >= 1)
					meal_printf ("%d >= n / 2 = %d observations concentrated in one point found.\r\n", dwZero) ;
				return FALSE ;
			}

			if (m_nTrace >= 1)
				meal_printf ("%d observations are exatly at the median.\r\n", dwZero) ;
			if (m_nTrace >= 0 && dwZero > 1)
				meal_warning ("The current L1median estimate is ident with more than one observation. The resulting l1median estimation might be incorrect. [CL1Median_VZ::Iter]") ;

			m_vTt.Reset (0) ;
			EO<if_C_ApaBdD>::VtMcVcVc_NC (*m_vTt, m_mX, m_mIsZero, m_vRowSums) ;

			m_vRt.Reset (0) ;
			EO<if_C_ApaBdD>::VtMcVcVc_NC (*m_vRt, m_mXc, m_mIsZero, m_vRowSums) ;

			double dTemp = 0 ;
			EO<if_C_Apa_inv_b>::SVcVc (dTemp, m_vRowSums, m_mIsZero) ;
			EO<SOP::a_divide>::VSc (*m_vTt, dTemp) ;

			dTemp = 0 ;
			EO<SOP::Apa_sqr_B>::SVc (dTemp, m_vRt) ;
			double edivr = dwZero / sqrt (dTemp) ;

			if (edivr > 1)
				EO<SOP::a_multiply>::VSc (*m_vMed, edivr) ;

			if (edivr < 1)
				EO<SOP::ApaBmC>::VScVc (*m_vMed, 1 - edivr, m_vTt) ;
		}
		else
		{
			m_vMed.Reset (0) ;
			EO<SOP::ApaBdC>::VtMcVc_NC (*m_vMed,m_mX, m_vRowSums) ;

			double dSum = 0 ;
			EO<SOP::Apa1dB>::SVc (dSum, m_vRowSums) ;

			EO<SOP::a_divide>::VSc (*m_vMed, dSum) ;
		}
		return TRUE ;
 	}

	CL1Median_VZ::CL1Median_VZ (int *pnParIn, int *pnParOut, double *pdParIn, double *pdX, double *pdMed, double *pdWeights)
		: m_dwN (pnParIn[0]), m_dwP (pnParIn[1]), m_dwMaxIt (pnParIn[2]), m_dwUseWeights (pnParIn[3])
		, m_nTrace (pnParIn[4])
		, m_nRetCode (pnParOut[0]), m_nIter (pnParOut[1])

		, m_dTol (pdParIn[0]), m_dZeroTol (pdParIn [1])

		, m_dwNHalf (m_dwN >> 1)
		, m_nEqs (0)

		, m_mX (pdX, m_dwN, m_dwP), m_mXc (m_dwN, m_dwP)
		, m_vMed (pdMed, m_dwP), m_vRt (m_dwP), m_vTt (m_dwP), m_vOldMed (m_dwP)
		//, m_vWeights (pdWeights, m_dwN)
		, m_vRowSums (m_dwN), m_vTemp (m_dwN)
		, m_mIsZero (m_dwN)
	{ 
		Calc (pdWeights) ;
	}

	CL1Median_VZ::CL1Median_VZ (int n, int p, int &nCode, int &nIter, double *pdParIn, double *pdX, double *pdMed, double *pdWeights)
		: m_dwN (n), m_dwP (p), m_dwMaxIt ((int) pdParIn[0]), m_dwUseWeights (0)
		, m_nTrace ((int) pdParIn[1])
		, m_nRetCode (nCode), m_nIter (nIter)

		, m_dTol (pdParIn[2]), m_dZeroTol (pdParIn [3])

		, m_dwNHalf (m_dwN >> 1)
		, m_nEqs (0)

		, m_mX (pdX, m_dwN, m_dwP), m_mXc (m_dwN, m_dwP)
		, m_vMed (pdMed, m_dwP), m_vRt (m_dwP), m_vTt (m_dwP), m_vOldMed (m_dwP)
		//, m_vWeights (pdWeights, m_dwN)
		, m_vRowSums (m_dwN), m_vTemp (m_dwN)
		, m_mIsZero (m_dwN)
	{ 
		Calc (pdWeights) ;
	}

	void CL1Median_VZ::Calc (double *pdWeights)
	{ 

		if (pdWeights)
			m_vWeights.Set (pdWeights, m_dwN) ;

		t_size i ;

		for (i = m_dwMaxIt - 1; i != NAI; i--)
		{
			m_vOldMed.Copy (m_vMed) ;
			if (!Iter ())
				break ;

			double dAbsDiff = 0, dAbsSum = 0 ; 
			EO<Apa_abs_c_Bpa_abs_DmC>::SSVcVc_NC (dAbsSum, dAbsDiff, m_vMed, m_vOldMed) ;

			if (m_nTrace >= 2)
			{
				if (m_nTrace >= 3)
				{
					meal_printf ("k=%3d rel.chg=%12.15g, m=(",m_dwMaxIt - i, dAbsDiff/dAbsSum),
					meal_printf (")\n") ;
				}
				else
					meal_printf (".") ;
			}

			if (dAbsDiff < m_dTol * dAbsSum)
				break ;
		}

        if(m_nTrace)
			meal_printf (" needed %d iterations (%d of which had y==x_k)\r\n", m_dwMaxIt - i, m_nEqs) ;

		m_nIter = m_dwMaxIt - i ;
	}

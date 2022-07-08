#include "PCAgrid.h"


#define CHECK_ORTH
#define GLOBAL_POWER	2	//	use 2 for maximizing (robust) variances, 1 for maximizing (robust) standard deviations.


#if GLOBAL_POWER == 1
		double ngpf (const double &d) {return d ; }
#else 
	#if GLOBAL_POWER == 2
		double ngpf (const double &d) {return d * d ; }
	//	#else: error, because gpf not defined.
	#endif
#endif



/////////////////
//	CsPCAGrid  //
/////////////////


	CsPCAGrid::CsPCAGrid (int *pnParamIn, int *pnParamOut, double *pdParamIn, double *pdData, double *pdLoadings, double *pdSDev, double *pdObj/*, double *pdMaxMaha*/, double *pdLambda, double *pdBackTransHD)
		:	CPCAGrid (pnParamIn, pnParamOut, pdParamIn, pdData, pdLoadings, pdSDev, pdObj/*, pdMaxMaha*/)
		,	m_nGloScatter (pnParamIn[9]), m_nSpeedUp (pnParamIn[11]), m_dwPHD (pnParamIn[10])
		,	m_dQ (pdParamIn[1]), m_dS (pdParamIn[2])
		,	m_bUseQ (m_dQ != 1.0), m_bUseS (m_dS != 1.0)
//		,	m_bUseQ (TRUE), m_bUseS (TRUE)
		,	m_vLambda (pdLambda, m_dwK), m_vTempP (m_dwP), m_vTempPSub (m_dwP)
		,	m_dGloScatter (1)
	{
		if (m_dwPHD)															//	always false for sPCAgrid
		{
			m_mBackTransHD.Set (pdBackTransHD, m_dwPHD, m_dwP) ;
			m_mBackProj.Require (m_dwP, m_dwPHD) ;
			m_vLoadHD.Require (m_dwPHD) ;
			m_vSumLoadOthers.Require (m_dwPHD) ;
		}
		else
		{
			m_mBackProj.Require (m_dwP, m_dwP) ;
			m_vSumLoadOthers.Require (m_dwP) ;
		}
		if (!m_dwCheckOrth && m_nGloScatter == 0)
			m_dGloScatter = ngpf (ApplyMethodMean (m_mX)) ;
//		else if (m_nGloScatter == -1)
//			m_dGloScatter = 1 ;
	}

	void CsPCAGrid::OnCalcPC ()
	{
		if (!m_dwCheckOrth && m_nGloScatter == 1)
			m_dGloScatter = ngpf (ApplyMethodMean (TempY ())) ;
//meal_printf ("gloscat %f\r\n", m_dGloScatter) ;

		m_vTempPSub.Reshape (m_dwPSub) ;
		m_dCurLambda = m_vLambda (m_dwCurK - m_dwkIni) ;

//meal_printf ("Glo scat for dim %d = %f\r\n", m_dwCurK, sqrt (m_dGloScatter)) ;

		// calc new backtransformation matrix.

		if (m_dwPHD)															//	always false sPCAgrid
			sme_matmult_R (m_mBackTransHD, m_mL.GetColRef (m_dwCurK, m_dwP), !m_mBackProj) ;
		else
			m_mBackProj.Copy_R (m_mL.GetColRef (m_dwCurK, m_dwP)) ;				//	2do: use a reference!
	}

	void CsPCAGrid::InitPenalty ()
	{
		m_vSumLoadOthers.Reset (0) ;
		EO<SOP::ApaBmC>::VMcVct (*m_vSumLoadOthers, m_mBackProj, m_vAfin) ;

		m_vSumLoadThis.Copy_R (m_mBackProj.GetColRef (m_dwCurP)) ;				//	2do: try simple assignment or use this column directly where m_vSumLoadThis is needed!
	}

																				//	computes the objective function of a direction
	double CsPCAGrid::CalcObj ( const double dCos, const double dSin			//	the direction to be examined
								, double &dScat									//	(output) the scatter of a matrix cbind (m_vSumLoadThis, m_vSumLoadOthers) projected onto a direction (dCos, dSin)
								, double &dScatOrth)							//	(output) the scatter of a matrix cbind (m_vSumLoadThis, m_vSumLoadOthers) projected onto a direction (dSin, dCos)
																				//			(-> maybe dSin and dCos are mixed up in this comment)
	{
		dScat = CalcProjScat (dCos, dSin) ;

		if (m_dwCheckOrth)														//	don't care about the experimental "m_dwCheckOrth" - option so far
		{																		//	it's always off for sPCAgrid
//			CalcMaha (dScat) ;


			dScatOrth = CalcProjScat (-dSin, dCos) ;
//			CalcMaha (dScatOrth) ;

			return ngpf (dScat / dScatOrth) + GetPenalty (dCos, dSin); ;

////			dCurObj /= ngpf (dScatOrth + dCurScat * 0.001) ;
////			dCurObj /= (m_dGloScatter + ngpf (dScatOrth)) ;
		}
		return																	//	the objective function's value:
			  ngpf (dScat)														//	the scatter (or variance) of the data projected onto the current direction
																				//	ngpf(x) returns either x, or x^2, depending on the definition of GLOBAL_POWER
																				//	NOTE: GLOBAL_POWER is now deprecated, as a "q" and an "s" parameter are available in "GetPenalty", which do provide the same functionality

			+ GetPenalty (dCos, dSin)											//	+ the penalty function for the current direction
			* m_dGloScatter;													//	* some normalization factor for the overall variance of the data in the currently considered subspace
	}

																				//	computes the sparseness - penalty term for a direction
	double CsPCAGrid::GetPenalty (const double& dCos, const double& dSin)		//	dCos & dSin specify the direction to be checked
	{
/*		//	computing the Sparseness penalty:

		r = c (dSin, dCos) ;
		z = cbind (m_vSumLoadOthers, m_vSumLoadThis)

		return sum ((abs (z %*% r))^q)^p ;	// the sparseness - criterion

*/

		ASSERT (!m_nSpeedUp) ;	//	the old penalty stuff has not been implemented here...

		if (m_dCurLambda == 0)
			return 0 ;

		double dRet = 0 ;

		if (m_bUseQ)															//	the q-parameter for the norm is != 1
		{
			const double adParams [] = {dCos, dSin, m_dQ} ;

			if (fabs (dCos) <= m_dZeroTol)										//	the value of dCos is zero
				EO<UOP::Apa_pow_abs_C_B>::SScVc (dRet, m_dQ, m_vSumLoadThis) ;	//	specific (faster) compuation for dCos == 0
																				//	computes dRet = sum (abs (m_vSumLoadOthers) ^ m_dQ)
/*		DESCRIPTION of matrix computation:

                                       A     B       C
    EO<UOP::Apa_pow_abs_C_B>::SScVc (dRet, m_dQ, m_vSumLoadThis)
    Apa_pow_abs_C_B     <==>     A += pow (abs (C), B)		 (for each value of C)

	pa  == plus assign (+=)
	p = plus
	m = multiply
	d = divid
	s = subtract

*/
			else if (fabs (dSin) <= m_dZeroTol)									//	the value of dSin is zero
				EO<UOP::Apa_pow_abs_C_B>::SScVc (dRet, m_dQ, m_vSumLoadOthers) ;//	specific (faster) computation for dSim == 0
																				//	computes dRet = sum (abs (m_vSumLoadOthers) ^ m_dQ)

			else																//	dCos != 0 && dSin != 0 -> no simple computation possible
				//EO<UOP::Apa_abs_BmDpCmE_>::SScScVcVc_NC (dRet, dCos, dSin, m_vSumLoadOthers, m_vSumLoadThis) ;
				EO<UOP::Apa_pow_abs_B0mCpb1mD_B2>::SSVcVc_NC (dRet, adParams, m_vSumLoadOthers, m_vSumLoadThis) ;
																				//	computes dRet = sum (abs (dCos * m_vSumLoadOthers + dSin * m_vSumLoadThis) ^ m_dQ)
		}
		else
		{
			if (fabs (dCos) <= m_dZeroTol)
				EO<UOP::Apa_abs_B>::SVc (dRet, m_vSumLoadThis) ;
																				//	computes dRet = sum (abs (m_vSumLoadThis))
			else if (fabs (dSin) <= m_dZeroTol)
				EO<UOP::Apa_abs_B>::SVc (dRet, m_vSumLoadOthers) ;
			else																//	computes dRet = sum (abs (m_vSumLoadOthers))
				EO<UOP::Apa_abs_BmDpCmE_>::SScScVcVc_NC (dRet, dCos, dSin, m_vSumLoadOthers, m_vSumLoadThis) ;
																				//	computes dRet = sum (abs (dCos * m_vSumLoadOthers + dSin * m_vSumLoadThis))
		}

		if (m_bUseS)															//	the s - parameter (for the norm) is != 1
			dRet = pow (dRet, m_dS) ;

		return - dRet * m_dCurLambda ;

		


//		dRet *= - m_dCurLambda ;
//		if (!m_dwCheckOrth)
//			dRet *= m_dGloScatter ;
//
//		return  dRet ;
	}

////////////////
//	CPCAGrid  //
////////////////
						

	CPCAGrid::CPCAGrid (int *pnParamIn, int *pnParamOut, double *pdParamIn, double *pdData, double *pdLoadings, double *pdSDev, double *pdObj/*, double *pdMaxMaha*/)
		: m_dwN (pnParamIn[0]), m_dwP (pnParamIn[1]), m_dwK (pnParamIn[2]), m_dwSplitCircle (pnParamIn[3]), m_dwMaxIter (pnParamIn[4]), m_dwMethod (pnParamIn[5]), m_dwTrace (pnParamIn[6]), m_dwkIni (pnParamIn[7]), m_dwCheckOrth (pnParamIn[8])
		, m_nReturn (pnParamOut [0])
		, m_dZeroTol (pdParamIn[0])
		, m_mX (pdData, m_dwN, m_dwP), m_mL (pdLoadings, m_dwP, m_dwP)//, m_mTempPP(m_dwP, m_dwP), m_mTempPN(m_dwN, m_dwP)
		, m_vAfin(m_dwP), m_vAfinBest (m_dwP), m_vScl(m_dwP), m_vYOpt (m_dwN), m_vSDev(pdSDev, m_dwP), m_vObj (pdObj, m_dwK)
		, m_vProj (m_dwN)//, m_vMaxMaha (pdMaxMaha, m_dwN)
		, m_vOrd (m_dwP)
		, m_dwCurK (0), m_dwCurP (0), m_dwPSub (0), m_dwTempYIdx (0)

		, m_pdProj (m_vProj)//, m_pdEndProj (m_vProj.GetDataEnd ())
			, m_pdCurLC (m_vYOpt), m_pdCurLCEnd (m_vYOpt.GetDataEnd ())
	{
		m_mY[0].Require (m_dwN, m_dwP) ;
		m_mY[1].Require (m_dwN, m_dwP) ;
//		m_vMaxMaha.Reset (0) ;
	}

	int CPCAGrid::Calc ()
	{
		if (m_dwK > m_dwP)
			return 1 ;		//	k > p

/*		if ((m_nSplitCircle & 1) == 0)	//
		//if (m_dwSplitCircle & 1)		//  only allow even values for splitcircle
			++ m_nSplitCircle ;*/

		if (m_dwkIni)
			sme_matmult_R (m_mX, m_mL.GetColRef (m_dwkIni, m_dwP), !TempY ()) ;
		else
		{
			TempY ().Copy (m_mX) ;
			SetDiag_sq (!m_mL) ;
			//m_mL.setdiag () ;													//	this MUST now happen in the calling R routine! // this has been changed when introducing the m_dwkIni argument // why not here?
		}

		for (m_dwCurK = m_dwkIni; m_dwCurK < m_dwK; m_dwCurK++)					//	for each PC which to be computed
		{
			m_dwPSub = m_dwP - m_dwCurK ;										//	dimensionality of the subspace

			OnCalcPC () ;

			if (m_dwPSub == 1)													//	only 1 dimension left -> return this direction
			{
				m_vSDev (m_dwCurK) = ApplyMethod (TempY ().GetColRef (0)) ;
				continue ;	//	break ;
			}

			m_vScl.Reshape (m_dwPSub) ;
			m_vOrd.Reshape (m_dwPSub) ;
			ApplyMethod (TempY (), m_vScl) ;									//	2do: m_vScl can be a temporary vector

			meal_sort_order_rev (m_vScl, m_vOrd, m_vScl.size ()) ;				//	gets the order (m_vOrd) of the dimensions regarding to their variance (m_vScl) in decreasing order

			m_dwCurP = m_vOrd(0) ;												//	index of the coloumn of x with biggest scatter

			m_vAfinBest.Reshape (m_dwPSub) ;
			m_vAfin.Reshape (m_dwPSub) ;

			m_vAfin.Reset (0) ;
			m_vAfin (m_dwCurP) = 1 ;

			CopyCol (*m_vYOpt, TempY (), m_dwCurP) ;							//	loads the loading with max scatter as the initial solution

			t_size i, j ;
			double dCurSplit ;

			double dScatBest = 0 ;

			double dObjBest = 0 ;
			for (i = 0; i <= m_dwMaxIter; i++)									//	the outer iteration, which subsequently decreases the gridsize ( = dCurSplit)
			{
				//double dScat, dObj, dSumAbsDelta = 0 ;
																				// 2do: check if it's better to use  * 0.5 each time?
				dCurSplit = pow (0.5, (double) i) ;								//	the current gridSize

				for (j = 0; j < m_dwPSub; j++)									//	for each loading in the current solution
				{
					m_dwCurP = m_vOrd (j) ;										//	m_dwCurP = the j-th largest coponent
																				//	the m_dwCurP-st loading in tje current solution is now altered in order to increase the objective function

					m_vCurY = TempY ().GetColRef (m_dwCurP) ;					//	2do: move this 2 rows down.
					m_pdCurY = m_vCurY ;

					const double dL = m_vAfin (m_dwCurP) ;						//	current loading 
					if (fabs (dL) == 1)
						continue ;
					RemoveLoading (/*i*/) ;

					m_dNL = dL ;
					GridPlane (dCurSplit) ;										//	increasing the objective function by trying several directions on a grid
					AddLoading (m_dNL, m_dNCL) ;								//	add the found loading tho the current solution

					//double dNL = dL, dNCL ;
					//GridPlane (dNL, dNCL, dScat, dObj, dCurSplit) ;
					//AddLoading (dNL, dNCL) ;
					//dSumAbsDelta += fabs (dL - dNL) ;
				}

				EO<SOP::a_divide>::VSc (*m_vAfin, norm2 (m_vAfin)) ;			//	2do: check norm of m_vAfin. should be 1 anyway!, if not it's sufficient to perform this normalization after this for loop, only with the m_vAfinBest - vector!

				if (!i || dObjBest <= m_dBestObj)								//	checking whether we've found a better solution -> if true store it.
				{
					dObjBest = m_dBestObj ;
					m_vAfinBest.Copy_NC (m_vAfin) ;
					dScatBest = m_dCurScat ;
				}

//meal_printf ("delta: %.22f ->", dSumAbsDelta) ;
/*				if (dSumAbsDelta <= m_dZeroTol)	//	no changes of any loading -> quit // doesn't make sense as we're operating on a raster
				//if (dCurSplit<= m_dZeroTol)	//	no changes of any loading -> quit
				{
//meal_printf ("stop iteration\n") ;
					if (m_dwTrace >= 3)
						meal_printf ("Calculation of PC %d stopped after %d loops\r\n", m_dwCurK + 1, i + 1) ;
					break ;
				}
//meal_printf ("continue iteration\n") ;
*/			}

			m_vSDev (m_dwCurK) = dScatBest ;		//	2do: use ptrs instead of vector access!
			m_vObj(m_dwCurK) = dObjBest ;
			BackTransform () ;
		}
		return 0 ;
	}

	void CPCAGrid::BackTransform ()
	{
		ASSERT_TEMPRANGE (0, 2) ;
		SMatD m_mTempPP (tempRef (0), m_dwPSub, m_dwPSub) ;
		SetDiag_sq (!m_mTempPP) ;

		int dwIdxRef = m_vOrd (0) ;

		set_neg (*m_vAfinBest) ;

		m_vAfinBest (dwIdxRef) += 1 ;

		double dNorm = norm2 (m_vAfinBest) ;

		if (dNorm > m_dZeroTol)
		{
			static const double dSqrt2 = sqrt ((double) 2.0) ;
			EO<SOP::a_divide>::VSc (*m_vAfinBest, dNorm / dSqrt2) ;
			EO<SOP::AsaBmC>::MVcVct (!m_mTempPP, m_vAfinBest, m_vAfinBest) ;
		}

		SMatD mProjSorted (tempRef (1), m_dwPSub, m_dwPSub) ;

		mProjSorted.CopyCol_Order_NC (m_mTempPP, *m_vOrd) ;						//	undo sorting

		SMatD mOldLoadings (tempRef (2), m_dwP, m_dwPSub) ;						//	2do: copying cols should be done in constructor (using GetColRef ()
		CopyCol (!mOldLoadings, m_mL, m_dwCurK, m_dwP) ;

		sme_matmult (mOldLoadings, mProjSorted, !m_mL.GetColRef (m_dwCurK, m_dwP)) ;
		sme_matmult_R (TempY (), mProjSorted.GetColRef (1, m_dwPSub), !TempYC ()) ;

		SwapTempY () ;
	}

	void CPCAGrid::RemoveLoading (/*int i*/)
	{
		const double &dL = m_vAfin (m_dwCurP) ;
		if (dL == 0)
			return ;

		const double dCL = sqrt (1.0 - sm_sqr (dL)) ;							//	current loading and complementary loading, such that dL^2 + dCL^2 == 1

		EO<UOP::Aa_AsDmB_dC>::VScScVc (*m_vYOpt, dL, dCL, m_vCurY) ;			//	m_vYOpt = m_vYOpt
		EO<SOP::a_divide>::VSc (*m_vAfin, dCL) ;
		m_vAfin(m_dwCurP) = 0 ;
	}

	void CPCAGrid::AddLoading (const double &dNL, const double &dNCL)
	{
		//	dNL		...	New Loading
		//	dNCL	... New Complement Loading

		EO<UOP::Aa_AmC_p_DmB>::VScScVc (*m_vYOpt, dNL, dNCL, m_vCurY) ;			//	m_vYOpt = m_vYOpt * dNCL + dNL * m_vCurY

		//	Aa_AmC_p_DmB = A = A*C + D*B
		//(*m_vYOpt, dNL, dNCL, m_vCurY)  == (A, B, C, D)

		EO<SOP::a_multiply>::VSc (*m_vAfin, dNCL) ;								//	m_vAfin = m_vAfin * dNCL
		m_vAfin (m_dwCurP) = dNL ;												//	m_vAfin[m_dwCurP] = dNL
	}

	double CPCAGrid::ApplyMethodMean (const SCMatD &m)
	{
		double dSd = 0 ;
		int i ;
		for (i = m.ncol () - 1; i != (int) -1; i--)
			dSd += sm_sqr(ApplyMethod (m.GetColRef (i))) ;
		return sqrt (dSd  / m.ncol ()) ;
	}

	void CPCAGrid::ApplyMethod (const SCMatD &m, SVecD &v)
	{
		v.Reshape (m.ncol ()) ;
		int i ;
		for (i = m.ncol () - 1; i != (int) -1; i--)
			v(i) = ApplyMethod (m.GetColRef (i)) ;
	}

	double CPCAGrid::ApplyMethod (const SCVecD &v)	//	2do: remove..
	{
		return ::ApplyMethod (v, m_dwMethod) ;
	}


																				//	projects the (2d) scores cbind (m_pdCurLC, m_pdCurY) onto the direction c(dCos, dSin) and computes a scale estimate of the projected data
	double CPCAGrid::CalcProjScat (const double dCos, const double dSin)
	{
		const double *pYOpt = m_pdCurLC, *pCurY = m_pdCurY ;
		double *pProj = m_pdProj ;

		while (pYOpt < m_pdCurLCEnd)											//	projecting the data
		{
			*pProj = *pYOpt * dCos + *pCurY  * dSin ;
			++pProj ;
			++pYOpt ;
			++pCurY ;
		}

		return ApplyMethod (m_vProj) ;											//	computing the scale estimate
	}

/*	void CPCAGrid::CalcMaha (const double dScat)
	{
		double *pProj = m_pdProj ;
		double *pdProjEnd = m_pdProj + m_dwN ;

		double *pdMaha = m_vMaxMaha ;

		while (pProj < pdProjEnd)
		{
			sm_setmax (*pdMaha, *pProj / dScat) ;
			++pdMaha ;
			++pProj ;
		}
	}
*/
	double CPCAGrid::CalcObj (const double dCos, const double dSin, double &dScat, double &dScatOrth)
	{
		dScat = CalcProjScat (dCos, dSin) ;

		if (m_dwCheckOrth)														//	always false for  PCAgrid
		{
//			CalcMaha (dScat) ;

			dScatOrth = CalcProjScat (dCos, -dSin) ;
//			CalcMaha (dScatOrth) ;

			return ngpf (dScat / dScatOrth) ;

////			dCurObj /= ngpf (dScatOrth + dCurScat * 0.001) ;
////			dCurObj /= (m_dGloScatter + ngpf (dScatOrth)) ;
		}

		return ngpf (dScat) ;													//	the PCAgrid objective function is equal to the scatter of the projected data
	}

	void CPCAGrid::EvalDirection (const double dCos, const double dSin)
	{
		double	dScat, dScatOrth ;
		double dCurObj = CalcObj (dCos, dSin, dScat, dScatOrth) ;

		if (dCurObj > m_dBestObj)
		{
			m_dBestObj = dCurObj ;
			m_dCurScat = dScat ;
			m_dCurScatOrth = dScatOrth ;

			m_dNL = dSin ;
			m_dNCL = dCos ;
		}
	}

	double CPCAGrid::CalcVarTrimmed (double dCos, double dSin, double dScat, double dScatOrth)
	{
		if (dScatOrth <= m_dZeroTol || dScat <= m_dZeroTol)
			return dScat ;

		const double *pYOpt = m_pdCurLC, *pCurY = m_pdCurY ;
		double dCurProj, dCurProjOrth ;

		dScat = 1 / dScat ;
		dScatOrth = 1 / dScatOrth ;

		double dS = 0, dSS = 0 ;
		int n = 0 ;

		while (pYOpt < m_pdCurLCEnd)
		{
			dCurProj = *pYOpt * dCos + *pCurY  * dSin ;
			dCurProjOrth = *pYOpt * dSin - *pCurY  * dCos ;

			if (sm_sqr (dCurProj) * dScat + sm_sqr (dCurProjOrth) * dScatOrth < 6)
			{
				dS += dCurProj ;
				dSS += sm_sqr (dCurProj) ;
				++n ;
			}

			++pYOpt ;
			++pCurY ;
		}

		return (dSS / n - (sm_sqr (dS / n))) * n / (n - 1.0) * 1.3178;	//	correction factor 1.3178 for 95% quantile in sqr maha distance (6)
	}

	double CPCAGrid::CalcScatTrimmed (double dCos, double dSin, double dScat, double dScatOrth)
	{
		if (dScatOrth <= m_dZeroTol || dScat <= m_dZeroTol)
			return dScat ;

		const double *pYOpt = m_pdCurLC, *pCurY = m_pdCurY ;
		double dCurProjOrth ;
		double *pProj = m_pdProj ;

		while (pYOpt < m_pdCurLCEnd)
		{
			dCurProjOrth = *pYOpt * dSin - *pCurY  * dCos ;

			if (sm_sqr (dCurProjOrth) / dScatOrth <= 3.841459)
			{
				*pProj = *pYOpt * dCos + *pCurY  * dSin ;
				++pProj ;
			}

			++pYOpt ;
			++pCurY ;
		}

//		return (dSS / n - (sm_sqr (dS / n))) * n / (n - 1.0) * 1.3178;	//	correction factor 1.3178 for 95% quantile in sqr maha distance (6)
		return ApplyMethod (SVecD (m_pdProj, pProj - m_pdProj)) ;
	}

	void CPCAGrid::GridPlane (double dCurSplit)
	{
		const double dASinNL = asin (m_dNL) ;

		const double dSm1 = (m_dwSplitCircle > 1) ? m_dwSplitCircle - 1 : 1 ;

		double dSplitFact = meal_PI () * dCurSplit ;

		InitPenalty () ;

		t_size i ;

		ASSERT_TEMPRANGE (11, 11) ;

		SVecD vProj (tempRef (11), m_dwN) ;

		m_dBestObj = meal_NegInf () ;

		if (m_dNL && fabs (m_dNL) < 1e-6)
			EvalDirection (1, 0) ;

		const t_size dwEnd = (dCurSplit == 1.0) ? m_dwSplitCircle - 1 : m_dwSplitCircle ;	//	dCurSplit means, that we're checking an angle of 180° ( = PI). thus the first and last checked point would be the same.

		for (i = 0; i < dwEnd; i++)
		{
			double dAngle = (i / dSm1 - 0.5) * dSplitFact + dASinNL;

			EvalDirection (cos (dAngle), sin (dAngle)) ;
		}
		if (m_dwCheckOrth)														//	always false for PCAgrid and sPCAgrid
			m_dCurScat = sqrt (CalcVarTrimmed (m_dNCL, m_dNL, m_dCurScat, m_dCurScatOrth)) ;

//			m_dCurScat = CalcScatTrimmed (m_dNCL, m_dNL, m_dCurScat, m_dCurScatOrth) ;
//		else
//			m_dCurScat = m_dCurScat ;
	}

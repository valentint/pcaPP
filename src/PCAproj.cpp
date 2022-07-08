#include "PCAproj.h"

	void vec_mult_mat_t_partial (double *pA, double const *pB, double const *pC, const int n, const int p, const int nDimN)
	{
		//	calculates
		//
		//	pA <- pB %*% t (pC[1:n, ])	
		//	with	pA: a result vector of length n
		//			pB: a vector of length p
		//			pC: a matrix of dimension nDimN x p
		//			with n <= nDimN

		THROW (n <= nDimN) ;
		const t_size dwJump = nDimN - n ;
		double const *const pEndC = pC + nDimN * p ;

		double * const pStartA = pA, * const pEndA = pA + n ;

		for (; pA < pEndA; ++pA)
			*pA = 0 ;

		while (pC < pEndC)
		{
			pA = pStartA ;
			while (pA < pEndA)
			{
				*pA += *pB * *pC ;
				++pA ;
				++pC ;
			}

			pC += dwJump ;
			++pB ;
		}
	}

	CPCAproj::CPCAproj (int *pnParIn, double *pdParIn, double *pdX, double *pdZ, double *pdL, double *pdSDev)
		: m_dwN (pnParIn[0]), m_dwP (pnParIn[1]), m_dwRealN (pnParIn[2]), m_dwK (pnParIn[3])
		, m_nScal (pnParIn[4]), m_nScores (pnParIn[5])
		, m_dZeroTol (pdParIn [0]), m_dCurLambda (0)
		, m_mX (pdX, m_dwN, m_dwP), m_mL (pdL, m_dwP, m_dwK), m_mA (m_dwN, m_dwP)
		, m_vSDev (pdSDev, m_dwK), m_vCurScore (m_dwN)
		, m_vHelpTF (m_dwN)
	{
		if (m_nScores)
			m_mZ.Set (pdZ, m_dwRealN, m_dwK) ;
	}

	void NULL1 (const SMatD &a)	//	takes a pxp loadings matrix, where the last column is not filled, and computes this last column
	{
		ASSERT_TEMPRANGE (0, 1) ;
		ASSERT (a.nrow () == a.ncol ()) ;

		t_size p = a.nrow () ;

		const SVecD &vLastCol = a.GetColRef (p-1) ;
		const SMatD mLpm1 (a.GetColRef (0, p - 1)) ;

		vLastCol.Reset (1) ;
		EO<UOP::Aa_As_sqrB>::VMc (*vLastCol, mLpm1) ;
		EO<SOP::a_sqrt>::V (*vLastCol) ;

		int n = which_max1 (vLastCol) ;

		SVecD vTemp1 (tempRef (0), p - 1), vTemp2 (tempRef (1), p) ;

		CopyRow (*vTemp1, mLpm1, n) ;				//	2do: implement SnVec -> no copy..
		vTemp2.Reset (0) ;
		EO<SOP::ApaBmC>::VMcVct_NC (*vTemp2, mLpm1, vTemp1) ;

 		vTemp2(n) *= -1.0 ;																	//	for not triggering element n in the next line, which is positive by default...
		EO<UOP::if_B_gr_0_AamA>::VVc (*vLastCol, vTemp2) ;									//	assigns the negative signs of vector vTemp2 to vector vCurEVec (assuming that vCurEVec has only positive elemts)
	}

	void CPCAproj::SetSingular (t_size dwK)
	{
		m_mZ.GetColRef (dwK, m_dwK).Reset (0) ;
		m_vSDev.GetDataRef (dwK, m_dwK).Reset (0) ;
		if (!dwK)
			SetDiag (!m_mL) ;
		else
			m_vSDev.GetDataRef (dwK, m_dwK).Reset (-1) ;	//	sets the sdev to -1, indicating, that the according columns of the Loadings Matrix are invalid!
			
	}

	void CPCAproj::Calc ()
	{
		SVecD vPcol (m_dwN), vVH (m_dwP), vHlp (m_dwN), vHlpS (vHlp) ;
		SVecD vCurA (tempRef (0), m_dwP) ;

		SVecD vCurScoreS (*m_vCurScore, m_dwRealN) ;

		//double *pdCurLambda = m_vSDev ;

		t_size i, j ;
		for (i = 0; i < m_dwK; i++)
		{
			const SVecD &vCurEVec = m_mL.GetColRef (i) ;

			vHlp.Reset (0) ;
			EO<SOP::Apa_sqr_B>::VMc (*vHlp, m_mX) ;												//R	vHlp <- rowSums (mY^2)
			m_dwShortN = 0 ;
			EO<UOP::aB_cA_C_le_D>::SVScVc (m_dwShortN, *m_vHelpTF, m_dZeroTol, vHlp) ;			//R	m_vHelpTF <- (vHlp>dZeroTol); m_dwShortN <- sum ()

			if (!m_dwShortN)	//	all observations seem to be concentrated in one point (in the center) when considering the current subspace
			{
				SetSingular (i) ;
				return ;
			}

			vHlpS.Reshape (m_dwShortN) ;
			m_mA.Reshape (m_dwShortN, m_dwP) ;

			EO<SOP::a_sqrt>::V (*vHlp) ;														//R	vHlp <- sqrt (vHlp)
			EO<SOP::divide>::MsMcVcVbc (!m_mA, m_mX, vHlp, m_vHelpTF) ;							//R	m_mA <- (m_mX / vHlp)[m_vHelpTF,]

			m_vCurScore.Reshape (m_dwRealN) ;

			if (i < m_dwP - 1)
			{
				t_size dwBestj = NAI;

				//for (j = 0; j < m_dwShortN; ++j)
				for (j = m_dwShortN - 1; j != NAI; --j)
				{
					CopyRow (*vCurA, m_mA, j) ;														//R	vCurA <- m_mA[j, ]
					vec_mult_mat_t_partial (m_vCurScore, vCurA, m_mX, m_dwRealN, m_dwP, m_dwN) ;	//R	m_vCurScore <- vCurA %*% m_mX[1:m_dwRealN,]	//RR	Y =A %*% t(y[1:n,]);		
					double dScat = ApplyMethod (m_vCurScore, m_nScal) ;								//R dScat = fscale (m_vCurScore)				//RR	pcol = apply (Y, 1, fs)
					if (dwBestj == NAI || m_dCurLambda < dScat)
					{
						dwBestj = j ;																												//RR	istar <- which.max (pcol)
						m_dCurLambda = dScat ;																										//RR	lambda[i] <- pcol[istar]
					}
				}

				CopyRow (*vCurEVec, m_mA, dwBestj) ;												//R	vCurEVec <- m_mA[dwBestj,]					//RR	vhelp <- A[istar,]
				m_vCurScore.Reshape (m_dwN) ;
				m_vCurScore.Reset (0) ;
				EO<SOP::ApaBmC>::VMcVct (*m_vCurScore, m_mX, vCurEVec) ;							//R m_vCurScore <- m_mX %*% vCurEVec			//RR	scorevec <- y%*%(A[istar,])
				Update (vCurEVec) ;

				if (m_nScores)
					Copy (*m_mZ.GetColRef (i), vCurScoreS) ;

				if (i < m_dwK - 1)
					EO<SOP::AsaBmC>::MVcVct (!m_mX, m_vCurScore, vCurEVec) ;						//R	m_mX <- m_mX - m_vCurScore %*% vCurEVec

				m_vSDev (i) = m_dCurLambda ;
//				*pdCurLambda = m_dCurLambda ;
//				++pdCurLambda ;
			}
			else
			{
				NULL1 (m_mL) ;	//	computes the last eigenvector (loadings vector) in m_mL

				m_vCurScore.Reshape (m_dwN) ;
				m_vCurScore.Reset (0) ;
				EO<SOP::ApaBmC>::VMcVct (*m_vCurScore, m_mX, vCurEVec) ;							//R m_vCurScore <- m_mX %*% vCurEVec			//RR	scorevec <- y%*%(A[istar,])
								//	2do -> pass this to BLAS

				//*pdCurLambda = 
				m_vSDev (i) = ApplyMethod (m_vCurScore, m_nScal) ;									//R	dCurLambda = fscale (m_vCurScore)

				if (m_nScores)
					Copy (*m_mZ.GetColRef (i), vCurScoreS) ;

			}
		}
	}

	CPCAprojU::CPCAprojU (int *pnParIn, double *pdParIn, double *pdX, double *pdZ, double *pdL, double *pdSDev)
		: CPCAproj (pnParIn, pdParIn, pdX, pdZ, pdL, pdSDev)
		, m_dwMaxIt (pnParIn[6]), m_dwMaxHalf (pnParIn[7]) {}

	void CPCAprojU::Update (const SVecD &vCurEVec)
	{
		ASSERT_TEMPRANGE (11, 12) ;

		t_size m, kk ;

		SVecN	vScoreSign (tempRef (0), m_dwShortN) ;
		SVecD	vVH (tempRef (11), m_dwP),
				vTScores (tempRef (12), m_dwN) ;

		for (m = m_dwMaxIt; m ; --m)
		{
			EO<SOP::sign>::VsVcVbc (*vScoreSign, m_vCurScore, m_vHelpTF) ;						//R	vScoreSign <- sign (vScoreSign [m_vHelpTF])

			vVH.Reset (0) ;																		//R	{
			EO<UOP::ApaBm_signC>::VtMcVc (*vVH, m_mA, vScoreSign) ;								//R		vVH <- t (m_mA) %*% vScoreSignS
			double dNewObj, dSqSum = 0 ;														//R		dSqSum <- sum (vVH^2)
			EO<SOP::Apa_sqr_B>::SVc (dSqSum, vVH) ;												//R	}

			//	...
			//		this part of the original R source has been merged with the following for loop
			//	...

			for (kk = 0; kk <= m_dwMaxHalf; ++kk)
			{
				if (kk)
				{
					dSqSum = 0 ;																//R	vVH <- (vVH + vCurEVec) / 2
					EO<UOP::Ba_BpC_d2_Apa_sq_B>::SVVc (dSqSum, *vVH, vCurEVec) ;				//R	dSqSum <- sum (vVH^2)
				}

				EO<SOP::a_divide>::VSc (*vVH, sqrt (dSqSum)) ;									//R	vVH <- vVH / sqrt (dSqSum)

				vTScores.Reset (0) ;															//R
				EO<SOP::ApaBmC>::VMcVct (*vTScores, m_mX, vVH) ;								//R	VTScores <- mX %*% vVH
				dNewObj = ApplyMethod (vTScores, m_nScal) ;										//R	dNewObj <- fscale (vTScores)

				if (dNewObj >= m_dCurLambda)
					break ;
			}

			if (dNewObj < m_dCurLambda)
				break ;

			Copy (*m_vCurScore, vTScores) ;														//R	m_vCurScore <- vTScores
			Copy (*vCurEVec, vVH) ;																//R	vCurEVec <- vVH
			m_dCurLambda = dNewObj;
		}
	}


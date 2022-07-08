#include "L1Median.h"

	double calObj (const SCMatD &mXc, const SCVecD &vMed)
	{
		ASSERT_TEMPRANGE (0, 0) ;
		SVecD vTemp (tempRef (0), mXc.nrow ()) ;

		vTemp.Reset (0) ;
		EO<SOP::Apa_sqr_BsC>::VMcVct_NC (*vTemp, mXc, vMed) ;

		double dSum = 0 ;
		EO<SOP::Apa_sqrt_B>::SVc (dSum, vTemp) ;
		return dSum ;
	}

	int l1median_HoCr (const SCMatD &mX, const SVecD &vdMedian, double dZeroTol, double dTol, int dwMaxit, int nTrace, int *pdwIterCount)
	{
		const double dLog2 = log ((double) 2.0) ;
		const int n = mX.nrow (), &p = mX.ncol () ;


		SVecD vdMedianOld (p), vdNorms (n), vdNormsOrdered, vdWeights (n), vdDelta (p) ;
		SMatD mXc (n, p) ;
		SVecN vinter (n) ;

		int nK, &k = pdwIterCount ? *pdwIterCount : nK, dwNStep, dwMaxHalf = 0 ;

		double * pdNorms, * const pdStartNorms = vdNorms, * const pdEndNorms = vdNorms.GetDataEnd (), * const pdNormsHalf = pdStartNorms + ((n + 1) >> 1) ;
		double * const pdWeights = vdWeights ;
		int *pnOrder, * const pnStartOrder = vinter ;

		double dND, dObj, dObjOld = 0 ;


		for (k = 0; k < dwMaxit; k++)
		{
			vdNorms.Reset (0) ;
			EO<CL1Median_VZ::AaCmD_BpaAmA>::MVMcVct (!mXc, *vdNorms, mX, vdMedian) ;

			if (!k)	//	in the first round we (cheaply) calculate the "old" (=initial) objective function's value, as half of it we've already calculated with vector "vdNorms"
				EO<SOP::Apa_sqrt_B>::SVc (dObjOld, vdNorms) ;

			sort_order (*vdNorms, *vinter) ;

			double dSumWeights = 0 ;

			pnOrder = pnStartOrder ;
			pdNorms = pdStartNorms ;


			while (pdNorms < pdEndNorms)
			{
//				if (*pdNorms <= dZeroTol)
				if (*pdNorms <= 0)
				{
					pdWeights[*pnOrder] = 0 ;

					if (pdNorms > pdNormsHalf)	//	there's more than half of the values concentrated at the current median estimation -> return this estimation
						return 3 ;
				}
				else
					dSumWeights += (pdWeights[*pnOrder] = pow (*pdNorms, -0.5)) ;
				++pdNorms ;
				++pnOrder ;
			}

			vdDelta.Reset (0) ;
			EO<SOP::ApaBmC>::VtMcVc (*vdDelta, mXc, vdWeights) ;

			dND = 0 ;
			EO<SOP::BdaC_Apa_sqr_B>::SVSc(dND, *vdDelta, dSumWeights) ;
			dND = sqrt (dND) ;

			if (nTrace >= 3)	meal_printf ("nd at %g in iteration %d (tol at %g)\r\n", dND, k, dTol) ;

			vdMedianOld.Copy_NC (vdMedian) ;
			EO<SOP::a_plus>::VVc_NC (*vdMedian, vdDelta) ;
			dObj = calObj (mX, vdMedian) ;

			if (dND < dTol)	//	converged
				return 0 ;

			dwMaxHalf = (int) ceil (log (dND / dTol) / dLog2) ;

			for (dwNStep = dwMaxHalf; dObj >= dObjOld; --dwNStep)	//	do step-halving
			{
				EO<SOP::a_divide>::VSc (*vdDelta, 2) ;
				EO<SOP::add>::VVcVc_NC (*vdMedian, vdMedianOld, vdDelta) ;
				dObj = calObj (mX, vdMedian) ;

				if (!dwNStep)
				{	
					vdMedian.Copy (vdMedianOld) ;
					dObj = dObjOld ;
					return 2 ;		//	step-halving failed
				}
			}

			dObjOld = dObj ;
		}

		return 1 ;	//	did not converge
	}

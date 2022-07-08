#ifndef INCLUDE_PCAPP_PCAGRID_H
#define INCLUDE_PCAPP_PCAGRID_H

#include "pcaPP.h"

	double ngpf (const double &d) ;

	class CPCAGrid
	{
	public:
		CPCAGrid (int *pnParamIn, int *pnParamOut, double *pdParamIn, double *pdData, double *pdLoadings, double *pdSDev, double *pdObj/*, double *pdMaxMaha*/) ;

		int Calc () ;

	protected:
//		void GridPlane (double &dNL, double &dNCL, double &dScat, double &dObj, double dCurSplit) ;
		void GridPlane (double dCurSplit) ;
		void EvalDirection (const double dCos, const double dSin) ;
		virtual double CalcObj (const double dCos, const double dSin, double &dScat, double &dScatOrth) ;
		double CalcProjScat (const double dCos, const double dSin) ;
		double CalcScatTrimmed (double dCos, double dSin, double dScat, double dScatOrth) ;
		double CalcVarTrimmed (double dCos, double dSin, double dScat, double dScatOrth) ;
		virtual void OnCalcPC () {}
		virtual void InitPenalty () {} 

//		void CalcMaha (const double dScat) ;

		double ApplyMethod (const SCVecD &v) ;
		void ApplyMethod (const SCMatD &m, SVecD &v) ;
		double ApplyMethodMean (const SCMatD &m) ;

		inline SMatD &TempY () { return m_mY[m_dwTempYIdx] ; }
		inline SMatD &TempYC () { return m_mY[1-m_dwTempYIdx] ; }
		void SwapTempY () { m_dwTempYIdx = 1 - m_dwTempYIdx ;}

		void BackTransform () ;

		void RemoveLoading (/*DWORD i*/) ;
		void AddLoading (const double &dNL, const double &dNCL) ;

		const t_size m_dwN, m_dwP, m_dwK, m_dwSplitCircle, m_dwMaxIter, m_dwMethod, m_dwTrace, m_dwkIni, m_dwCheckOrth ;
		int &m_nReturn ;
		const double m_dZeroTol ;

		SMatD m_mX, m_mL, m_mY[2] ; //, m_mTempPP, m_mTempPN ;
		SVecD m_vAfin, m_vAfinBest, m_vScl, m_vYOpt, m_vSDev, m_vObj ; //, m_vTempN, m_vTempN2, m_vTempN3;
		SVecD m_vCurY, m_vProj ;//, m_vMaxMaha ;
		SVecN m_vOrd ;

		t_size m_dwCurK, m_dwCurP ;// iteration variables
		t_size m_dwPSub, m_dwTempYIdx ;

		///

		double m_dBestObj, m_dCurScat, m_dCurScatOrth, m_dNL, m_dNCL ;
		double * const m_pdProj, * const m_pdCurLC, * const m_pdCurLCEnd, *m_pdCurY ;	//	, * const m_pdEndProj
	} ;

	class CsPCAGrid : public CPCAGrid
	{
	public:
		CsPCAGrid (int *pnParamIn, int *pnParamOut, double *pdParamIn, double *pdData, double *pdLoadings, double *pdSDev, double *pdObj/*, double *pdMaxMaha*/, double *pdLambda, double *pdBackTransHD) ;

	protected:
		virtual void OnCalcPC () ;
		double GetPenalty (const double& dCos, const double& dSin) ;
		virtual double CalcObj (const double dCos, const double dSin, double &dScat, double &dScatOrth) ;
		virtual void InitPenalty () ;

		const int m_nGloScatter, m_nSpeedUp ;
		const t_size m_dwPHD ;
		const double m_dQ, m_dS ;
		const BOOL m_bUseQ, m_bUseS ;

		SMatD m_mBackTransHD, m_mBackProj ;
		SVecD m_vLambda, m_vLoadHD, m_vTempP, m_vTempPSub, m_vSumLoadOthers, m_vSumLoadThis ;

		double m_dGloScatter, m_dCurLambda ;
		double m_dLoadSumThis, m_dLoadSumOther ;
	} ;

#endif  //  #ifndef INCLUDE_PCAPP_PCAGRID_H

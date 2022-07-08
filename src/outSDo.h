#include "pcaPP.h"

class CSDoOut
{
	public:
		CSDoOut (int *pnParIn, double *pdX, double *pdMaxMaha, int *pnNChanged) ;

		void Calc () ;

	protected:
		void CalcCenter () ;

		void IterObs () ;
		void IterDiffObs (int n) ;
		void IterRand (int n) ;
		void IterRandDiffObs (int n) ;

		int DoDir (const SCVecD &vLoad) ;

		const t_size m_dwN, m_dwP, m_dwIterMethod, m_dwIterParam, m_dwCenterMethod, m_dwScatterMethod, m_dwReset ;

		SMatD m_mX ;
		SVecD m_vMaxMaha ;

		const t_size m_dwNDir ;
		SVecD m_vXProj, m_vCurDir ;
//		SVecN m_vChanged ;
		int * const m_pnNChanged ;
		double * const m_pdDiff ;


		double *m_pdXProj, * const m_pdEndXProj, *m_pdMaxMaha ;
} ;

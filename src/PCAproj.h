#include "pcaPP.h"


	class CPCAproj
	{
	public:
		CPCAproj (int *pnParIn, double *pdParIn, double *pdX, double *pdZ, double *pdL, double *pdSDev) ;

		void Calc () ;


	protected:
		void SetSingular (t_size dwK) ;

		virtual void Update (const SVecD &vCurEVec) {}

		const t_size m_dwN, m_dwP, m_dwRealN, m_dwK ;
		t_size m_dwShortN ;
		const int m_nScal, m_nScores ;

		const double m_dZeroTol ;
		double m_dCurLambda ;

		SMatD m_mX, m_mL, m_mZ, m_mA ;
		SVecD m_vSDev, m_vCurScore ;
		SVecN m_vHelpTF ;

	} ;

	class CPCAprojU : public CPCAproj
	{
	public:
		CPCAprojU (int *pnParIn, double *pdParIn, double *pdX, double *pdZ, double *pdL, double *pdSDev) ;

	protected:

		const t_size m_dwMaxIt, m_dwMaxHalf ;

		virtual void Update (const SVecD &vCurEVec) ;

	} ;


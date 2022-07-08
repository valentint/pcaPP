#ifndef ES_PCAPP_H
#define ES_PCAPP_H

#ifdef ES_DEV_ENV
	#include "../../../SMat/smat.h"
#else
	#include "smat.h"
#endif
	
	class UOP																		//	User Operators
	{
	public:
		class Aa_AsDmB_dC			{ CALC_4_1(void) { a = (a-d*b)/c ; }	} ;		//	adding and removing loadings for pcaPP::sPCAGrid...
		class Aa_AmC_p_DmB			{ CALC_4_1(void) { a = a * c + d * b ; }	} ;
		class Apa_abs_BmDpCmE_		{ CALC_5_1(void) { a += fabs (b * d + c * e) ; }	} ;

		class Apa_pow_abs_B0mCpb1mD_B2	{ CALC_4_1(void) { a += pow (fabs (b[0] * c + b[1] * d), b[2]) ; }	} ;


		class Apa_abs_B				{ CALC_2_1(void) { a += fabs (b) ; }	} ;

		class Apa_pow_abs_C_B		{ CALC_3_1(void) { a += pow (fabs (c), b) ; }	} ;

		class aB_cA_C_le_D			{ CALC_4_2(void) { if (c < d) {b = 1 ; a += 1 ;} else b = 0 ; }	} ;	//	b = c < d; a += c < d
		class Ba_BpC_d2_Apa_sq_B	{ CALC_3_2(void) { b = (b + c) / 2.0; a += sm_sqr (b) ; } } ;
		class ApaBm_signC			{ CALC_3_1(void) { if (c < 0) a-= b ; else a += b ; }	} ;

		class Aa_As_sqrB			{ CALC_2_1(void) { a -= sm_sqr (b) ; } } ;
		class if_B_gr_0_AamA		{ CALC_2_1(void) { if (b > 0) a = -a ; } } ;

//		class Apa_sq_B_B_setsign_C	{ CALC_3_2(void) { a += sm_sqr (b) ; if (c < 0) b = -b ; } } ;		//	2do: delete
//		class Aa_AmB_pC				{ CALC_3_1(void) { a = (a * b) + c ; }	} ;
//		class Aa_AsC_dB				{ CALC_3_1(void) { a = (a - c) / b ; }	} ;
	} ;

	double ApplyMethod (const SCVecD &v, const int nMethod) ;
	double ApplyMethod_V (const SVVecD &v, const int nMethod) ;

#endif	//	#ifndef ES_PCAPP_H


#ifdef ES_DEV_ENV
	#include "../../../SMat/smat.h"
#else
	#include "smat.h"
#endif
#include "defint64.h"

	void qn		(double &dQn, const	double *pX, const int n) ;
	void qn_V	(double &dQn,		double *pX, const int n) ;

	void qn_nc (double &dQn, const double *pX, const int n) ;
	double qn_raw (double *pY, const int n) ;

	double qn (const SVDataD &a) ;
	double qn (const SCDataD &a) ;

	double qn_corrN (const int n, const double dQnCNorm = 2.21914446598508) ;

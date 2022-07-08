#include "R_package.h"
#include "math.h"

void Hess_Sub (int p, double *pdX_i, double *pdMu, double *pdHess, double *pdTempP)
{
	int l, k ;

	double dNorm = 0 ;
	for (l = p - 1; l != -1; l--)
	{
		double &dCur = pdTempP [l] = pdX_i[l] - pdMu[l] ;
		dNorm += dCur * dCur ;
	}

	dNorm = 1 / sqrt (dNorm) ;
	double dNorm3 = pow (dNorm, 3.0) ;

	for (l = p - 1; l != -1; l--)
	{
		pdHess [l * p + l] += dNorm ;

		//for (k = p - 1; k != -1; k--)
		for (k = l; k != -1; k--)
			pdHess [l * p + k] -= pdTempP [l] * pdTempP[k] * dNorm3 ;
	}

}

void Hess (int p, int n, double *pdX, double *pdMu, double *pdHess, double *pdTempP1, double *pdTempP2)
{
	int i, j ;

	for (i = p - 1; i != -1; i--)
		for (j = p - 1; j != -1; j--)
			pdHess [i + j * p] = 0 ;

	for (i = n - 1; i != -1; i--)
	{
		for (j = p - 1; j != -1; j--)
			pdTempP2 [j] = pdX[i + j * n] ;

		Hess_Sub (p, pdTempP2, pdMu, pdHess, pdTempP1) ;
	}

	for (i = p - 1; i != -1; i--)
		for (j = i - 1; j != -1; j--)
			pdHess [i + j * p] = pdHess [i * p + j] ;
}

void Hess_Sub_R (int *pnPar, double *pdX_i, double *pdMu, double *pdHess)
{
	const int &p = pnPar[0]  ;
	double *pdTempP = new double [p] ;

	Hess_Sub (pnPar [0], pdX_i, pdMu, pdHess, pdTempP) ;
// VT::04.12.2017
// Causes an warning in clang++: "'delete' applied to a pointer that was allocated with 'new[]'; did you mean 'delete[]'? [-Wmismatched-new-delete]"
// 	- replace delete by delete[]
//	delete pdTempP ;
	delete[] pdTempP ;

}

void Hess_R (int *pnPar, double *pdX, double *pdMu, double *pdHess)
{
	double *pdTempP1 = new double [pnPar[0]], *pdTempP2 = new double [pnPar[0]] ;
	Hess (pnPar[0], pnPar[1], pdX, pdMu, pdHess, pdTempP1, pdTempP2) ;

// VT::04.12.2017
// -  as above
	delete[] pdTempP1 ;
	delete[] pdTempP2 ;
}

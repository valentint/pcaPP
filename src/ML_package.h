#ifdef MATLAB_MEX_FILE



#ifdef ES_DEV_ENV
//	#include "../../../RDev/R.Inc.h"
	#include "../../../SMat/smat.def.h"
	#include "../../../SMat/smat.h"
	#include "../../../MLDev/ML_meal.h"
#else
//	#include "R.Inc.h"
	#include "smat.def.h"
	#include "smat.h"
	#include "ML_meal.h"
#endif





#define MEX_FUNCTION_EXPORTS

#ifdef MEX_FUNCTION_EXPORTS
#define MEX_FUNCTION_API __declspec(dllexport)
#else
#define MEX_FUNCTION_API __declspec(dllimport)
MEX_FUNCTION_API void mexFunction(int nlhs, mxArray* plhs[], int nrhs, mxArray* prhs[]);
#endif



	void ML_l1median_HoCr (int nlhs, mxArray* plhs[], int nrhs, mxArray *prhs[]) ;
	void ML_l1median_VaZh (int nlhs, mxArray* plhs[], int nrhs, mxArray *prhs[]) ;
	void ML_qn (int nlhs, mxArray* plhs[], int nrhs, mxArray *prhs[]) ;
	void ML_PCAgrid (int nlhs, mxArray* plhs[], int nrhs, mxArray *prhs[]) ;
	void ML_sPCAgrid (int nlhs, mxArray* plhs[], int nrhs, mxArray *prhs[]) ;
	void ML_PCAprojU (int nlhs, mxArray* plhs[], int nrhs, mxArray *prhs[]) ;
	void ML_PCAproj (int nlhs, mxArray* plhs[], int nrhs, mxArray *prhs[]) ;



#endif

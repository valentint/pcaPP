#include "qnn.h"


	double pull (double const * const pA, const int n, int k)
	{
		ASSERT_TEMPRANGE (0, 0) ;
		double * pB = tempRef (0, pB, n) ;
		Copy (pB, pA, n) ;

		double ax, buffer ;
		int l = 0, lr = n - 1, jnc, j ;

		while (l < lr)					//	20	//	2do: put into pull fct
		{
			ax = pB[k] ;

			jnc = l ;
			j = lr ;

			while (jnc <= j)			//	30
			{
				while (pB[jnc] < ax)	//	40
					++jnc ;

				while (pB[j] > ax)		//	50
					--j ;

				if(jnc <= j)			//	60
				{
					sm_swap (pB [jnc], pB[j], buffer) ;
					++jnc ;
					--j ;
				}
			}							//	70

			if (j < k)
				l = jnc ;
			if (k < jnc)
				lr = j ;
		}

		return pB [k] ;
	}

	double whimed (double * const pA, int * const pIw, int n)
	{
		ASSERT_TEMPRANGE (1, 2) ;
		int i ; //= (n + 1) >> 1 ;

		double * pAcand = tempRef (2, pAcand, n) ;
		int * pIWcand = tempRef (1, pIWcand, n) ;

		int64_t nWtotal ; 
		sum (pIw, n, nWtotal) ;
		if (!nWtotal)
			return meal_NaN () ;

		int64_t nTemp, nWrest = 0,nWleft= 0, nWmid = 0,nWright = 0 ;
		int nKcand, nn = n ;

		while (1)
		{
			double dTrial = pull(pA, nn, (nn >> 1)) ;

			nWleft = nWright = nWmid = 0 ;
			for (i = 0; i < nn; ++i)
			{
				if (pA[i] < dTrial)
					nWleft += pIw [i] ;
				else if (pA [i] > dTrial)
					nWright += pIw [i] ;
				else
					nWmid += pIw [i] ;
			}

			nTemp = nWrest + nWleft ;
			if ((nTemp << 1)> nWtotal)
			{
				nKcand = 0 ;
				for (i = 0; i < nn; ++i)
				{
				  if (pA [i] < dTrial)
				  {
					  pAcand [nKcand] = pA [i] ;
					  pIWcand [nKcand] = pIw [i] ;
					  ++nKcand ;
				  }
				}

				nn = nKcand ;
			}
			else if (((nTemp + nWmid) << 1) > nWtotal)
				return dTrial ;
			else
			{
				nKcand = 0 ;
				for (i = 0; i < nn; ++i)
				{
					if(pA [i] > dTrial)
					{
						pAcand [nKcand] = pA [i] ;
						pIWcand [nKcand] = pIw [i] ;
						++nKcand ;
					}
				}
				nn = nKcand ;
				nWrest += nWleft + nWmid ;
			}
			Copy (pA, pAcand, nn) ;
			Copy (pIw, pIWcand, nn) ;
		}
	}

#ifdef _MSC_VER
	#define NO_INLINE								//	MS compilers don't make problems here when using the helper functions..
#else
	#define NO_INLINE	 __attribute__ ((noinline))	//	other compilers (e.g. MinGW) are not supposed to inline these functions...
#endif

													//	workaround of some compiler optimization - issue:
													//	when computing  a-b < x, whereas a-b == x the resulting value was 
													//	sometimes "true" on windows machines (using MS VS6.0 and Mingw, various versions)
													//	in some very rare occasions this caused the qn algo to end in an infinite loop.
	BOOL NO_INLINE  isgr_s (const double &a, const double &b) 	//
	{												//	the same code worked without problems on linux.
		return a > b ;								//	an assumption is, that some compilers "optimize" the expression 
	}												//		
													//						"a-b < x"
	BOOL NO_INLINE isle_s (const double &a, const double &b)	//				to
	{												//			"a - b - x < 0" or "a - x - b < 0"  (or sth similar)
		return a < b ;								//
	}												//	which then gives numerical problems.
													//	
													//	QUICKFIX:
													//	Thus by using functions "isgr_s" and "isle_s" which are not allowed to be 
													//	inlined, this optimization is avoided and the algorithm runs smoothly. \o/
													//
													//	SOLUTION:
													//	Directly turning off the corresponding optimization for these lines.
													//		("#pragma optimize ("", off")" didn't help so far)

	double qn_raw (double *pY, const int n)
	{
		TEMP_GUARD ;
		ASSERT_TEMPRANGE (3, 8) ;
		const int ns1 = n - 1 ;

		double * const pWork = tempArray<double> (8, n) ;
		int * const pLeft = tempArray<int> (7, n), * const pRight = tempArray<int> (6, n), * const pWeight = tempArray<int> (5, n), * const pQ = tempArray<int> (4, n), * const pP = tempArray<int> (3, n) ;

		tempArray<double> (0, n) ;
		tempArray<double> (1, n) ;
		tempArray<double> (2, n) ;
		int i ;
		double dTrial ;

		const int64_t h = n/2+1 ;
		int64_t 
			k = h * (h-1) >> 1 ,
			jhelp = (n*((int64_t)n+1)) >> 1,
			knew = k + jhelp,
			nL = jhelp,
			nR = ((int64_t) n) * n,
			dwSumQ, dwSumP, j ;

		meal_sort (pY, n) ;

		for (i = n - 1; i != -1; --i)
		{
			pLeft[i] = n - i ;
			pRight [i] = n ;
		}

		while (nR - nL > n)												 //F 200   continue	
		{
			j = 0 ;
			for (i = 1; i < n; ++i)
			{
				if (pLeft[i] < pRight [i])
				{
					pWeight [j] = pRight[i] - pLeft[i] ;				//F	weight(j)=right(i)-left(i)+1
					jhelp = pLeft[i] + (pWeight[j] >> 1) ;				//F	jhelp=left(i)+weight(j)/2
					pWork[j] = pY[i] - pY[n - jhelp - 1] ;				//F	work(j)=y(i)-y(n+1-jhelp)
					++j ;
				}
			}
			dTrial = whimed (pWork, pWeight, int (j)) ;					//F	trial=whimed(work,weight,j-1)

			dwSumP = dwSumQ = j = 0 ;									//F j = n + 1

			for (i = n - 1; i != -1; --i)
			{

				while (j < n && isle_s ((pY[i] - pY[ns1 - j]), dTrial))	//F	if ((j.lt.n).and.((y(i)-y(n-j)).lt.trial)) then
					++j ;
//				while (j < n && (pY[i] - pY[ns1 - j]) < dTrial)			//F	if ((j.lt.n).and.((y(i)-y(n-j)).lt.trial)) then
//					++j ;

				pP[i] = int (j) ;
				dwSumP += int (j) ;										//F	sumP+P(i)
			}
			j = n ;
			for (i = 0; i < n; i++)
			{
				while (isgr_s(pY[i] - pY[n - j], dTrial))				//F	if ((y(i)-y(n-j+2)).gt.trial) then
					--j ;
//				while (pY[i] - pY[n - j] > dTrial)						//F	if ((y(i)-y(n-j+2)).gt.trial) then
//					--j ;
				pQ[i] = int (j) ;
				dwSumQ += int (j) ;
			}

			if (knew <= dwSumP)
			{
				Copy (pRight, pP, n) ;
				nR = dwSumP ;
			}
			else if (knew > dwSumQ)
			{
				Copy (pLeft, pQ, n) ;
				nL = dwSumQ ;
			}
			else
				return dTrial ;
		}
		int jj ;

		j = 0 ; //j=1
		for (i = 1; i < n; ++i)											//F	do 90 i=2,n
		{
			if (pLeft[i] < pRight[i])									//F	if (left(i).le.right(i)) then
			{
				for (jj = pLeft[i]; jj < pRight[i]; ++jj)				//F	do 100 jj=left(i),right(i)
				{
					pWork[j] = pY[i] - pY[ns1-jj] ;						//F	work(j)=y(i)-y(n-jj+1)
					++j ;
				}
			}															//F	100
		}																//F	90

		return pull (pWork, int (j), int (knew-nL - 1)) ;				//F	Qn=pull(work,j-1,knew-nL)
	}

	double qn_corrN (const int n, const double dQnCNorm)
	{
		if (n <= 9)
		{
			static const double adCorrFact [] = {0.400, 0.993, 0.514, 0.845, 0.612, 0.859, 0.670, 0.874} ;
			return dQnCNorm * adCorrFact[n - 2] ;
		}

		if (n & 1)				
			return dQnCNorm * n / (n + 1.4) ;	//	odd
		return dQnCNorm * n / (n + 3.8) ;		//	even
	}

	void qn_nc (double &dQn, const double *pX, const int n)
	{
		TEMP_GUARD ;
		ASSERT_TEMPRANGE (9, 9) ;

		double * const pY = tempArray<double> (9, n) ;
		Copy (pY, pX, n) ;

		dQn = qn_raw (pY, n) ;
	}

	void qn_V (double &dQn, double *pX, const int n)
	{
		dQn = qn_raw (pX, n) ; 
		dQn *= qn_corrN (n) ;
	}

	void qn (double &dQn, const double *pX, const int n)
	{
		qn_nc (dQn, pX, n) ;
		dQn *= qn_corrN (n) ;
	}

/*
	double qn (const double *pX, t_size n)
	{
		double dQn ;
		qn (dQn, pX, n) ;
		return dQn ;
	}

	EXPORT void ex_pull (int *pnParIn, double *pnParOut, double const * const pA)
	{
		pnParOut[0] = pull (pA, pnParIn[0], pnParIn[1]) ;
	}

	EXPORT void ex_whimed (int *pnParIn, double *pnParOut, double * const pA, int * const pIw)
	{
		pnParOut[0] =  whimed (pA, pIw, pnParIn[0]) ;
	}
*/


/*
	void qn (const double *x, int *npLength, double *pdQn)								//	used by sPCApp: 2do: change scale function definition ...
	{
		qn (*pdQn, x, *npLength) ;
	}
*/

	double qn (const SVDataD &a)
	{
		double dRet ;
		qn (dRet, a.GetData (), a.size ()) ;
		return dRet ;
	}

	double qn (const SCDataD &a)
	{
		double dRet; 
		qn (dRet, a.GetData (), a.size ()) ;
		return dRet ;
	}

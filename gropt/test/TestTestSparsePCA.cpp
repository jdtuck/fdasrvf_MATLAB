
#include "TestTestSparsePCA.h"

#if !defined(MATLAB_MEX_FILE) && defined(TESTTESTSPARSEPCA)

std::map<integer *, integer> *CheckMemoryDeleted;

int main(void)
{
	init_genrand(0);
	CheckMemoryDeleted = new std::map<integer *, integer>;

	integer p = 100, r = 5, n = 20;
	double epsilon = 1e-4;
	double mu = 1e-6;

	double *B = new double[p * n + r + 2 * p * r];
	double *Dsq = B + n * p;
	double *X = Dsq + r;
	double *Xopt = X + p * r;
	for (integer i = 0; i < p * n; i++)
		B[i] = genrand_gaussian();

	integer minpn = (p > n) ? n : p;
	integer maxpn = (p > n) ? p : n;
	double *Bcopy = new double[p * n + minpn + p * minpn + minpn * n];
	double *S = Bcopy + p * n;
	double *U = S + minpn;
	double *Vt = U + p * minpn;

	integer *iwork = new integer[8 * minpn];
	integer lwork = -1;
	double workoptsize = 0;
	integer np = n * p, inc = 1, unref = 1, info;
	dcopy_(&np, B, &inc, Bcopy, &inc);
	char *jobs = const_cast<char *> ("s");
	// obtain optimal size for work
	std::cout << "start obtaining optimal size for dgesdd" << std::endl;
	dgesdd_(jobs, &p, &n, Bcopy, &p, S, U, &p, Vt, &minpn, &workoptsize, &lwork, iwork, &info);
	lwork = static_cast<integer> (workoptsize);
	std::cout << "the optimal size is " << lwork << std::endl;
	double *work = new double[lwork];
	dgesdd_(jobs, &p, &n, Bcopy, &p, S, U, &p, Vt, &minpn, work, &lwork, iwork, &info);
	delete[] iwork;
	delete[] work;
	std::cout << "info of dgesdd:" << info << std::endl;
	if (info != 0)
	{ // dgesdd fails, use dgesvd instead
		dcopy_(&np, B, &inc, Bcopy, &inc);
		lwork = -1;
		std::cout << "start obtaining optimal size for dgesvd" << std::endl;
		dgesvd_(jobs, jobs, &p, &n, Bcopy, &p, S, U, &p, Vt, &minpn, &workoptsize, &lwork, &info);
		lwork = static_cast<integer> (workoptsize);
		std::cout << "the optimal size is " << lwork << std::endl;
		work = new double[lwork];
		dgesvd_(jobs, jobs, &p, &n, Bcopy, &p, S, U, &p, Vt, &minpn, work, &lwork, &info);
		std::cout << "info of dgesvd:" << info << std::endl;
		delete[] work;
		if (info != 0)
		{
			std::cout << "dgesvd fails!" << std::endl;
			return 0;
		}
	}
	for (integer i = 0; i < r; i++)
		Dsq[i] = S[i] * S[i];
	for (integer i = 0; i < p * r; i++)
	{
		X[i] = U[i];
	}
	delete[] Bcopy;

	testTestSparsePCA(B, Dsq, p, n, r, epsilon, mu, X, Xopt);
	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			std::cout << "Global address:" << iter->first << ", sharedtimes:" << iter->second << std::endl;
	}
	delete CheckMemoryDeleted;
	delete[] B;

#ifdef _WIN64
#ifdef _DEBUG
	_CrtDumpMemoryLeaks();
#endif
#endif
	return 0;
}
#endif

#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs < 3)
    {
        mexErrMsgTxt("The number of arguments should be at least three.\n");
    }
	double *B, *D, *X, *Xopt;
	B = mxGetPr(prhs[0]);
	D = mxGetPr(prhs[1]);
	X = mxGetPr(prhs[2]);
	double epsilon = mxGetScalar(prhs[3]);
	double mu = mxGetScalar(prhs[4]);
	/* dimensions of input matrices */
	integer p, n, r;
	p = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);
	r = mxGetM(prhs[1]);
    if(mxGetN(prhs[1]) != 1)
    {
        mexErrMsgTxt("The size of matrix Dsq is not correct.\n");
    }
    if(mxGetM(prhs[2]) != p || mxGetN(prhs[2]) != r)
    {
        mexErrMsgTxt("The size of matrix X is not correct.\n");
    }
    
	std::cout << "(p, n, r):" << p << "," << n << "," << r << std::endl;

	/*create output matrix*/
	plhs[0] = mxCreateDoubleMatrix(p, r, mxREAL);
	Xopt = mxGetPr(plhs[0]);

	init_genrand(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;
	testTestSparsePCA(B, D, p, n, r, epsilon, mu, X, Xopt);
	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			std::cout << "Global address:" << iter->first << ", sharedtimes:" << iter->second << std::endl;
	}
	delete CheckMemoryDeleted;
	return;
}

#endif

void testTestSparsePCA(double *B, double *Dsq, integer p, integer n, integer r, double epsilon, double mu, double *X, double *Xopt)
{
	// reduce n to r
	ObliqueVariable ObliqueX(p, r);

	Oblique Mani(p, r);
	Mani.CheckParams();
	//Mani.CheckIntrExtr(&ObliqueX);
	//Mani.CheckRetraction(&ObliqueX);
	//Mani.CheckDiffRetraction(&ObliqueX);
	//Mani.CheckLockingCondition(&ObliqueX);
	//Mani.CheckcoTangentVector(&ObliqueX);
	//Mani.CheckIsometryofVectorTransport(&ObliqueX);
	//Mani.CheckIsometryofInvVectorTransport(&ObliqueX);
	//Mani.CheckVecTranComposeInverseVecTran(&ObliqueX);
	//Mani.CheckTranHInvTran(&ObliqueX);
	//Mani.CheckHaddScaledRank1OPE(&ObliqueX);

	double *ObliqueXptr = ObliqueX.ObtainWriteEntireData();
	for (integer i = 0; i < ObliqueX.Getlength(); i++)
		ObliqueXptr[i] = X[i];
	ObliqueTestSparsePCA Prob(B, Dsq, epsilon, mu, p, n, r);
	Prob.SetDomain(&Mani);

	//Prob.CheckGradHessian(&ObliqueX);

	// test LRBFGS
	std::cout << "********************************Check all line search algorithm in LRBFGS*************************************" << std::endl;
	for (integer i = 0; i < 1; i++)
	{
		LRBFGS *LRBFGSsolver = new LRBFGS(&Prob, &ObliqueX);
		LRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		LRBFGSsolver->DEBUG = ITERRESULT;
		LRBFGSsolver->Max_Iteration = 10000;
		LRBFGSsolver->OutputGap = 1000;
// 		LRBFGSsolver->Stop_Criterion = FUN_REL;
// 		LRBFGSsolver->Tolerance = 1e-4;
		LRBFGSsolver->CheckParams();
		LRBFGSsolver->Run();
		const Element *xopt = LRBFGSsolver->GetXopt();
		std::cout << "[0:0.0005)  :" << GetNumBetweenC1andC2(xopt, 0, 0.0005) << std::endl;
		std::cout << "[0.0005:0.1):" << GetNumBetweenC1andC2(xopt, 0.0005, 0.1) << std::endl;
		std::cout << "[0.1:1)    :" << GetNumBetweenC1andC2(xopt, 0.1, 1) << std::endl;
		if (Xopt != nullptr)
		{
			const double *xoptptr = xopt->ObtainReadData();
			for (integer i = 0; i < xopt->Getlength(); i++)
				Xopt[i] = xoptptr[i];
		}
		delete LRBFGSsolver;
	}

	// 	// test RTRNewton
	// 	std::cout << "********************************Check RTRNewton*************************************" << std::endl;
	// 	RTRNewton RTRNewtonsolver(&Prob, &ObliqueX);
	// 	RTRNewtonsolver.DEBUG = ITERRESULT;
	// 	RTRNewtonsolver.OutputGap = 5;
	// 	RTRNewtonsolver.Max_Iteration = 10000;
	// 	RTRNewtonsolver.CheckParams();
	// 	RTRNewtonsolver.Run();
	// 	const Element *xopt = RTRNewtonsolver.GetXopt();
	// 	std::cout << "[0:0.0005)  :" << GetNumBetweenC1andC2(xopt, 0, 0.0005) << std::endl;
	// 	std::cout << "[0.0005:0.1):" << GetNumBetweenC1andC2(xopt, 0.0005, 0.1) << std::endl;
	// 	std::cout << "[0.1:1)    :" << GetNumBetweenC1andC2(xopt, 0.1, 1) << std::endl;

	//// test RCG
	//std::cout << "********************************Check all Formulas in RCG*************************************" << std::endl;
	//for (integer i = 0; i < RCGMETHODSLENGTH; i++)
	//{
	//	Mani.HasHHR = true;
	//	RCG *RCGsolver = new RCG(&Prob, &ObliqueX);
	//	RCGsolver->RCGmethod = static_cast<RCGmethods> (i);
	//	RCGsolver->LineSearch_LS = STRONGWOLFE;
	//	RCGsolver->LS_beta = 0.1;
	//	RCGsolver->DEBUG = ITERRESULT;
	//	RCGsolver->CheckParams();
	//	RCGsolver->Run();
	//	const Element *xopt = RCGsolver->GetXopt();
	//std::cout << "[0:0.0005)  :" << GetNumBetweenC1andC2(xopt, 0, 0.0005) << std::endl;
	//std::cout << "[0.0005:0.1):" << GetNumBetweenC1andC2(xopt, 0.0005, 0.1) << std::endl;
	//std::cout << "[0.1:1)    :" << GetNumBetweenC1andC2(xopt, 0.1, 1) << std::endl;
	//	delete RCGsolver;
	//}
};

integer GetNumBetweenC1andC2(const Element *x, double c1, double c2)
{
	const double *xptr = x->ObtainReadData();
	integer length = x->Getlength();
	integer result = 0;
	for (integer i = 0; i < length; i++)
	{
		if (fabs(xptr[i]) >= c1 && fabs(xptr[i]) < c2)
		{
			result++;
		}
	}
	return result;
};

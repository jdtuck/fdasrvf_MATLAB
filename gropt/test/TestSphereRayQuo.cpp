
#include "TestSphereRayQuo.h"

#if !defined(MATLAB_MEX_FILE) && defined(TESTSPHERERAYQUO)

std::map<integer *, integer> *CheckMemoryDeleted;

int main(void)
{
	init_genrand(0);

	// size of the Sphere
	integer n = 12;

	// Generate the matrices in the problem.
	double *B = new double[n * n + 1];
	double *D = B + n * n;
	for (integer i = 0; i < n; i++)
	{
		for (integer j = i; j < n; j++)
		{
			B[i + j * n] = genrand_gaussian();
			B[j + i * n] = B[i + j * n];
		}
	}
	D[0] = 1;

	CheckMemoryDeleted = new std::map<integer *, integer>;

	testSphereRayQuo(B, D, n);
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
    if(nrhs < 2)
    {
        mexErrMsgTxt("The number of arguments should be at least two.\n");
    }
    double Dv = 1;
	double *B, *D = &Dv, *X, *Xopt;
	B = mxGetPr(prhs[0]);
	X = mxGetPr(prhs[1]);
    
	/* dimensions of input matrices */
	integer n;
	n = mxGetM(prhs[0]);
    
    if(mxGetN(prhs[0]) != n)
    {
        mexErrMsgTxt("The size of matrix is not correct.\n");
    }
    if(mxGetM(prhs[1]) != n || mxGetN(prhs[1]) != 1)
    {
        mexErrMsgTxt("The size of the initial X is not correct!\n");
    }
    
	std::cout << "n:" << n << std::endl;

	/*create output matrix*/
	plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
	Xopt = mxGetPr(plhs[0]);

	init_genrand(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;
	testSphereRayQuo(B, D, n, X, Xopt);
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

void testSphereRayQuo(double *B, double *D, integer n, double *X, double *Xopt)
{
	// choose a random seed
	unsigned tt = (unsigned)time(NULL);
	tt = 0;
	init_genrand(tt);

	// Obtain an initial iterate by taking the Q factor of qr decomposition
	// lapack is used
	SphereVariable SphereX(n);
	if (X == nullptr)
	{
		SphereX.RandInManifold();
	}
	else
	{
		double *SphereXptr = SphereX.ObtainWriteEntireData();
		for (integer i = 0; i < n; i++)
			SphereXptr[i] = X[i];
	}

	// Generate the matrices in the Brockett problem with p = 1
	// In this case, it is the Rayleigh Quotient problem on the Sphere.

	// Define the manifold
	Sphere Domain(n);
	Domain.ChooseSphereParamsSet3();

	// Define the Brockett problem with p = 1
	// In this case, it is the Rayleigh Quotient problem on the Sphere.
	StieBrockett Prob(B, D, n, 1);
	Prob.SetDomain(&Domain);

	Domain.CheckParams();

	//Domain.CheckIntrExtr(&SphereX);
	//Domain.CheckRetraction(&SphereX);
	//Domain.CheckDiffRetraction(&SphereX);
	//Domain.CheckLockingCondition(&SphereX);
	//Domain.CheckcoTangentVector(&SphereX);
	//Domain.CheckIsometryofVectorTransport(&SphereX);
	//Domain.CheckIsometryofInvVectorTransport(&SphereX);
	//Domain.CheckVecTranComposeInverseVecTran(&SphereX);
	//Domain.CheckTranHInvTran(&SphereX);
	//Domain.CheckHaddScaledRank1OPE(&SphereX);

	// test RSD
	std::cout << "********************************Check all line search algorithm in RSD*****************************************" << std::endl;
	for (integer i = 0; i < LSALGOLENGTH; i++)
	{
		RSD *RSDsolver = new RSD(&Prob, &SphereX);
		RSDsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		RSDsolver->DEBUG = FINALRESULT;
		RSDsolver->CheckParams();
		RSDsolver->Run();
		delete RSDsolver;
	}

	// test RNewton
	std::cout << "********************************Check all line search algorithm in RNewton*************************************" << std::endl;
	for (integer i = 0; i < LSALGOLENGTH; i++)
	{
		RNewton *RNewtonsolver = new RNewton(&Prob, &SphereX);
		RNewtonsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		RNewtonsolver->DEBUG = FINALRESULT;
		RNewtonsolver->CheckParams();
		RNewtonsolver->Run();
		delete RNewtonsolver;
	}

	// test RCG
	std::cout << "********************************Check all Formulas in RCG*************************************" << std::endl;
	for (integer i = 0; i < RCGMETHODSLENGTH; i++)
	{
		RCG *RCGsolver = new RCG(&Prob, &SphereX);
		RCGsolver->RCGmethod = static_cast<RCGmethods> (i);
		RCGsolver->LineSearch_LS = STRONGWOLFE;
		RCGsolver->LS_beta = 0.1;
		RCGsolver->DEBUG = FINALRESULT;
		RCGsolver->CheckParams();
		RCGsolver->Run();
		delete RCGsolver;
	}

	// test RBroydenFamily
	std::cout << "********************************Check all line search algorithm in RBroydenFamily*************************************" << std::endl;
	for (integer i = 0; i < LSALGOLENGTH; i++)
	{
		RBroydenFamily *RBroydenFamilysolver = new RBroydenFamily(&Prob, &SphereX);
		RBroydenFamilysolver->LineSearch_LS = static_cast<LSAlgo> (i);
		RBroydenFamilysolver->DEBUG = FINALRESULT;
		RBroydenFamilysolver->CheckParams();
		RBroydenFamilysolver->Run();
		delete RBroydenFamilysolver;
	}

	// test RWRBFGS
	std::cout << "********************************Check all line search algorithm in RWRBFGS*************************************" << std::endl;
	for (integer i = 0; i < LSALGOLENGTH; i++)
	{
		RWRBFGS *RWRBFGSsolver = new RWRBFGS(&Prob, &SphereX);
		RWRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		RWRBFGSsolver->DEBUG = FINALRESULT; //ITERRESULT;//
		RWRBFGSsolver->CheckParams();
		RWRBFGSsolver->Run();
		delete RWRBFGSsolver;
	}

	// test RBFGS
	std::cout << "********************************Check all line search algorithm in RBFGS*************************************" << std::endl;
	for (integer i = 0; i < LSALGOLENGTH; i++)
	{
		RBFGS *RBFGSsolver = new RBFGS(&Prob, &SphereX);
		RBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		RBFGSsolver->DEBUG = FINALRESULT;
		RBFGSsolver->CheckParams();
		RBFGSsolver->Run();
		delete RBFGSsolver;
	}

	// test LRBFGS
	std::cout << "********************************Check all line search algorithm in LRBFGS*************************************" << std::endl;
	for (integer i = 0; i < LSALGOLENGTH; i++)
	{
		LRBFGS *LRBFGSsolver = new LRBFGS(&Prob, &SphereX);
		LRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		LRBFGSsolver->DEBUG = FINALRESULT;
		LRBFGSsolver->CheckParams();
		LRBFGSsolver->Run();
		delete LRBFGSsolver;
	}

	// test RTRSD
	std::cout << "********************************Check RTRSD*************************************" << std::endl;
	RTRSD RTRSDsolver(&Prob, &SphereX);
	RTRSDsolver.DEBUG = FINALRESULT;
	RTRSDsolver.CheckParams();
	RTRSDsolver.Run();

	// test RTRNewton
	std::cout << "********************************Check RTRNewton*************************************" << std::endl;
	RTRNewton RTRNewtonsolver(&Prob, &SphereX);
	RTRNewtonsolver.DEBUG = FINALRESULT;
	RTRNewtonsolver.CheckParams();
	RTRNewtonsolver.Run();

	// test RTRSR1
	std::cout << "********************************Check RTRSR1*************************************" << std::endl;
	RTRSR1 RTRSR1solver(&Prob, &SphereX);
	RTRSR1solver.DEBUG = FINALRESULT;
	RTRSR1solver.CheckParams();
	RTRSR1solver.Run();

	// test LRTRSR1
	std::cout << "********************************Check LRTRSR1*************************************" << std::endl;
	LRTRSR1 LRTRSR1solver(&Prob, &SphereX);
	LRTRSR1solver.DEBUG = FINALRESULT;
	LRTRSR1solver.CheckParams();
	LRTRSR1solver.Run();

	// Check gradient and Hessian
	Prob.CheckGradHessian(&SphereX);
	const Variable *xopt = RTRNewtonsolver.GetXopt();
	Prob.CheckGradHessian(xopt);
    
	if (Xopt != nullptr)
	{
		const double *xoptptr = xopt->ObtainReadData();
		for (integer i = 0; i < n; i++)
			Xopt[i] = xoptptr[i];
	}
};

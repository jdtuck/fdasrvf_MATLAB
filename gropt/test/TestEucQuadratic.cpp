
#include "TestEucQuadratic.h"

#if !defined(MATLAB_MEX_FILE) && defined(TESTEUCQUADRATIC)

std::map<integer *, integer> *CheckMemoryDeleted;

int main(void)
{
	init_genrand(0);

	// size of the domain
	integer dim = 10;

	// Generate the matrices in the Euclidean Quadratic problem.
	// Use blas to obtain a positive definite matrix by M = Temp * Temp^T
	double *M = new double[dim * dim];
	double *Temp = new double[dim * dim];
	for (integer i = 0; i < dim * dim; i++)
		Temp[i] = genrand_gaussian();
	char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
	double one = 1, zero = 0;
	integer N = dim;
	std::cout << "start" << std::endl;
	dgemm_(transn, transt, &N, &N, &N, &one, Temp, &N, Temp, &N, &zero, M, &N);
	std::cout << "end" << std::endl;

	delete[] Temp;

	CheckMemoryDeleted = new std::map<integer *, integer>;

	testEucQuadratic(M, dim);
	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			std::cout << "Global address:" << iter->first << ", sharedtimes:" << iter->second << std::endl;
	}
	delete CheckMemoryDeleted;
	delete[] M;

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
    
	double *M, *X, *Xopt;
	M = mxGetPr(prhs[0]);
	X = mxGetPr(prhs[1]);
	/* dimensions of input matrices */
	integer dim;
	dim = mxGetM(prhs[0]);
    if(mxGetN(prhs[0]) != dim)
    {
        mexErrMsgTxt("The size of matrix is not correct.\n");
    }
    if(mxGetM(prhs[1]) != dim || mxGetN(prhs[1]) != 1)
    {
        mexErrMsgTxt("The size of the initial X is not correct!\n");
    }
    
	std::cout << "dim:" << dim << std::endl;

	/*create output matrix*/
	plhs[0] = mxCreateDoubleMatrix(dim, 1, mxREAL);
	Xopt = mxGetPr(plhs[0]);

	init_genrand(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;
	testEucQuadratic(M, dim, X, Xopt);
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

void testEucQuadratic(double *M, integer dim, double *X, double *Xopt)
{
	// choose a random seed
	unsigned tt = (unsigned)time(NULL);
	std::cout << "tt:" << tt << std::endl;
	tt = 0;
	init_genrand(tt);

	// Obtain an initial iterate
	EucVariable EucX(dim, 1);
	if (X == nullptr)
	{
		EucX.RandInManifold();
	}
	else
	{
		double *EucXptr = EucX.ObtainWriteEntireData();
		for (integer i = 0; i < dim; i++)
			EucXptr[i] = X[i];
	}

	// Define the manifold
	Euclidean Domain(dim);

	// Define the problem
	EucQuadratic Prob(M, dim);
	Prob.SetDomain(&Domain);

	// test RSD
	std::cout << "********************************Check all line search algorithm in RSD*****************************************" << std::endl;
	for (integer i = 0; i < LSALGOLENGTH; i++)
	{
		RSD *RSDsolver = new RSD(&Prob, &EucX);
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
		RNewton *RNewtonsolver = new RNewton(&Prob, &EucX);
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
		RCG *RCGsolver = new RCG(&Prob, &EucX);
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
		RBroydenFamily *RBroydenFamilysolver = new RBroydenFamily(&Prob, &EucX);
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
		RWRBFGS *RWRBFGSsolver = new RWRBFGS(&Prob, &EucX);
		RWRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		RWRBFGSsolver->DEBUG = FINALRESULT;
		RWRBFGSsolver->CheckParams();
		RWRBFGSsolver->Run();
		delete RWRBFGSsolver;
	}

	// test RBFGS
	std::cout << "********************************Check all line search algorithm in RBFGS*************************************" << std::endl;
	for (integer i = 0; i < LSALGOLENGTH; i++)
	{
		RBFGS *RBFGSsolver = new RBFGS(&Prob, &EucX);
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
		LRBFGS *LRBFGSsolver = new LRBFGS(&Prob, &EucX);
		LRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		LRBFGSsolver->DEBUG = FINALRESULT;
		LRBFGSsolver->CheckParams();
		LRBFGSsolver->Run();
		delete LRBFGSsolver;
	}

	std::cout << "********************************Check RTRSD*************************************" << std::endl;
	RTRSD RTRSDsolver(&Prob, &EucX);
	std::cout << std::endl;
	RTRSDsolver.DEBUG = FINALRESULT;
	RTRSDsolver.CheckParams();
	RTRSDsolver.Run();

	std::cout << "********************************Check RTRNewton*************************************" << std::endl;
	RTRNewton RTRNewtonsolver(&Prob, &EucX);
	std::cout << std::endl;
	RTRNewtonsolver.DEBUG = FINALRESULT;
	RTRNewtonsolver.CheckParams();
	RTRNewtonsolver.Run();

	std::cout << "********************************Check RTRSR1*************************************" << std::endl;
	RTRSR1 RTRSR1solver(&Prob, &EucX);
	std::cout << std::endl;
	RTRSR1solver.DEBUG = FINALRESULT;
	RTRSR1solver.CheckParams();
	RTRSR1solver.Run();

	std::cout << "********************************Check LRTRSR1*************************************" << std::endl;
	LRTRSR1 LRTRSR1solver(&Prob, &EucX);
	std::cout << std::endl;
	LRTRSR1solver.DEBUG = FINALRESULT;
	LRTRSR1solver.CheckParams();
	LRTRSR1solver.Run();

	// Check gradient and Hessian
	Prob.CheckGradHessian(&EucX);
	const Variable *xopt = RTRNewtonsolver.GetXopt();
	Prob.CheckGradHessian(xopt);
    
	if (Xopt != nullptr)
	{
		const double *xoptptr = xopt->ObtainReadData();
		for (integer i = 0; i < dim; i++)
			Xopt[i] = xoptptr[i];
	}
};

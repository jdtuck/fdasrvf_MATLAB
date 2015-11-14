
#include "TestEucFrechetMean.h"

#if !defined(MATLAB_MEX_FILE) && defined(TESTEUCFRECHETMEAN)

std::map<integer *, integer> *CheckMemoryDeleted;

int main(void)
{
	init_genrand(0);

	// size of the samples and number of samples
	integer num = 6;
	integer dim = 4;

	// Generate the matrices in the Euclidean Frechetmean problem.
	double *Weights, *Data;
	double sumW;
	Weights = new double[num + num * dim];
	Data = Weights + num;
	for (integer i = 0; i < num; i++)
		Weights[i] = genrand_real1() * 1e2 + 1e-4;
	Weights[0] = 1e-2;
	Weights[1] = 1e4;
	sumW = 0;
	for (integer i = 0; i < num; i++)
		sumW += Weights[i];
	for (integer i = 0; i < num; i++)
		Weights[i] /= sumW;
	for (integer i = 0; i < num * dim; i++)
		Data[i] = genrand_gaussian();

	CheckMemoryDeleted = new std::map<integer *, integer>;

	testEucFrechetMean(Data, Weights, num, dim);
	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			std::cout << "Global address:" << iter->first << ", sharedtimes:" << iter->second << std::endl;
	}
	delete CheckMemoryDeleted;
	delete[] Weights;

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
    
	double *Data, *Weights, *X, *Xopt;
	Data = mxGetPr(prhs[0]);
	Weights = mxGetPr(prhs[1]);
	X = mxGetPr(prhs[2]);
	/* dimensions of input matrices */
	integer dim, num;
	dim = mxGetM(prhs[0]);
	num = mxGetN(prhs[0]);
    if(mxGetM(prhs[2]) != dim || mxGetN(prhs[2]) != 1)
    {
        mexErrMsgTxt("The size of initial X is not correct!\n");
    }
    if(mxGetM(prhs[1]) != num || mxGetN(prhs[1]) != 1)
    {
        mexErrMsgTxt("The size of the weights is not correct!\n");
    }
    
	std::cout << "(dim, num):" << dim << "," << num << std::endl;

	/*create output matrix*/
	plhs[0] = mxCreateDoubleMatrix(dim, 1, mxREAL);
	Xopt = mxGetPr(plhs[0]);

	init_genrand(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;
	testEucFrechetMean(Data, Weights, num, dim, X, Xopt);
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

void testEucFrechetMean(double *Data, double *Weights, integer num, integer dim, double *X, double *Xopt)
{
	// choose a random seed
	init_genrand(0);// (unsigned)time(NULL));


	// Obtain an initial iterate
	EucVariable EucX(dim);
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
	EucFrechetMean Prob(Weights, Data, num, dim);
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

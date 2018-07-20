
#include "TestStieBrockett.h"

#if !defined(MATLAB_MEX_FILE) && defined(TESTSTIEBROCKETT)

std::map<integer *, integer> *CheckMemoryDeleted;

int main(void)
{
	init_genrand(0);

	// size of the Stiefel manifold
	integer n = 12, p = 8;

	// Generate the matrices in the Brockett problem.
	double *B = new double[n * n + p];
	double *D = B + n * n;
	for (integer i = 0; i < n; i++)
	{
		for (integer j = i; j < n; j++)
		{
			B[i + j * n] = genrand_gaussian();
			B[j + i * n] = B[i + j * n];
		}
	}
	for (integer i = 0; i < p; i++)
		D[i] = static_cast<double> (i + 1);

	CheckMemoryDeleted = new std::map<integer *, integer>;

	testStieBrockett(B, D, n, p);
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

void testStieBrockett(double *B, double *D, integer n, integer p, double *X, double *Xopt)
{
	// choose a random seed
	unsigned tt = (unsigned)time(NULL);
	tt = 0;
	init_genrand(tt);

	// Obtain an initial iterate by taking the Q factor of qr decomposition
	// lapack is used
	StieVariable StieX(n, p);
	if (X == nullptr)
	{
		StieX.RandInManifold();
	}
	else
	{
		double *StieXptr = StieX.ObtainWriteEntireData();
		for (integer i = 0; i < n * p; i++)
			StieXptr[i] = X[i];
	}

	// Define the manifold
	Stiefel Domain(n, p);

	// Define the Brockett problem
	StieBrockett Prob(B, D, n, p);
	Prob.SetDomain(&Domain);

	Domain.CheckParams();

	//Domain.CheckIntrExtr(&StieX);
	//Domain.CheckRetraction(&StieX);
	//Domain.CheckDiffRetraction(&StieX);
	//Domain.CheckLockingCondition(&StieX);
	//Domain.CheckcoTangentVector(&StieX);
	//Domain.CheckIsometryofVectorTransport(&StieX);
	//Domain.CheckIsometryofInvVectorTransport(&StieX);
	//Domain.CheckVecTranComposeInverseVecTran(&StieX);
	//Domain.CheckTranHInvTran(&StieX);
	//Domain.CheckHaddScaledRank1OPE(&StieX);

	//// Check gradient and Hessian
	//Prob.CheckGradHessian(&StieX);

	//// test RSD
	//std::cout << "********************************Check all line search algorithm in RSD*****************************************" << std::endl;
	//for (integer i = 0; i < LSALGOLENGTH; i++)
	//{
	//	RSD *RSDsolver = new RSD(&Prob, &StieX);
	//	RSDsolver->LineSearch_LS = static_cast<LSAlgo> (i);
	//	RSDsolver->DEBUG = FINALRESULT;
	//	RSDsolver->CheckParams();
	//	RSDsolver->Run();
	//	delete RSDsolver;
	//}

	//// test RNewton
	//std::cout << "********************************Check all line search algorithm in RNewton*************************************" << std::endl;
	//for (integer i = 0; i < LSALGOLENGTH; i++)
	//{
	//	RNewton *RNewtonsolver = new RNewton(&Prob, &StieX);
	//	RNewtonsolver->LineSearch_LS = static_cast<LSAlgo> (i);
	//	RNewtonsolver->DEBUG = FINALRESULT;
	//	RNewtonsolver->CheckParams();
	//	RNewtonsolver->Run();
	//	delete RNewtonsolver;
	//}

	//// test RCG
	//std::cout << "********************************Check all Formulas in RCG*************************************" << std::endl;
	//for (integer i = 0; i < RCGMETHODSLENGTH; i++)
	//{
	//	RCG *RCGsolver = new RCG(&Prob, &StieX);
	//	RCGsolver->RCGmethod = static_cast<RCGmethods> (i);
	//	RCGsolver->LineSearch_LS = STRONGWOLFE;
	//	RCGsolver->LS_beta = 0.1;
	//	RCGsolver->DEBUG = FINALRESULT;
	//	RCGsolver->CheckParams();
	//	RCGsolver->Run();
	//	delete RCGsolver;
	//}

	//// test RBroydenFamily
	//std::cout << "********************************Check all line search algorithm in RBroydenFamily*************************************" << std::endl;
	//for (integer i = 0; i < LSALGOLENGTH; i++)
	//{
	//	RBroydenFamily *RBroydenFamilysolver = new RBroydenFamily(&Prob, &StieX);
	//	RBroydenFamilysolver->LineSearch_LS = static_cast<LSAlgo> (i);
	//	RBroydenFamilysolver->DEBUG = FINALRESULT;
	//	RBroydenFamilysolver->CheckParams();
	//	RBroydenFamilysolver->Run();
	//	delete RBroydenFamilysolver;
	//}

	//// test RWRBFGS
	//std::cout << "********************************Check all line search algorithm in RWRBFGS*************************************" << std::endl;
	//for (integer i = 0; i < LSALGOLENGTH; i++)
	//{
	//	RWRBFGS *RWRBFGSsolver = new RWRBFGS(&Prob, &StieX);
	//	RWRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
	//	RWRBFGSsolver->DEBUG = FINALRESULT; //ITERRESULT;//
	//	RWRBFGSsolver->CheckParams();
	//	RWRBFGSsolver->Run();
	//	delete RWRBFGSsolver;
	//}

	//// test RBFGS
	//std::cout << "********************************Check all line search algorithm in RBFGS*************************************" << std::endl;
	//for (integer i = 0; i < LSALGOLENGTH; i++)
	//{
	//	RBFGS *RBFGSsolver = new RBFGS(&Prob, &StieX);
	//	RBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
	//	RBFGSsolver->DEBUG = FINALRESULT;
	//	RBFGSsolver->CheckParams();
	//	RBFGSsolver->Run();
	//	delete RBFGSsolver;
	//}

	// test LRBFGS
	std::cout << "********************************Check all line search algorithm in LRBFGS*************************************" << std::endl;
	for (integer i = 0; i < 1; i++)//LSALGOLENGTH
	{
		LRBFGS *LRBFGSsolver = new LRBFGS(&Prob, &StieX);
		LRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		LRBFGSsolver->DEBUG = FINALRESULT; //ITERRESULT;// 
		LRBFGSsolver->CheckParams();
		LRBFGSsolver->Run();
		delete LRBFGSsolver;
	}

	//// test RTRSD
	//std::cout << "********************************Check RTRSD*************************************" << std::endl;
	//RTRSD RTRSDsolver(&Prob, &StieX);
	//RTRSDsolver.DEBUG = FINALRESULT;
	//RTRSDsolver.CheckParams();
	//RTRSDsolver.Run();

	// test RTRNewton
	std::cout << "********************************Check RTRNewton*************************************" << std::endl;
	RTRNewton RTRNewtonsolver(&Prob, &StieX);
	RTRNewtonsolver.DEBUG = FINALRESULT;
	RTRNewtonsolver.CheckParams();
	RTRNewtonsolver.Run();
	//
	//// test RTRSR1
	//std::cout << "********************************Check RTRSR1*************************************" << std::endl;
	//RTRSR1 RTRSR1solver(&Prob, &StieX);
	//RTRSR1solver.DEBUG = FINALRESULT;
	//RTRSR1solver.CheckParams();
	//RTRSR1solver.Run();

	//// test LRTRSR1
	//std::cout << "********************************Check LRTRSR1*************************************" << std::endl;
	//LRTRSR1 LRTRSR1solver(&Prob, &StieX);
	//LRTRSR1solver.DEBUG = FINALRESULT;
	//LRTRSR1solver.CheckParams();
	//LRTRSR1solver.Run();

	// Check gradient and Hessian
	Prob.CheckGradHessian(&StieX);
	const Variable *xopt = RTRNewtonsolver.GetXopt();
	Prob.CheckGradHessian(xopt);
    
	//const Variable *xopt = RTRNewtonsolver.GetXopt();
	if (Xopt != nullptr)
	{
		const double *xoptptr = xopt->ObtainReadData();
		for (integer i = 0; i < n * p; i++)
			Xopt[i] = xoptptr[i];
	}
}

#endif

#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs < 5)
    {
        mexErrMsgTxt("The number of arguments should be at least five.\n");
    }
	double *B, *D, *X, *Xopt;
	B = mxGetPr(prhs[0]);
	D = mxGetPr(prhs[1]);
	X = mxGetPr(prhs[2]);
	/* dimensions of input matrices */
	integer p, n, HasHHR;
	n = mxGetM(prhs[0]);
	p = mxGetM(prhs[1]);
    
    if(mxGetN(prhs[0]) != n)
    {
        mexErrMsgTxt("The size of matrix B is not correct.\n");
    }
    if(mxGetN(prhs[1]) != 1)
    {
        mexErrMsgTxt("The size of the D is not correct!\n");
    }
    if(mxGetM(prhs[2]) != n || mxGetN(prhs[2]) != p)
    {
        mexErrMsgTxt("The size of the initial X is not correct!\n");
    }
	HasHHR = static_cast<integer> (mxGetScalar(prhs[3]));
    
	std::cout << "(n, p):" << n << "," << p << std::endl;

	/*create output matrix*/
	plhs[0] = mxCreateDoubleMatrix(n, p, mxREAL);
	Xopt = mxGetPr(plhs[0]);

	init_genrand(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;
//	testStieBrockett(B, D, n, p, X, Xopt);

	// Obtain an initial iterate by taking the Q factor of qr decomposition
	StieVariable StieX(n, p);
	double *StieXptr = StieX.ObtainWriteEntireData();
	for (integer i = 0; i < n * p; i++)
		StieXptr[i] = X[i];

	// Define the manifold
	Stiefel Domain(n, p);

	// Define the Brockett problem
	StieBrockett Prob(B, D, n, p);
	Prob.SetDomain(&Domain);

	Domain.SetHasHHR(HasHHR != 0);
	//Domain.CheckParams();
	ParseSolverParamsAndOptimizing(prhs[4], &Prob, &StieX, plhs);
    
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

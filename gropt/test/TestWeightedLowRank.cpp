#include "TestWeightedLowRank.h"

#if !defined(MATLAB_MEX_FILE) && defined(TESTWEIGHTEDLOWRANK)

std::map<integer *, integer> *CheckMemoryDeleted;

int main(void)
{
	long seed = static_cast<long> (time(NULL));
	//seed = 1417791199;//---
	seed = 0;
	std::cout << "seed:" << seed << std::endl;
	init_genrand(seed);
	
	CheckMemoryDeleted = new std::map<integer *, integer>;
    
	testWeightedLowRank();
	
	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			std::cout << "Global address:" << iter->first << ", sharedtimes:" << iter->second << std::endl;
	}
	delete CheckMemoryDeleted;
	
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
	init_genrand(0);
	
	CheckMemoryDeleted = new std::map<integer *, integer>;
	testWeightedLowRank();
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

void testWeightedLowRank(void)
{
	integer m = 5, n = 4, r = 2;
	LowRank Domain(m, n, r);
	Domain.SetHasHHR(true);
	LowRankVariable InitialX(m, n, r);
	InitialX.RandInManifold();
	//Domain.CheckParams();
	//Domain.CheckIntrExtr(&InitialX);
	//Domain.CheckRetraction(&InitialX);
	//Domain.CheckcoTangentVector(&InitialX);
	//Domain.CheckDiffRetraction(&InitialX, false);
	//Domain.CheckIsometryofVectorTransport(&InitialX);

	//Domain.CheckLockingCondition(&InitialX);
	//Domain.CheckIsometryofInvVectorTransport(&InitialX);
	//Domain.CheckVecTranComposeInverseVecTran(&InitialX);
	//Domain.CheckTranHInvTran(&InitialX);

	//InitialX.Print("initialX:");
	//return;
	// Generate the matrices in the Low rank approximation problem.
    integer mn = m * n, mr = m * r, nr = n * r;
	double *A = new double[mn];
    double *A_U = new double[mr];
	double *A_V = new double[nr];
    for (integer i = 0; i < m * r; i++)
    {
       A_U[i] = genrand_gaussian();
    }
	for (integer i = 0; i < n * r; i++)
	{
		A_V[i] = genrand_gaussian();
	}
	char *transn = const_cast<char *> ("n");
  	char *transt = const_cast<char *> ("t");
	double one = 1, zero = 0;
	dgemm_(transn, transt, &m, &n, &r, &one, A_U, &m, A_V, &n, &zero, A, &m);
	delete []A_U;
	delete []A_V;
    
    double *W_temp = new double[mn * mn];
	double *W = new double[mn * mn];
    for (integer i = 0; i < mn; i++)
    {
        for (integer j = 0; j < mn; j++)
        {
            W_temp[i + j * mn] = genrand_gaussian();
			//if (i == j)
			//{
			//	W_temp[i + i * mn] = 1;
			//}
			//else
			//{
			//	W_temp[i + j * mn] = 0;
			//}
        }
    }
    dgemm_(transn, transt, &mn, &mn, &mn, &one, W_temp, &mn, W_temp, &mn, &zero, W, &mn);
	delete []W_temp;
	
	Stiefel mani1(m, r);
	Euclidean mani2(r, r);
	Stiefel mani3(n, r);
	
    WeightedLowRank Prob(A, W, m, n, r);
    Prob.SetDomain(&Domain);
	
	//Prob.CheckGradHessian(&InitialX);
    
	LRBFGS *RSDsolver = new LRBFGS(&Prob, &InitialX);
	RSDsolver->LineSearch_LS = ARMIJO;
	//RSDsolver->LS_beta = 0.01;
	//RSDsolver->RCGmethod = DAI_YUAN;
	RSDsolver->DEBUG = ITERRESULT;
	RSDsolver->OutputGap = 100;
	RSDsolver->Max_Iteration = 10000;
	RSDsolver->CheckParams();
	RSDsolver->Run();
	Prob.CheckGradHessian(&InitialX);//--
	Prob.CheckGradHessian(RSDsolver->GetXopt());//--

	delete RSDsolver;
	
	delete []A;
    delete []W;
};



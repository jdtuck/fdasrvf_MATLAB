
#ifndef TESTSIMPLEEXAMPLE_CPP
#define TESTSIMPLEEXAMPLE_CPP

#include "ForDebug.h"
#include <iostream>
#include "randgen.h"
#include "Manifold.h"
#include "Problem.h"
#include "SolversLS.h"
#include <ctime>

#include "EucVariable.h"
#include "EucVector.h"
#include "EucFrechetMean.h"
#include "EucQuadratic.h"

#include "StieBrockett.h"
#include "StieVector.h"
#include "StieVariable.h"
#include "Stiefel.h"

#include "RSD.h"
#include "RNewton.h"
#include "RCG.h"
#include "RBroydenFamily.h"
#include "RWRBFGS.h"
#include "RBFGS.h"
#include "LRBFGS.h"

#include "SolversTR.h"
#include "RTRSD.h"
#include "RTRNewton.h"
#include "RTRSR1.h"
#include "LRTRSR1.h"

#include "def.h"

#ifdef MATLAB_MEX_FILE
#define TESTSIMPLEEXAMPLE
int main(void);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    main();
	return;
}
#endif

#ifdef TESTSIMPLEEXAMPLE
int main(void)
{
	// choose a random seed
	unsigned tt = (unsigned)time(NULL);
	init_genrand(tt);

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

	// Obtain an initial iterate
	StieVariable StieX(n, p);
	StieX.RandInManifold();

	// Define the Stiefel manifold
	Stiefel Domain(n, p);

	// Define the Brockett problem
	StieBrockett Prob(B, D, n, p);

	// Set the domain of the problem to be the Stiefel manifold
	Prob.SetDomain(&Domain);

	// output the parameters of the manifold of domain
	Domain.CheckParams();


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

	//// test LRBFGS
	//std::cout << "********************************Check all line search algorithm in LRBFGS*************************************" << std::endl;
	//for (integer i = 0; i < 1; i++)//LSALGOLENGTH
	//{
	//	LRBFGS *LRBFGSsolver = new LRBFGS(&Prob, &StieX);
	//	LRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
	//	LRBFGSsolver->DEBUG = FINALRESULT; //ITERRESULT;// 
	//	LRBFGSsolver->CheckParams();
	//	LRBFGSsolver->Run();
	//	delete LRBFGSsolver;
	//}

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

	delete[] B;

	return 0;
}

#endif // end of TESTSIMPLEEXAMPLE
#endif // end of TESTSIMPLEEXAMPLE_CPP

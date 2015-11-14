
#ifndef TESTPRODUCTEXAMPLE_CPP
#define TESTPRODUCTEXAMPLE_CPP

#include "ForDebug.h"
#include <iostream>
#include "randgen.h"
#include "Manifold.h"
#include "Problem.h"
#include "SolversLS.h"
#include <ctime>

#include "StieSumBrockett.h"
#include "StieVector.h"
#include "StieVariable.h"
#include "Stiefel.h"
#include <ProductElement.h>
#include <ProductManifold.h>

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
#define TESTPRODUCTEXAMPLE
int main(void);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    main();
	return;
}
#endif

#ifdef TESTPRODUCTEXAMPLE

int main(void)
{
	// choose a random seed
	unsigned tt = (unsigned)time(NULL);
	init_genrand(tt);

	// size of the Stiefel manifold
	integer n = 12, p = 8, m = 6, q = 2;

	// Generate the matrices in the Brockett problem.
	double *B1 = new double[n * n * 2 + p * 2 + m * m + q];
	double *B2 = B1 + n * n;
	double *B3 = B2 + n * n;
	double *D1 = B3 + m * m;
	double *D2 = D1 + p;
	double *D3 = D2 + p;

	for (integer i = 0; i < n; i++)
	{
		for (integer j = i; j < n; j++)
		{
			B1[i + j * n] = genrand_gaussian();
			B1[j + i * n] = B1[i + j * n];

			B2[i + j * n] = genrand_gaussian();
			B2[j + i * n] = B2[i + j * n];
		}
	}
	for (integer i = 0; i < m; i++)
	{
		for (integer j = i; j < m; j++)
		{
			B3[i + j * m] = genrand_gaussian();
			B3[j + i * m] = B3[i + j * m];
		}
	}
	for (integer i = 0; i < p; i++)
	{
		D1[i] = static_cast<double> (i + 1);
		D2[i] = D1[i];
	}
	for (integer i = 0; i < q; i++)
	{
		D3[i] = static_cast<double> (i + 1);
	}

	// number of manifolds in product of manifold
	integer numofmanis = 2; // two kinds of manifolds
	integer numofmani1 = 2; // the first one has two
	integer numofmani2 = 1; // the second one has one

	// Obtain an initial iterate
	StieVariable StieX1(n, p);
	StieVariable StieX2(m, q);
	ProductElement ProdX(numofmanis, &StieX1, numofmani1, &StieX2, numofmani2);
	ProdX.RandInManifold();

	// Define the Stiefel manifold
	Stiefel mani1(n, p);
	Stiefel mani2(m, q);
	ProductManifold Domain(numofmanis, &mani1, numofmani1, &mani2, numofmani2);

	// Define the Brockett problem
	StieSumBrockett Prob(B1, D1, B2, D2, B3, D3, n, p, m, q);

	// Set the domain of the problem to be the Stiefel manifold
	Prob.SetDomain(&Domain);

	// output the parameters of the manifold of domain
	Domain.CheckParams();

	//// test RSD
	//std::cout << "********************************Check all line search algorithm in RSD*****************************************" << std::endl;
	//for (integer i = 0; i < LSALGOLENGTH; i++)
	//{
	//	RSD *RSDsolver = new RSD(&Prob, &ProdX);
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
	//	RNewton *RNewtonsolver = new RNewton(&Prob, &ProdX);
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
	//	RCG *RCGsolver = new RCG(&Prob, &ProdX);
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
	//	RBroydenFamily *RBroydenFamilysolver = new RBroydenFamily(&Prob, &ProdX);
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
	//	RWRBFGS *RWRBFGSsolver = new RWRBFGS(&Prob, &ProdX);
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
	//	RBFGS *RBFGSsolver = new RBFGS(&Prob, &ProdX);
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
	//	LRBFGS *LRBFGSsolver = new LRBFGS(&Prob, &ProdX);
	//	LRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
	//	LRBFGSsolver->DEBUG = FINALRESULT; //ITERRESULT;// 
	//	LRBFGSsolver->CheckParams();
	//	LRBFGSsolver->Run();
	//	delete LRBFGSsolver;
	//}

	//// test RTRSD
	//std::cout << "********************************Check RTRSD*************************************" << std::endl;
	//RTRSD RTRSDsolver(&Prob, &ProdX);
	//RTRSDsolver.DEBUG = FINALRESULT;
	//RTRSDsolver.CheckParams();
	//RTRSDsolver.Run();

	// test RTRNewton
	std::cout << "********************************Check RTRNewton*************************************" << std::endl;
	RTRNewton RTRNewtonsolver(&Prob, &ProdX);
	RTRNewtonsolver.DEBUG = FINALRESULT;
	RTRNewtonsolver.CheckParams();
	RTRNewtonsolver.Run();

	//// test RTRSR1
	//std::cout << "********************************Check RTRSR1*************************************" << std::endl;
	//RTRSR1 RTRSR1solver(&Prob, &ProdX);
	//RTRSR1solver.DEBUG = FINALRESULT;
	//RTRSR1solver.CheckParams();
	//RTRSR1solver.Run();

	//// test LRTRSR1
	//std::cout << "********************************Check LRTRSR1*************************************" << std::endl;
	//LRTRSR1 LRTRSR1solver(&Prob, &ProdX);
	//LRTRSR1solver.DEBUG = FINALRESULT;
	//LRTRSR1solver.CheckParams();
	//LRTRSR1solver.Run();

	// Check gradient and Hessian
	Prob.CheckGradHessian(&ProdX);
	const Variable *xopt = RTRNewtonsolver.GetXopt();
	Prob.CheckGradHessian(xopt);

	delete[] B1;

	return 0;
}

#endif // end of TESTPRODUCTEXAMPLE
#endif // end of TESTPRODUCTEXAMPLE_CPP

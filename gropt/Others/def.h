#ifndef DEF_H
#define DEF_H

#define INITIALBLOCKSIZE 64

#define MATLAB_MEX_FILE //-------

	// test examples
//#define TESTEUCFRECHETMEAN
//#define TESTEUCQUADRATIC
//#define TESTPRODUCT
//#define TESTSPHERERAYQUO
#define TESTSTIEBROCKETT
//#define TESTSTIESOFTICA
//#define TESTTESTSPARSEPCA
//#define TESTWEIGHTEDLOWRANK
//#define TESTELASTICCURVESRO

//#define TESTSIMPLEEXAMPLE
//#define TESTPRODUCTEXAMPLE

#ifndef MATLAB_MEX_FILE
    // blas and lapack related
    #include <dgemm.h>
    #include <dgetrf.h>
    #include <dgetrs.h>
    #include <dgemv.h>
    #include <dcopy.h>
    #include <ddot.h>
    #include <dscal.h>
    #include <daxpy.h>
    #include <dger.h>
    #include <dgeqp3.h>
    #include <dorgqr.h>
    #include <dormqr.h>
    #include <dtrsm.h>
    #include <dlarfx.h>
	#include <ddot.h>
	#include <dgesdd.h>
	#include <dgesvd.h>
	#include <dsymv.h>
	#include <dgetri.h>
	#include <dlapmt.h>
	#include <dgees.h>
#endif // end of ifndef MATLAB_MEX_FILE

#ifdef _WIN64
      //define something for Windows (64-bit only)
// Test memory leaking
#ifdef _DEBUG
#define DEBUG_CLIENTBLOCK   new( _CLIENT_BLOCK, __FILE__, __LINE__)
//#define CHECKMEMORYDELETED
#else
#define DEBUG_CLIENTBLOCK
#endif

#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>

#ifdef _DEBUG
#define new DEBUG_CLIENTBLOCK
#endif

#elif _WIN32
   //define something for Windows (32-bit and 64-bit, this part is common)
#elif __APPLE__
    #include "TargetConditionals.h"
    #if TARGET_IPHONE_SIMULATOR
         // iOS Simulator
    #elif TARGET_OS_IPHONE
        // iOS device
    #elif TARGET_OS_MAC
        // Other kinds of Mac OS
    #else
        // Unsupported platform
    #endif
#elif __linux
    // linux
#elif __unix // all unices not caught above
    // Unix
#elif __posix
    // POSIX
#endif // end of checking platforms


#ifdef MATLAB_MEX_FILE
	#include "mex.h"
    #include "blas.h"
    #include "lapack.h"
    #define integer ptrdiff_t
#define dgemm_ dgemm
#define dgetrf_ dgetrf
#define dgetrs_ dgetrs
#define dgemv_ dgemv
#define dcopy_ dcopy
#define ddot_ ddot
#define dscal_ dscal
#define daxpy_ daxpy
#define dger_ dger
#define dgeqp3_ dgeqp3
#define dorgqr_ dorgqr
#define dormqr_ dormqr
#define dtrsm_ dtrsm
#define dlarfx_ dlarfx
#define dgesdd_ dgesdd
#define dgesvd_ dgesvd
#define dsymv_ dsymv
#define dgetri_ dgetri
#define dgees_ dgees
#endif // end of ifdef MATLAB_MEX_FILE

#include "ForDebug.h"
#include <climits>
#include <limits>
#include "Timer.h"

#ifdef __INTEL_COMPILER
const class {
public:
	template<class T> // convertible to any type
	operator T*(void) const // of null non-member
	{
		return 0;
	} // pointer...
	template<class C, class T> // or any type of null
	operator T C::*(void) const // member pointer...
	{
		return 0;
	}
private:
	void operator&(void) const; // whose address can't be taken
} nullptr = {};
#endif // end of __GNUC__

#include <map>
#include <string>

typedef std::map<std::string, double> PARAMSMAP;

#define PI 3.14159265358979323846264

#endif // end of DEF_H

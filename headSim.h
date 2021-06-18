#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <math.h>
#include <time.h>
#include <vector>
#include <assert.h>
#include <stdio.h>
#include <ctime>
#include <cstdlib> // for exit function
#include <tr1/random>
#include <numeric>
			

#ifndef RANDOMC_H
#define RANDOMC_H


const int TL=11;
const int mA=50;

using namespace std;
using namespace std::tr1;

const double PI=3.14159265;
const double TWOPI=6.28318531;
const long long int LONG_MAX=3147483647000;


struct totals
	{
	int Jww[TL]; int fJww[TL]; int mJww[TL]; int bJww[TL]; int fJwd[TL]; int mJwd[TL]; int bJwd[TL]; int Jdd[TL]; int Jwr[TL]; int fJwr[TL]; int mJwr[TL]; int bJwr[TL]; int Jrr[TL]; int Jdr[TL]; int mbJrr[TL]; int mbJdr[TL];
	int Mww[mA]; int Mwd[mA]; int Mdd[mA]; int Mwr[mA]; int Mrr[mA]; int Mdr[mA];
	int Vww[mA]; int fVww[mA]; int mVww[mA]; int bVww[mA]; int fVwd[mA]; int mVwd[mA]; int bVwd[mA]; int Vdd[mA]; int Vwr[mA]; int fVwr[mA]; int mVwr[mA];  int bVwr[mA];  int Vrr[mA]; int Vdr[mA];
	int Fwwww[mA]; int Fwwwd[mA]; int Fwwdd[mA]; int Fwwwr[mA]; int Fwwrr[mA]; int Fwwdr[mA];
	int fFwwww[mA]; int fFwwwd[mA]; int fFwwdd[mA]; int fFwwwr[mA]; int fFwwrr[mA]; int fFwwdr[mA];
	int mFwwww[mA]; int mFwwwd[mA]; int mFwwdd[mA]; int mFwwwr[mA]; int mFwwrr[mA]; int mFwwdr[mA];
	int bFwwww[mA]; int bFwwwd[mA]; int bFwwdd[mA]; int bFwwwr[mA]; int bFwwrr[mA]; int bFwwdr[mA];
	int fFwdww[mA]; int fFwdwd[mA]; int fFwddd[mA]; int fFwdwr[mA]; int fFwdrr[mA]; int fFwddr[mA];
	int mFwdww[mA]; int mFwdwd[mA]; int mFwddd[mA]; int mFwdwr[mA]; int mFwdrr[mA]; int mFwddr[mA];
	int bFwdww[mA]; int bFwdwd[mA]; int bFwddd[mA]; int bFwdwr[mA]; int bFwdrr[mA]; int bFwddr[mA];
	int Fddww[mA]; int Fddwd[mA]; int Fdddd[mA]; int Fddwr[mA]; int Fddrr[mA]; int Fdddr[mA];
	int fFwrww[mA]; int fFwrwd[mA]; int fFwrdd[mA]; int fFwrwr[mA]; int fFwrrr[mA]; int fFwrdr[mA];
	int mFwrww[mA]; int mFwrwd[mA]; int mFwrdd[mA]; int mFwrwr[mA]; int mFwrrr[mA]; int mFwrdr[mA];
	int bFwrww[mA]; int bFwrwd[mA]; int bFwrdd[mA]; int bFwrwr[mA]; int bFwrrr[mA]; int bFwrdr[mA];
	int Fwrww[mA]; int Fwrwd[mA]; int Fwrdd[mA]; int Fwrwr[mA]; int Fwrrr[mA]; int Fwrdr[mA];
	int Frrww[mA]; int Frrwd[mA]; int Frrdd[mA]; int Frrwr[mA]; int Frrrr[mA]; int Frrdr[mA];
	int Fdrww[mA]; int Fdrwd[mA]; int Fdrdd[mA]; int Fdrwr[mA]; int Fdrrr[mA]; int Fdrdr[mA];
	int MwwT; int MwdT; int MddT; int MwrT; int MrrT; int MdrT;
	int VwwT; int fVwwT; int mVwwT; int bVwwT; int fVwdT; int mVwdT; int bVwdT; int VddT; int VwrT; int fVwrT; int mVwrT; int bVwrT; int VrrT; int VdrT;
	int FwwT; int fFwwT; int mFwwT; int bFwwT; int fFwdT; int mFwdT; int bFwdT; int FddT; int FwrT; int fFwrT; int mFwrT; int bFwrT; int FrrT; int FdrT;
	int JTot,MTot,VTot,FTot;
	double mate_rate;
	int neggs,keepeggs,r2;
	int toFwwww, toFwwwd, toFwwdd, toFwwwr, toFwwrr, toFwwdr;
	int tofFwwww, tofFwwwd, tofFwwdd, tofFwwwr, tofFwwrr, tofFwwdr;
	int tomFwwww, tomFwwwd, tomFwwdd, tomFwwwr, tomFwwrr, tomFwwdr;
	int tobFwwww, tobFwwwd, tobFwwdd, tobFwwwr, tobFwwrr, tobFwwdr;
	int tofFwdww, tofFwdwd, tofFwddd, tofFwdwr, tofFwdrr, tofFwddr;
	int tomFwdww, tomFwdwd, tomFwddd, tomFwdwr, tomFwdrr, tomFwddr;
	int tobFwdww, tobFwdwd, tobFwddd, tobFwdwr, tobFwdrr, tobFwddr;
	int toFddww, toFddwd, toFdddd, toFddwr, toFddrr, toFdddr;
	int toFwrww, toFwrwd, toFwrdd, toFwrwr, toFwrrr, toFwrdr;
	int tofFwrww, tofFwrwd, tofFwrdd, tofFwrwr, tofFwrrr, tofFwrdr;
	int tomFwrww, tomFwrwd, tomFwrdd, tomFwrwr, tomFwrrr, tomFwrdr;
	int tobFwrww, tobFwrwd, tobFwrdd, tobFwrwr, tobFwrrr, tobFwrdr;
	int toFrrww, toFrrwd, toFrrdd, toFrrwr, toFrrrr, toFrrdr;
	int toFdrww, toFdrwd, toFdrdd, toFdrwr, toFdrrr, toFdrdr;
	};	



struct initials
{
	int startNum;
	int driver_time;
	int NumDriver;
	};	


struct Times
{
	int interval;
	int maxT;
	int NumRuns;
};

struct Pars
	{
	double numDsx;
	double Fgamma,Mgamma,beta,theta,delta,em,ef;
	double omegaww,omegafww,omegamww,omegabww,omegafwd,omegamwd,omegabwd,omegawr,omegafwr,omegamwr,omegabwr,omegadd,omegarr,omegadr;
	double play,xi,Frho,Mrho;
	double lamM[2000];double lamF[2000];double kM[2000];double kF[2000];
	double lamMStart,lamFStart,kMStart,kFStart;
	double lamMEnd,lamFEnd,kMEnd,kFEnd;
	int NumEggs,NeControl;
	int densdepend;
	char control;
	double FracEggs;
	double fwwww[14]; double fwwwd[14]; double fwwdd[14]; double fwwwr[14]; double fwwrr[14]; double fwwdr[14];
	double ffwwww[14]; double ffwwwd[14]; double ffwwdd[14]; double ffwwwr[14]; double ffwwrr[14]; double ffwwdr[14];
	double mfwwww[14]; double mfwwwd[14]; double mfwwdd[14]; double mfwwwr[14]; double mfwwrr[14]; double mfwwdr[14];
	double bfwwww[14]; double bfwwwd[14]; double bfwwdd[14]; double bfwwwr[14]; double bfwwrr[14]; double bfwwdr[14];
	double ffwdww[14]; double ffwdwd[14]; double ffwddd[14]; double ffwdwr[14]; double ffwdrr[14]; double ffwddr[14];
	double mfwdww[14]; double mfwdwd[14]; double mfwddd[14]; double mfwdwr[14]; double mfwdrr[14]; double mfwddr[14];
	double bfwdww[14]; double bfwdwd[14]; double bfwddd[14]; double bfwdwr[14]; double bfwdrr[14]; double bfwddr[14];
	double fddww[14]; double fddwd[14]; double fdddd[14]; double fddwr[14]; double fddrr[14]; double fdddr[14];
	double fwrww[14]; double fwrwd[14]; double fwrdd[14]; double fwrwr[14]; double fwrrr[14]; double fwrdr[14];
	double ffwrww[14]; double ffwrwd[14]; double ffwrdd[14]; double ffwrwr[14]; double ffwrrr[14]; double ffwrdr[14];
	double mfwrww[14]; double mfwrwd[14]; double mfwrdd[14]; double mfwrwr[14]; double mfwrrr[14]; double mfwrdr[14];
	double bfwrww[14]; double bfwrwd[14]; double bfwrdd[14]; double bfwrwr[14]; double bfwrrr[14]; double bfwrdr[14];
	double frrww[14]; double frrwd[14]; double frrdd[14]; double frrwr[14]; double frrrr[14]; double frrdr[14];
	double fdrww[14]; double fdrwd[14]; double fdrdd[14]; double fdrwr[14]; double fdrrr[14]; double fdrdr[14];
	};

		void RunOnceInt(double);
		void record(int);
		void RunMaxT();
		void RunNReps(int);
		void initiate();
		void clearAd();
		void clearJ();
		void OneStep(int);
		void JuvGetOlder();
		void VirginsMate();
		void FemaleCount();
		void LayEggs();
		void JuvEmerge();
		void AdultsDie(int);
		void ComputeTotals();
		void UpdateComp(int);
		void UpdateMate();
		void SetFertility();
		int random_poisson(double);
int randNegBin(double,double);
		int random_binomial(int,double);
		double random_normal(double,double);
		int* random_multinom(int,int[6]);
		int* random_multinom_d(int,double[11]);
		int* random_multinomEqualProb(int,int);
		int* random_multinom_var(int,int,double*,double);


// Define 32 bit signed and unsigned integers.
// GieIange these definitions, if necessary, to match a particular platform
#if defined(_WIN16) || defined(__MSDOS__) || defined(_MSDOS)
   // 16 bit systems use long int for 32 bit integer
   typedef long int           int32;   // 32 bit signed integer
   typedef unsigned long int  uint32;  // 32 bit unsigned integer
#else
   // Most other systems use int for 32 bit integer
   typedef int                int32;   // 32 bit signed integer
   typedef unsigned int       uint32;  // 32 bit unsigned integer
#endif

// Define 64 bit signed and unsigned integers, if possible
#if (defined(__WINDOWS__) || defined(_WIN32)) && (defined(_MSC_VER) || defined(__INTEL_COMPILER))
   // Microsoft and other compilers under Windows use __int64
   typedef __int64            int64;   // 64 bit signed integer
   typedef unsigned __int64   uint64;  // 64 bit unsigned integer
   #define INT64_DEFINED               // Remember that int64 is defined
#elif defined(__unix__) && (defined(_M_IX86) || defined(_M_X64))
   // Gnu and other compilers under Linux etc. use long long
   typedef long long          int64;   // 64 bit signed integer
   typedef unsigned long long uint64;  // 64 bit unsigned integer
   #define INT64_DEFINED               // Remember that int64 is defined
#else
   // 64 bit integers not defined
   // You may include definitions for other platforms here
#endif


/***********************************************************************
System-specific user interface functions
***********************************************************************/

void EndOfProgram(void);               // System-specific exit code (userintf.cpp)

void FatalError(char * ErrorText);     // System-specific error reporting (userintf.cpp)


/***********************************************************************
Define random number generator classes
***********************************************************************/

class CRandomMersenne {                // Encapsulate random number generator
#if 0
   // Define constants for type MT11213A:
#define MERS_N   351
#define MERS_M   175
#define MERS_R   19
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   17
#define MERS_A   0xE4BD75F5
#define MERS_B   0x655E5280
#define MERS_C   0xFFD58000
#else
   // or constants for type MT19937:
#define MERS_N   624
#define MERS_M   397
#define MERS_R   31
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   18
#define MERS_A   0x9908B0DF
#define MERS_B   0x9D2C5680
#define MERS_C   0xEFC60000
#endif
public:
   CRandomMersenne(uint32 seed) {      // Constructor
      RandomInit(seed); LastInterval = 0;}
   void RandomInit(uint32 seed);       // Re-seedSetDirectory["~/Dropbox/YDriveBurkina/codesB/DoubleSex/Mosaic"];
   void RandomInitByArray(uint32 seeds[], int length); // Seed by more than 32 bits
   int IRandom (int min, int max);     // Output random integer
   int IRandomX(int min, int max);     // Output random integer, exact
   double Random();                    // Output random float
   uint32 BRandom();                   // Output random bits
private:
   void Init0(uint32 seed);            // Basic initialization procedure
   uint32 mt[MERS_N];                  // State vector
   int mti;                            // Index into mt
   uint32 LastInterval;                // Last interval length for IRandomX
   uint32 RLimit;                      // Rejection limit used by IRandomX
   enum TArch {LITTLE_ENDIAN1, BIG_ENDIAN1, NONIEEE}; // Definition of architecture
   TArch Architecture;                 // Conversion to float depends on architecture
};


class CRandomMother {             // Encapsulate random number generator
public:
   void RandomInit(uint32 seed);       // Initialization
   int IRandom(int min, int max);      // Get integer random number in desired interval
   double Random();                    // Get floating point random number
   uint32 BRandom();                   // Output random bits
   CRandomMother(uint32 seed) {   // Constructor
      RandomInit(seed);}
protected:
   uint32 x[5];                        // History buffer
};

#endif

// body file: mersenne.cpp

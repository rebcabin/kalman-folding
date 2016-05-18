#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#include <vecLib/clapack.h>
typedef __CLPK_doublereal doublereal;
typedef __CLPK_integer integer;
#else
typedef double doublereal;
typedef int integer;
// or typedef long int integer; on some systems
extern "C" void F77NAME(dspevd)
   (char*, char*, integer*, doublereal*, doublereal*,
    doublereal*, integer*, doublereal*, integer*,
    integer*, integer*, integer*);
#endif

integer eigensolve(vector<doublereal> &H, 
                   vector<doublereal> &Eval, 
                   vector<doublereal> &Evec)
{
   // Solve the eigenvalue problem with LAPACK's dsepvd routine

   integer N = Eval.size();
   assert(H.size() == N*(N+1)/2);
   assert(Evec.size() == N*N);

   integer info;
   char jobz='V';
   char uplo='U';
   vector<doublereal> work(1+6*N+N*N);
   integer lwork = work.size();
   vector<integer> iwork(3+5*N);
   integer liwork = iwork.size();

   F77NAME(dspevd)(&jobz,&uplo,&N,&(H[0]),&(Eval[0]),&(Evec[0]),&N,
           &(work[0]),&lwork,&(iwork[0]),&liwork,&info);

   return info;
}

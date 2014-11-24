/**************************************************************************
  Module:  Stress
  Purpose: Utility routines for the solid stress calculations
 ***************************************************************************/
#include <float.h>
#include <math.h>
#define MIN(a, b) ( (a < b) ? a : b)
#define MAX(a, b) ( (a > b) ? a : b)
#include "omp.h" 
typedef double real8;
#if 1
/************************************************************************
 * Function  : StressJaumRot
 * Purpose   : Increment the stress by an amount corresponding to
 *             the rotation terms in the Jaumann stress rate.
 ************************************************************************/

void StressJaumRot(real8 *newSxx,real8 *newSyy,real8 *newSzz,real8 *newTxy,real8 *newTxz,real8 *newTyz,const real8 *sxx,const real8 *syy,const real8 *txy,const real8 *txz,const real8 *tyz,const real8 *wxx,const real8 *wyy,const real8 *wzz,const int *zoneset,real8 deltaTime,int length)
{
  int index;
  int i;
  
#pragma omp parallel for private (index,i) firstprivate (deltaTime,length)
  for (i = 0; i <= length - 1; i += 1) {
    index = zoneset[i];
    newSxx[i] = sxx[zoneset[i]] + deltaTime * (- 2. * txy[zoneset[i]] * wzz[zoneset[i]] + 2. * txz[zoneset[i]] * wyy[zoneset[i]]);
    newSyy[i] = syy[zoneset[i]] + deltaTime * (2. * txy[zoneset[i]] * wzz[zoneset[i]] - 2. * tyz[zoneset[i]] * wxx[zoneset[i]]);
    newSzz[i] = -sxx[zoneset[i]] - syy[zoneset[i]] + deltaTime * (- 2. * txz[zoneset[i]] * wyy[zoneset[i]] + 2. * tyz[zoneset[i]] * wxx[zoneset[i]]);
    newTxy[i] = txy[zoneset[i]] + deltaTime * (wzz[zoneset[i]] * (sxx[zoneset[i]] - syy[zoneset[i]]) + wyy[zoneset[i]] * tyz[zoneset[i]] - wxx[zoneset[i]] * txz[zoneset[i]]);
    newTyz[i] = tyz[zoneset[i]] + deltaTime * (wxx[zoneset[i]] * (2. * syy[zoneset[i]] + sxx[zoneset[i]]) + wzz[zoneset[i]] * txz[zoneset[i]] - wyy[zoneset[i]] * txy[zoneset[i]]);
    newTxz[i] = txz[zoneset[i]] + deltaTime * (wyy[zoneset[i]] * (-syy[zoneset[i]] - 2. * sxx[zoneset[i]]) + wxx[zoneset[i]] * txy[zoneset[i]] - wzz[zoneset[i]] * tyz[zoneset[i]]);
  }
}
/************************************************************************
 * Function  : StressIncrIsoElas
 * Purpose   : Increment the deviatoric stresses assuming
 *             isotropic elastic behavior.  Also update the
 *             corresponding value of  sqrt(3/2) * J2.
 ************************************************************************/

void StressIncrIsoElas(real8 *newSxx,real8 *newSyy,real8 *newSzz,real8 *newTxy,real8 *newTxz,real8 *newTyz,real8 *fun2j,const real8 *dxx,const real8 *dyy,const real8 *dzz,const real8 *dxy,const real8 *dxz,const real8 *dyz,const real8 *shearMod,const int *zoneset,real8 deltaTime,int length)
{
  int i;
  int index;
  real8 twoj;
  real8 twoDelta = 2. * deltaTime;
  
#pragma omp parallel for private (index,twoj,i) firstprivate (length,twoDelta)
  for (i = 0; i <= length - 1; i += 1) {
    index = zoneset[i];
    newSxx[i] = newSxx[i] + twoDelta * shearMod[zoneset[i]] * dxx[zoneset[i]];
    newSyy[i] = newSyy[i] + twoDelta * shearMod[zoneset[i]] * dyy[zoneset[i]];
    newSzz[i] = newSzz[i] + twoDelta * shearMod[zoneset[i]] * dzz[zoneset[i]];
    newTxy[i] = newTxy[i] + twoDelta * shearMod[zoneset[i]] * dxy[zoneset[i]];
    newTyz[i] = newTyz[i] + twoDelta * shearMod[zoneset[i]] * dyz[zoneset[i]];
    newTxz[i] = newTxz[i] + twoDelta * shearMod[zoneset[i]] * dxz[zoneset[i]];
    twoj = newSxx[i] * newSxx[i] + newSyy[i] * newSyy[i] + newSzz[i] * newSzz[i] + 2. * (newTxy[i] * newTxy[i] + newTxz[i] * newTxz[i] + newTyz[i] * newTyz[i]);
    fun2j[i] = sqrt(1.5 * twoj);
  }
}
/************************************************************************
 * Function  : StressScale
 * Purpose   : 
 ************************************************************************/

void StressScale(real8 *Sxx,real8 *Syy,real8 *Szz,real8 *Txy,real8 *Txz,real8 *Tyz,real8 *fun2j,const real8 *shearMod,const real8 *shearRatio,const int *zoneset,int length,int Softening)
{
  int i;
  real8 twoj;
  if (Softening == 0) {
/* Do Nothing */
  }
  else {
    
#pragma omp parallel for private (twoj,i) firstprivate (length)
    for (i = 0; i <= length - 1; i += 1) {
      int index = zoneset[i];
      real8 rat = shearRatio[i] * shearMod[zoneset[i]];
      Sxx[i] = Sxx[i] * rat;
      Syy[i] = Syy[i] * rat;
      Szz[i] = Szz[i] * rat;
      Txy[i] = Txy[i] * rat;
      Tyz[i] = Tyz[i] * rat;
      Txz[i] = Txz[i] * rat;
      twoj = Sxx[i] * Sxx[i] + Syy[i] * Syy[i] + Szz[i] * Szz[i] + 2. * (Txy[i] * Txy[i] + Txz[i] * Txz[i] + Tyz[i] * Tyz[i]);
      fun2j[i] = sqrt(1.5 * twoj);
    }
  }
}
/************************************************************************
 * Function  : StressCheckYieldJ2
 * Purpose   : Correct stress with radial return if current (trial) stress
 *             lies outside the current (mises) yield surface.
 ************************************************************************/

void StressCheckYieldJ2(real8 *newSxx,real8 *newSyy,real8 *newSzz,real8 *newTxy,real8 *newTxz,real8 *newTyz,const real8 *yield,const real8 *fun2j,const int *zoneset,int length)
{
  int i;
  int index;
  real8 scale;
  
#pragma omp parallel for private (index,scale,i) firstprivate (length)
  for (i = 0; i <= length - 1; i += 1) {
    index = zoneset[i];
    if (fun2j[i] > yield[zoneset[i]]) {
      scale = yield[zoneset[i]] / fun2j[i];
      newSxx[i] = newSxx[i] * scale;
      newSyy[i] = newSyy[i] * scale;
      newSzz[i] = newSzz[i] * scale;
      newTxy[i] = newTxy[i] * scale;
      newTxz[i] = newTxz[i] * scale;
      newTyz[i] = newTyz[i] * scale;
    }
  }
}
/************************************************************************
 * Function  : StressCheckEpsFail
 * Purpose   : 
 ************************************************************************/

void StressCheckEpsFail(real8 *newSxx,real8 *newSyy,real8 *newSzz,real8 *newTxy,real8 *newTxz,real8 *newTyz,real8 *eps,real8 eps_failure_model,const int *zoneset,int length)
{
  int i;
  int index;
  
#pragma omp parallel for private (index,i) firstprivate (eps_failure_model,length)
  for (i = 0; i <= length - 1; i += 1) {
    index = zoneset[i];
    if (eps[zoneset[i]] > eps_failure_model) {
      newSxx[i] = 0.0;
      newSyy[i] = 0.0;
      newSzz[i] = 0.0;
      newTxy[i] = 0.0;
      newTxz[i] = 0.0;
      newTyz[i] = 0.0;
      eps[zoneset[i]] = eps_failure_model * 1.01;
    }
  }
}
#endif 
/************************************************************************
 * Function  : StressZero
 * 
 * Purpose   : 
 ************************************************************************/

void StressZero(real8 *newSxx,real8 *newSyy,real8 *newSzz,real8 *newTxy,real8 *newTxz,real8 *newTyz,const real8 *fun2j,const real8 *shearMod,real8 eosvmax,real8 stresscut,const int *zoneset,const real8 *vc,int length)
{
  int i;
  int index;
/* This value 1.e-20 is used to prevent underflow. It is NOT a
     cuttoff. DO NOT TOUCH THIS VALE. */
  real8 stress2 = stresscut * 1.e-20;
  real8 nstres2 = -stress2;
  
#pragma omp parallel for private (index,i) firstprivate (length,stress2)
  for (i = 0; i <= length - 1; i += 1) {
    index = zoneset[i];
    if (shearMod[zoneset[i]] == 0.0 || fun2j[i] < stresscut || vc[i] >= eosvmax) {
      newSxx[i] = 0.0;
      newSyy[i] = 0.0;
      newSzz[i] = 0.0;
      newTxy[i] = 0.0;
      newTxz[i] = 0.0;
      newTyz[i] = 0.0;
    }
#if 1
    if (newSxx[i] < stress2 && newSxx[i] > nstres2) {
      newSxx[i] = 0.0;
    }
    if (newSyy[i] < stress2 && newSyy[i] > nstres2) {
      newSyy[i] = 0.0;
    }
    if (newSzz[i] < stress2 && newSzz[i] > nstres2) {
      newSzz[i] = 0.0;
    }
    if (newTxy[i] < stress2 && newTxy[i] > nstres2) {
      newTxy[i] = 0.0;
    }
    if (newTxz[i] < stress2 && newTxz[i] > nstres2) {
      newTxz[i] = 0.0;
    }
    if (newTyz[i] < stress2 && newTyz[i] > nstres2) {
      newTyz[i] = 0.0;
    }
#endif
  }
}
#if 1
/************************************************************************
 * Function  : StressStrainWork
 * Purpose   : 
 ************************************************************************/

void StressStrainWork(real8 *deltz,real8 *delts,const real8 *newSxx,const real8 *newSyy,const real8 *newSzz,const real8 *newTxy,const real8 *newTxz,const real8 *newTyz,const real8 *sxx,const real8 *syy,const real8 *txy,const real8 *txz,const real8 *tyz,const real8 *dxx,const real8 *dyy,const real8 *dzz,const real8 *dxy,const real8 *dxz,const real8 *dyz,real8 deltaTime,const int *zoneset,const real8 *vc,const real8 *vnewc,int length)
{
  int i;
  int index;
  real8 quarterDelta = 0.25 * deltaTime;
  real8 szz;
  
#pragma omp parallel for private (index,szz,i) firstprivate (length,quarterDelta)
  for (i = 0; i <= length - 1; i += 1) {
    index = zoneset[i];
    szz = -sxx[zoneset[i]] - syy[zoneset[i]];
    deltz[zoneset[i]] += quarterDelta * (vnewc[i] + vc[i]) * (dxx[zoneset[i]] * (sxx[zoneset[i]] + newSxx[i]) + dyy[zoneset[i]] * (syy[zoneset[i]] + newSyy[i]) + dzz[zoneset[i]] * (szz + newSzz[i]) + 2. * dxy[zoneset[i]] * (txy[zoneset[i]] + newTxy[i]) + 2. * dxz[zoneset[i]] * (txz[zoneset[i]] + newTxz[i]) + 2. * dyz[zoneset[i]] * (tyz[zoneset[i]] + newTyz[i]));
    delts[i] += quarterDelta * (vnewc[i] + vc[i]) * (dxx[zoneset[i]] * sxx[zoneset[i]] + dyy[zoneset[i]] * syy[zoneset[i]] + dzz[zoneset[i]] * szz + 2. * dxy[zoneset[i]] * txy[zoneset[i]] + 2. * dxz[zoneset[i]] * txz[zoneset[i]] + 2. * dyz[zoneset[i]] * tyz[zoneset[i]]);
  }
}
/************************************************************************
 * Function  : StressFailureModel
 * Purpose   : 
 ************************************************************************/

void StressFailureModel(real8 *newSxx,real8 *newSyy,real8 *newSzz,real8 *newTxy,real8 *newTxz,real8 *newTyz,real8 *eps,real8 epsFailModel,const int *zoneset,int length)
{
  if (epsFailModel > 0.0) {
    StressCheckEpsFail(newSxx,newSyy,newSzz,newTxy,newTxz,newTyz,eps,epsFailModel,zoneset,length);
  }
}
/************************************************************************
 * Function  : StressCalcShearRatio
 * Purpose   : 
 ************************************************************************/

void StressCalcShearRatio(const real8 *shearMod,real8 *shearRatio,const int *zoneset,int length)
{
  int i;
  int index;
  
#pragma omp parallel for private (index,i) firstprivate (length)
  for (i = 0; i <= length - 1; i += 1) {
    index = zoneset[i];
    if (shearMod[zoneset[i]] > 0) {
      shearRatio[i] = 1. / shearMod[zoneset[i]];
    }
    else {
      shearRatio[i] = 0.0;
    }
  }
}
/************************************************************************
 * Function  : StressStrainHeat
 * Purpose   : 
 ************************************************************************/

void StressStrainHeat(const real8 *deltz,real8 *deltzh,real8 *deltrh,const real8 *shearMod,const real8 *shearRatio,const real8 *shearDer,const real8 *newSxx,const real8 *newSyy,const real8 *newSzz,const real8 *newTxy,const real8 *newTxz,const real8 *newTyz,const real8 *sxx,const real8 *syy,const real8 *txy,const real8 *txz,const real8 *tyz,real8 deltaTime,const int *zoneset,const real8 *vc,const real8 *vnewc,int length)
{
  real8 shearr;
  real8 sheari;
  real8 avgMod;
  int nz;
  int i;
/* Quiet the compiler - unused argument */
  deltaTime = deltaTime;
  
#pragma omp parallel for private (shearr,sheari,avgMod,nz,i) firstprivate (length)
  for (i = 0; i <= length - 1; i += 1) {
    nz = zoneset[i];
    shearr = 0.5 * shearRatio[i];
    if (shearMod[zoneset[i]] > 0.0) {
      sheari = 0.5 / shearMod[zoneset[i]];
      deltrh[zoneset[i]] = 0.25 * (vnewc[i] + vc[i]) * ((newSxx[i] * sheari - sxx[zoneset[i]] * shearr) * (sxx[zoneset[i]] + newSxx[i]) + (newSyy[i] * sheari - syy[zoneset[i]] * shearr) * (syy[zoneset[i]] + newSyy[i]) + (newSzz[i] * sheari + (syy[zoneset[i]] + sxx[zoneset[i]]) * shearr) * (newSzz[i] - sxx[zoneset[i]] - syy[zoneset[i]]) + 2. * (newTxy[i] * sheari - txy[zoneset[i]] * shearr) * (txy[zoneset[i]] + newTxy[i]) + 2. * (newTxz[i] * sheari - txz[zoneset[i]] * shearr) * (txz[zoneset[i]] + newTxz[i]) + 2. * (newTyz[i] * sheari - tyz[zoneset[i]] * shearr) * (tyz[zoneset[i]] + newTyz[i]));
    }
    else {
      deltrh[zoneset[i]] = - 0.25 * (vnewc[i] + vc[i]) * (sxx[zoneset[i]] * (sxx[zoneset[i]] + newSxx[i]) + syy[zoneset[i]] * (syy[zoneset[i]] + newSyy[i]) - (syy[zoneset[i]] + sxx[zoneset[i]]) * (newSzz[i] - sxx[zoneset[i]] - syy[zoneset[i]]) + 2. * txy[zoneset[i]] * (txy[zoneset[i]] + newTxy[i]) + 2. * txz[zoneset[i]] * (txz[zoneset[i]] + newTxz[i]) + 2. * tyz[zoneset[i]] * (tyz[zoneset[i]] + newTyz[i])) * shearr;
    }
    deltzh[zoneset[i]] = deltz[zoneset[i]] - deltrh[zoneset[i]];
    avgMod = 0.5 * shearMod[zoneset[i]];
    if (shearRatio[i] > 0.0) {
      avgMod = avgMod + 0.5 / shearRatio[i];
    }
    if (avgMod > 0.0) {
      deltrh[zoneset[i]] = shearDer[i] * deltrh[zoneset[i]] / avgMod;
    }
    else {
      deltrh[zoneset[i]] = 0.0;
    }
  }
}
//************************************************************************
///
/// @note The principal components of the total stress tensor are returned.
/// @param pstress Principal stresses, from largest to smallest
/// @param sx xx-component of deviatoric stress tensor
/// @param sy yy-component of deviatoric stress tensor
/// @param txy xy-component of deviatoric stress tensor
/// @param txz xz-component of deviatoric stress tensor
/// @param tyz yz-component of deviatoric stress tensor
/// @param p Pressure
/// @param acut Magnitude of cutoff for a
/// @pre  none
/// @post none
///

void StressCalculatePrincipalValues(real8 pstress[3],const real8 sx,const real8 sy,const real8 txy,const real8 txz,const real8 tyz,const real8 p,const real8 acut)
{
  real8 third = 1. / 3.0;
  real8 sz = -(sx + sy);
/* Determine principle values */
  real8 txy2 = txy * txy;
  real8 txz2 = txz * txz;
  real8 tyz2 = tyz * tyz;
  real8 a = (sx * sy + sy * sz + sz * sx - (txy2 + txz2 + tyz2)) / 3.0;
/* If a is at all positive here, it will blow up the sqrt
       * below... originally hardwired to -1e-13 */
  a = (a < -acut?a : -acut);
  real8 det = sx * sy * sz + 2. * txy * txz * tyz - sx * tyz2 - sy * txz2 - sz * txy2;
  real8 r = sqrt(-pow(a,3.0));
  real8 arg = det / 2. / r;
/* Cap the value of arg to (-1 <= arg <= 1) */
  arg = (arg > - 1.?arg : - 1.);
  arg = (arg < 1.?arg : 1.);
  real8 theta = acos(arg);
  real8 coef = 2. * pow(r,third);
/* Calculate eigenvalues of 3x3 matrix */
  real8 p1 = coef * cos(theta / 3.0);
  real8 p2 = coef * cos((theta + 2. * 3.14159265358979323846) / 3.0);
  real8 p3 = coef * cos((theta + 4.0 * 3.14159265358979323846) / 3.0);
  if (p2 > p1) {
    real8 tem = p1;
    p1 = p2;
    p2 = tem;
  }
  if (p3 > p1) {
    real8 tem = p1;
    p1 = p3;
    p3 = tem;
  }
  if (p3 > p2) {
    real8 tem = p2;
    p2 = p3;
    p3 = tem;
  }
/* Return principal components of total stress */
  pstress[0] = p1 - p;
  pstress[1] = p2 - p;
  pstress[2] = p3 - p;
}
#endif 

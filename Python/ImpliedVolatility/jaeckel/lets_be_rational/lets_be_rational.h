//
// This source code resides at www.jaeckel.org/LetsBeRational.7z .
//
// ======================================================================================
// Copyright (c) 2013-2023 Peter Jaeckel.
//
// Permission to use, copy, modify, and distribute this software is freely granted,
// provided that this notice is preserved.
//
// WARRANTY DISCLAIMER
// The Software is provided "as is" without warranty of any kind, either express or implied,
// including without limitation any implied warranties of condition, uninterrupted use,
// merchantability, fitness for a particular purpose, or non-infringement.
// ======================================================================================
//
#ifndef   LETS_BE_RATIONAL_H
#define   LETS_BE_RATIONAL_H

#include "dllmain.h"

// As defined in equation (1.1) in "Let's Be Rational".
extern "C" DLL_EXPORT double Black(double F, double K, double sigma, double T, double q /* q=+/-1 */);
// As defined in equation (2.3) in "Let's Be Rational": beta(x,s) := B(F,K,sigma,T,theta=+/-1)/sqrt(F*K) with x=ln(F/K) and s=sigmasqrtT
extern "C" DLL_EXPORT double NormalisedBlack(double x, double s, double q /* q=+/-1 */);

extern "C" DLL_EXPORT double ImpliedBlackVolatility(double price, double F, double K, double T, double q /* q=+/-1 */);

// The input 'beta' is taken as a "normalised Black" price as defined in equation (2.3) in "Let's Be Rational".
extern "C" DLL_EXPORT double NormalisedImpliedBlackVolatility(double beta, double x, double q /* q=+/-1 */);

//    bbar(x,s,theta)          :=   b_m_a_x(x,theta)   -  b(x,s,theta)
//                       =   e_x/_sqr*Phi(-x/s-s/2) + e_minus_x/_sqr*Phi(x/s-s/2)                                            |     for both theta = +/-1
// Same for calls and puts, i.e., no dependency on theta = +/-1.
extern "C" DLL_EXPORT double ComplementaryNormalisedBlack(double x, double s);

#define EXPOSE_VEGA_AND_MORE

#if defined(EXPOSE_VEGA_AND_MORE) || !defined(SWIG)
// dBlack(F,K,sigma,T)/dsigma  [no dependency on the call/put flag theta=+/-1]
extern "C" DLL_EXPORT double Vega(double F, double K, double sigma, double T);
// dbeta(x,s)/ds          [no dependency on the call/put flag theta=+/-1]   with x=ln(F/K) and s=sigmasqrtT
extern "C" DLL_EXPORT double NormalisedVega(double x, double s);
// d_sqrBlack(F,K,sigma,T)/dsigma_sqr  [no dependency on the call/put flag theta=+/-1]
extern "C" DLL_EXPORT double Volga(double F, double K, double sigma, double T);
// d_sqrbeta(x,s)/ds_sqr          [no dependency on the call/put flag theta=+/-1]   with x=ln(F/K) and s=sigmasqrtT
extern "C" DLL_EXPORT double NormalisedVolga(double x, double s);
// The attainable *relative* accuracy of beta = b(s) when s has *relative* accuracy epsilon is (to lowest order) (|s*b'(s)/b(x)|+1)*epsilon --- see the source code for a detailed derivation.
// The attainable *relative* accuracy of x = b_minus_1(beta) when beta has *relative* accuracy epsilon is (to lowest order) (|b(s)/(s*b'(s))|+1)*epsilon .
// This function returns (s*db(x,s)/ds)/b(x,s,theta=+/-1). In order to get the accuracy limit of implied volatility calculations, take (1+1/BlackAccuracyFactor(x,s,theta))*DBL_EPSILON.
extern "C" DLL_EXPORT double BlackAccuracyFactor(double x /* = ln(F/K) */, double s /* = sigmasqrtT */, double q /* q=+/-1 */);
// The attainable *relative* accuracy of x = b_minus_1(beta) when beta has *relative* accuracy epsilon is (to lowest order) (|b(s)/(s*b'(s))|+1)*epsilon .
extern "C" DLL_EXPORT double ImpliedVolatilityAttainableAccuracy(double x, double s, double q /* q=+/-1 */);

// DBL_EPSILON
extern "C" DLL_EXPORT double DblEpsilon();
// DBL_MIN
extern "C" DLL_EXPORT double DblMin();
// DBL_MAX
extern "C" DLL_EXPORT double DblMax();
#endif

//#define ENABLE_CHANGING_THE_MAXIMUM_ITERATION_COUNT
#ifdef ENABLE_CHANGING_THE_MAXIMUM_ITERATION_COUNT
extern "C" DLL_EXPORT int set_implied_volatility_maximum_iterations(int n);
#endif

//#define ENABLE_CHANGING_THE_HOUSEHOLDER_METHOD_ORDER
#ifdef ENABLE_CHANGING_THE_HOUSEHOLDER_METHOD_ORDER
extern "C" DLL_EXPORT int set_implied_volatility_householder_method_order(int m);
#endif

//#define ENABLE_SWITCHING_THE_OUTPUT_TO_ITERATION_COUNT
#ifdef ENABLE_SWITCHING_THE_OUTPUT_TO_ITERATION_COUNT
extern "C" DLL_EXPORT int set_implied_volatility_output_type(int type);
#endif

#endif // LETS_BE_RATIONAL_H

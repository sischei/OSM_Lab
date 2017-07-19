/*
 * Code Author: Miroslav Stoyanov
 *
 * Copyright (C) 2016  Miroslav Stoyanov
 *
 * This file is part of
 * Toolkit for Adaptive Stochastic Modeling And Non-Intrusive Approximation
 *              a.k.a. TASMANIAN
 *
 * TASMANIAN is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * TASMANIAN is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with TASMANIAN.  If not, see <http://www.gnu.org/licenses/>
 *
 */

#ifndef __TASMANIAN_SPARSE_GRID_H
#define __TASMANIAN_SPARSE_GRID_H

// ------------ C Interface for TasmanianSparseGrid -------------- //
// NOTE: you need to explicitly call the constructor and destructor
//       void * grid is a pointer to a C++ class
void* tsgConstructTasmanianSparseGrid();
void tsgDestructTasmanianSparseGrid( void * grid );

void tsgCopyGrid( void * destination, void * source );

const char* tsgGetVersion();
const char* tsgGetLicense();

void tsgWrite( void * grid, const char* filename );
int tsgRead( void * grid, const char* filename );

void tsgMakeGlobalGrid( void * grid, int dimensions, int outputs, int depth, const char * sType, const char * sRule, int * anisotropic_weights, double alpha, double beta, const char* custom_filename );
void tsgMakeSequenceGrid( void * grid, int dimensions, int outputs, int depth, const char * sType, const char * sRule, int * anisotropic_weights );
void tsgMakeLocalPolynomialGrid( void * grid, int dimensions, int outputs, int depth, int order, const char * sRule );
void tsgMakeWaveletGrid( void * grid, int dimensions, int outputs, int depth, int order );

void tsgUpdateGlobalGrid( void * grid, int depth, const char * sType, const int *anisotropic_weights );
void tsgUpdateSequenceGrid( void * grid, int depth, const char * sType, const int *anisotropic_weights );

double tsgGetAlpha( void * grid );
double tsgGetBeta( void * grid );
int tsgGetOrder( void * grid );
int tsgGetNumDimensions( void * grid );
int tsgGetNumOutputs( void * grid );
const char * tsgGetRule( void * grid );
const char * tsgGetCustomRuleDescription( void * grid );

int tsgGetNumLoaded( void * grid );
int tsgGetNumNeeded( void * grid );
int tsgGetNumPoints( void * grid );

double* tsgGetLoadedPoints( void * grid );
double* tsgGetNeededPoints( void * grid );
double* tsgGetPoints( void * grid );

double* tsgGetQuadratureWeights( void * grid );
double* tsgGetInterpolationWeights( void * grid, const double *x );

void tsgLoadNeededPoints( void * grid, const double *vals );

void tsgEvaluate( void * grid, const double *x, double *y );
void tsgIntegrate( void * grid, double *q );

// add an batch-eval function using OpenMP
void tsgBatchEvaluate( void * grid, const double *x, int num_x, double *y );
double* tsgBatchGetInterpolationWeights( void * grid, const double *x, int num_x );

int tsgIsGlobal( void * grid );
int tsgIsSequence( void * grid );
int tsgIsLocalPolynomial( void * grid );
int tsgIsWavelet( void * grid );

void tsgSetDomainTransform( void * grid, const double a[], const double b[] );
bool tsgIsSetDomainTransfrom( void * grid );
void tsgClearDomainTransform( void * grid );
void tsgGetDomainTransform( void * grid, double a[], double b[] );

void tsgSetAnisotropicRefinement( void * grid, const char * sType, int min_growth, int output );
int* tsgEstimateAnisotropicCoefficients( void * grid, const char * sType, int output, int *num_coefficients );
void tsgSetGlobalSurplusRefinement( void * grid, double tolerance, int output );
void tsgSetLocalSurplusRefinement( void * grid, double tolerance, const char * sRefinementType, int output );
void tsgClearRefinement( void * grid );

int* tsgGetGlobalPolynomialSpace( void * grid, int interpolation, int *num_indexes );

void tsgPrintStats( void *grid );

#endif

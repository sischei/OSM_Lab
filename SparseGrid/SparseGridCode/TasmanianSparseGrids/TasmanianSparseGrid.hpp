/*
 * Copyright (c) 2017, Miroslav Stoyanov
 *
 * This file is part of
 * Toolkit for Adaptive Stochastic Modeling And Non-Intrusive ApproximatioN: TASMANIAN
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
 *    and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
 *    or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * UT-BATTELLE, LLC AND THE UNITED STATES GOVERNMENT MAKE NO REPRESENTATIONS AND DISCLAIM ALL WARRANTIES, BOTH EXPRESSED AND IMPLIED.
 * THERE ARE NO EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT,
 * COPYRIGHT, TRADEMARK, OR OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL ACCOMPLISH THE INTENDED RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE.
 * THE USER ASSUMES RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS, CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING FROM OR ARISING OUT OF,
 * IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE SOFTWARE.
 */

#ifndef __TASMANIAN_SPARSE_GRID_HPP
#define __TASMANIAN_SPARSE_GRID_HPP

#include "tsgEnumerates.hpp"

#include "tsgGridGlobal.hpp"
#include "tsgGridSequence.hpp"
#include "tsgGridLocalPolynomial.hpp"
#include "tsgGridWavelet.hpp"

#include <iomanip> // only needed for printStats()

namespace TasGrid{

class TasmanianSparseGrid{
public:
        TasmanianSparseGrid();
        TasmanianSparseGrid( const TasmanianSparseGrid &source );
        ~TasmanianSparseGrid();

        static const char* getVersion();
        static const char* getLicense();

        void write( const char *filename ) const;
        bool read( const char *filename );

        void write( std::ofstream &ofs ) const;
        bool read( std::ifstream &ifs );

        void makeGlobalGrid( int dimensions, int outputs, int depth, TypeDepth type, TypeOneDRule rule, const int *anisotropic_weights = 0, double alpha = 0.0, double beta = 0.0, const char* custom_filename = 0 );
        void makeSequenceGrid( int dimensions, int outputs, int depth, TypeDepth type, TypeOneDRule rule, const int *anisotropic_weights = 0 );
        void makeLocalPolynomialGrid( int dimensions, int outputs, int depth, int order = 1, TypeOneDRule rule = rule_localp );
        void makeWaveletGrid( int dimensions, int outputs, int depth, int order = 1 );
        void copyGrid( const TasmanianSparseGrid *source );

        void updateGlobalGrid( int depth, TypeDepth type, const int *anisotropic_weights = 0 );
        void updateSequenceGrid( int depth, TypeDepth type, const int *anisotropic_weights = 0 );

        double getAlpha() const;
        double getBeta() const;
        int getOrder() const;

        int getNumDimensions() const;
        int getNumOutputs() const;
        TypeOneDRule getRule() const;
        const char* getCustomRuleDescription() const; // used only for Global Grids with rule_customtabulated

        int getNumLoaded() const;
        int getNumNeeded() const;
        int getNumPoints() const; // returns the number of loaded points unless no points are loaded, then returns the number of needed points

        double* getLoadedPoints() const;
        double* getNeededPoints() const;
        double* getPoints() const; // returns the loaded points unless no points are loaded, then returns the needed points

        double* getQuadratureWeights() const;
        double* getInterpolationWeights( const double x[] ) const;

        void loadNeededPoints( const double *vals );

        void evaluate( const double x[], double y[] ) const;
        void integrate( double q[] ) const;

        bool isGlobal() const;
        bool isSequence() const;
        bool isLocalPolynomial() const;
        bool isWavelet() const;

        void setDomainTransform( const double a[], const double b[] ); // set the ranges of the box
        bool isSetDomainTransfrom() const;
        void clearDomainTransform();
        void getDomainTransform( double a[], double b[] ) const;

        void setAnisotropicRefinement( TypeDepth type, int min_growth, int output );
        int* estimateAnisotropicCoefficients( TypeDepth type, int output );
        void setSurplusRefinement( double tolerance, int output );
        void setSurplusRefinement( double tolerance, TypeRefinement criteria, int output = -1 ); // -1 indicates using all outputs
        void clearRefinement();

        int* getGlobalPolynomialSpace( bool interpolation, int &num_indexes ) const;

        void printStats() const;

        // WARNING: the functions below are mostly for debugging and research purposes
        //          modifying the returned pointers will result in undefined behavior
        void removePointsBySurplus( double tolerance, int output = -1 );

        double* evalHierarchicalFunctions( const double x[] ) const;
        void setHierarchicalCoefficients( const double c[] );

        const double* getSurpluses() const;
        const int* getPointsIndexes() const;
        const int* getNeededIndexes() const;

protected:
        void clear();

        void mapCanonicalToTransformed( int num_dimensions, int num_points, TypeOneDRule rule, double x[] ) const;
        void mapTransformedToCanonical( int num_dimensions, TypeOneDRule rule, double x[] ) const;
        double getQuadratureScale( int num_dimensions, TypeOneDRule rule ) const;

private:
        BaseCanonicalGrid *base;

        GridGlobal *global;
        GridSequence *sequence;
        GridLocalPolynomial *pwpoly;
        GridWavelet *wavelet;

        double *domain_transform_a, *domain_transform_b;
};


// ------------ C Interface for use with Python ctypes and potentially other C codes -------------- //
// ------------ if you need a C header, use the included TasmanianSparseGrids.h file -------------- //
/*extern "C" void* tsgConstructTasmanianSparseGrid();
extern "C" void tsgDestructTasmanianSparseGrid( void * grid );

extern "C" void tsgCopyGrid( void * destination, void * source );

extern "C" const char* tsgGetVersion();
extern "C" const char* tsgGetLicense();

extern "C" void tsgWrite( void * grid, const char* filename );
extern "C" int tsgRead( void * grid, const char* filename );

extern "C" void tsgMakeGlobalGrid( void * grid, int dimensions, int outputs, int depth, const char * sType, const char * sRule, int * anisotropic_weights, double alpha, double beta, const char* custom_filename );
extern "C" void tsgMakeSequenceGrid( void * grid, int dimensions, int outputs, int depth, const char * sType, const char * sRule, int * anisotropic_weights );
extern "C" void tsgMakeLocalPolynomialGrid( void * grid, int dimensions, int outputs, int depth, int order, const char * sRule );
extern "C" void tsgMakeWaveletGrid( void * grid, int dimensions, int outputs, int depth, int order );

extern "C" void tsgUpdateGlobalGrid( void * grid, int depth, const char * sType, const int *anisotropic_weights );
extern "C" void tsgUpdateSequenceGrid( void * grid, int depth, const char * sType, const int *anisotropic_weights );

extern "C" double tsgGetAlpha( void * grid );
extern "C" double tsgGetBeta( void * grid );
extern "C" int tsgGetOrder( void * grid );
extern "C" int tsgGetNumDimensions( void * grid );
extern "C" int tsgGetNumOutputs( void * grid );
extern "C" const char * tsgGetRule( void * grid );
extern "C" const char * tsgGetCustomRuleDescription( void * grid );

extern "C" int tsgGetNumLoaded( void * grid );
extern "C" int tsgGetNumNeeded( void * grid );
extern "C" int tsgGetNumPoints( void * grid );

extern "C" double* tsgGetLoadedPoints( void * grid );
extern "C" double* tsgGetNeededPoints( void * grid );
extern "C" double* tsgGetPoints( void * grid );

extern "C" double* tsgGetQuadratureWeights( void * grid );
extern "C" double* tsgGetInterpolationWeights( void * grid, const double *x );

extern "C" void tsgLoadNeededPoints( void * grid, const double *vals );

extern "C" void tsgEvaluate( void * grid, const double *x, double *y );
extern "C" void tsgIntegrate( void * grid, double *q );

// add an batch-eval function using OpenMP
extern "C" void tsgBatchEvaluate( void * grid, const double *x, int num_x, double *y );
extern "C" double* tsgBatchGetInterpolationWeights( void * grid, const double *x, int num_x );

extern "C" int tsgIsGlobal( void * grid );
extern "C" int tsgIsSequence( void * grid );
extern "C" int tsgIsLocalPolynomial( void * grid );
extern "C" int tsgIsWavelet( void * grid );

extern "C" void tsgSetDomainTransform( void * grid, const double a[], const double b[] );
extern "C" bool tsgIsSetDomainTransfrom( void * grid );
extern "C" void tsgClearDomainTransform( void * grid );
extern "C" void tsgGetDomainTransform( void * grid, double a[], double b[] );

extern "C" void tsgSetAnisotropicRefinement( void * grid, const char * sType, int min_growth, int output );
extern "C" int* tsgEstimateAnisotropicCoefficients( void * grid, const char * sType, int output, int *num_coefficients );
extern "C" void tsgSetGlobalSurplusRefinement( void * grid, double tolerance, int output );
extern "C" void tsgSetLocalSurplusRefinement( void * grid, double tolerance, const char * sRefinementType, int output );
extern "C" void tsgClearRefinement( void * grid );
extern "C" void tsgRemovePointsBySurplus( void * grid, double tolerance, int output );

extern "C" double* tsgEvalHierarchicalFunctions( void * grid, const double *x );
extern "C" double* tsgBatchEvalHierarchicalFunctions( void * grid, const double *x, int num_x );
extern "C" void tsgSetHierarchicalCoefficients( void * grid, const double *c );
extern "C" const double* tsgGetSurpluses( void *grid );

extern "C" int* tsgGetGlobalPolynomialSpace( void * grid, int interpolation, int *num_indexes );

extern "C" void tsgPrintStats( void *grid );

// for Python purposes only, those functions delete int and double pointers and are callable from Python ctypes
// the purpose is to avoid memory leaks in the Python interface, native C and C++ code should not call those functions
extern "C" void tsgDeleteDoubles( double * p );
extern "C" void tsgDeleteInts( int * p );*/

}

#endif

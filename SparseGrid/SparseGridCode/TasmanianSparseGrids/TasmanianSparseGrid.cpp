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

#ifndef __TASMANIAN_SPARSE_GRID_CPP
#define __TASMANIAN_SPARSE_GRID_CPP

#include "TasmanianSparseGrid.hpp"

namespace TasGrid{

const char* TasmanianSparseGrid::getVersion(){
        return "4.0";
}
const char* TasmanianSparseGrid::getLicense(){
        return "BSD 3-Clause";
}

TasmanianSparseGrid::TasmanianSparseGrid() : base(0), global(0), sequence(0), pwpoly(0), wavelet(0), domain_transform_a(0), domain_transform_b(0){}
TasmanianSparseGrid::TasmanianSparseGrid( const TasmanianSparseGrid &source ) : base(0), global(0), sequence(0), pwpoly(0), wavelet(0), domain_transform_a(0), domain_transform_b(0)
{  copyGrid( &source );  }
TasmanianSparseGrid::~TasmanianSparseGrid(){
        clear();
}

void TasmanianSparseGrid::clear(){
        if ( global != 0 ){ delete global; global = 0; }
        if ( sequence != 0 ){ delete sequence; sequence = 0; }
        if ( pwpoly != 0 ){ delete pwpoly; pwpoly = 0; }
        if ( wavelet != 0 ){ delete wavelet; wavelet = 0; }
        if ( domain_transform_a != 0 ){ delete[] domain_transform_a; domain_transform_a = 0; }
        if ( domain_transform_b != 0 ){ delete[] domain_transform_b; domain_transform_b = 0; }
        base = 0;
}

void TasmanianSparseGrid::write( const char *filename ) const{
        std::ofstream ofs; ofs.open( filename );
        write( ofs );
        ofs.close();
}
bool TasmanianSparseGrid::read( const char *filename ){
        std::ifstream ifs; ifs.open( filename );
        bool isGood = read( ifs );
        ifs.close();
        return isGood;
}

void TasmanianSparseGrid::write( std::ofstream &ofs ) const{
        ofs << "TASMANIAN SG " << getVersion() << endl;
        ofs << "WARNING: do not edit this manually" << endl;
        if ( global != 0 ){
                ofs << "global" << endl;
                global->write( ofs );
        }else if ( sequence != 0 ){
                ofs << "sequence" << endl;
                sequence->write( ofs );
        }else if ( pwpoly != 0 ){
                ofs << "localpolynomial" << endl;
                pwpoly->write( ofs );
        }else if ( wavelet != 0 ){
                ofs << "wavelet" << endl;
                wavelet->write( ofs );
        }else{
                ofs << "empty" << endl;
        }
        if ( domain_transform_a != 0 ){
                ofs << "custom" << endl;
                ofs << std::scientific; ofs.precision( 17 );
                for( int j=0; j<base->getNumDimensions(); j++ ){
                        ofs << domain_transform_a[j] << " " << domain_transform_b[j] << endl;
                }
        }else{
                ofs << "canonical" << endl;
        }
        ofs << "TASMANIAN SG end" << endl;
}
bool TasmanianSparseGrid::read( std::ifstream &ifs ){
        std::string T;
        ifs >> T;  if ( !(T.compare("TASMANIAN") == 0) ){  cerr << "ERROR: wrong file format, first word in not 'TASMANIAN'" << endl; return false; }
        ifs >> T;  if ( !(T.compare("SG") == 0) ){  cerr << "ERROR: wrong file format, second word in not 'SG'" << endl; return false; }
        getline( ifs, T ); T.erase(0,1); if ( !(T.compare( getVersion() ) == 0) ){  cerr << "WARNING: Version mismatch, possibly undefined behavior!" << endl; }
        getline( ifs, T ); if ( !(T.compare("WARNING: do not edit this manually") == 0) ){ cerr << "ERROR: wrong file format, did not match 'WARNING: do not edit this manually'" << endl; return false; }
        ifs >> T;
        if ( T.compare( "global" ) == 0 ){
                clear();
                global = new GridGlobal();
                global->read( ifs );
                base = global;
                getline( ifs, T );
        }else if ( T.compare( "sequence" ) == 0 ){
                clear();
                sequence = new GridSequence();
                sequence->read( ifs );
                base = sequence;
                getline( ifs, T );
        }else if ( T.compare( "localpolynomial" ) == 0 ){
                clear();
                pwpoly = new GridLocalPolynomial();
                pwpoly->read( ifs );
                base = pwpoly;
                getline( ifs, T );
        }else if ( T.compare( "wavelet" ) == 0 ){
                clear();
                wavelet = new GridWavelet();
                wavelet->read( ifs );
                base = wavelet;
                getline( ifs, T );
        }else if ( T.compare( "empty" ) == 0 ){
                clear();
        }else{
                cerr << "ERROR: wrong file format!" << endl; return false;
        }
        getline( ifs, T );
        bool read_last_line = true; // for version 3.1 compatibility
        if ( T.compare( "custom" ) == 0 ){
                domain_transform_a = new double[ base->getNumDimensions() ];
                domain_transform_b = new double[ base->getNumDimensions() ];
                for( int j=0; j<base->getNumDimensions(); j++ ){
                        ifs >> domain_transform_a[j] >> domain_transform_b[j];
                }
                getline( ifs, T );
        }else if ( T.compare( "canonical" ) == 0 ){
                // do nothing
        }else if ( T.compare("TASMANIAN SG end") == 0 ){
                // for compatibility with version 3.0 and the missing domain transform
                read_last_line = false;
        }else{
                cerr << "ERROR: wrong file format! Domain is not specified!" << endl; return false;
        }
        if ( read_last_line ){
                getline( ifs, T );
                if ( !(T.compare("TASMANIAN SG end") == 0) ){ cerr << "WARNING: file did not end with 'TASMANIAN SG end', this may result in undefined behavior" << endl; }
        }
        return true;
}

void TasmanianSparseGrid::makeGlobalGrid( int dimensions, int outputs, int depth, TypeDepth type, TypeOneDRule rule, const int *anisotropic_weights, double alpha, double beta, const char* custom_filename ){
        clear();
        global = new GridGlobal();
        global->makeGrid( dimensions, outputs, depth, type, rule, anisotropic_weights, alpha, beta, custom_filename );
        base = global;
}
void TasmanianSparseGrid::makeSequenceGrid( int dimensions, int outputs, int depth, TypeDepth type, TypeOneDRule rule, const int *anisotropic_weights ){
        if ( outputs < 1 ){
                cerr << "ERROR: makeSequenceGrid is called with zero outputs, for zero outputs use makeGlobalGrid instead" << endl;
                return;
        }
        if ( OneDimensionalMeta::isSequence( rule ) ){
                clear();
                sequence = new GridSequence();
                sequence->makeGrid( dimensions, outputs, depth, type, rule, anisotropic_weights );
                base = sequence;
        }else{
                cerr << "ERROR: makeSequenceGrid is called with rule " << OneDimensionalMeta::getIORuleString(rule) << " which is not a sequence rule" << endl;
        }
}
void TasmanianSparseGrid::makeLocalPolynomialGrid( int dimensions, int outputs, int depth, int order, TypeOneDRule rule ){
        if ( (rule != rule_localp) && (rule != rule_localp0) && (rule != rule_semilocalp) ){
                cout << "ERROR: makeLocalPolynomialGrid is called with rule " << OneDimensionalMeta::getIORuleString(rule) << " which is not a local polynomial rule" << endl;
                cout << "       use either " << OneDimensionalMeta::getIORuleString(rule_localp) << " or " << OneDimensionalMeta::getIORuleString(rule_semilocalp)
                                 << " or " << OneDimensionalMeta::getIORuleString(rule_localp0)  << endl;
                return;
        }
        if ( order < -1 ){
                cout << "ERROR: makeLocalPolynomialGrid is called with order " << order << ", but the order cannot be less than -1." << endl;
                return;
        }
        clear();
        pwpoly = new GridLocalPolynomial();
        pwpoly->makeGrid( dimensions, outputs, depth, order, rule );
        base = pwpoly;
}
void TasmanianSparseGrid::makeWaveletGrid( int dimensions, int outputs, int depth, int order ){
        if ( (order != 1) && (order != 3) ){
                cout << "ERROR: makeWaveletGrid is called with order " << order << ", but wavelets are implemented only for orders 1 and 3." << endl;
                return;
        }
        clear();
        wavelet = new GridWavelet();
        wavelet->makeGrid( dimensions, outputs, depth, order );
        base = wavelet;
}
void TasmanianSparseGrid::copyGrid( const TasmanianSparseGrid *source ){
        clear();
        if ( source->global != 0 ){
                global = new GridGlobal( *(source->global) );
                base = global;
        }else if ( source->sequence != 0 ){
                sequence = new GridSequence( *(source->sequence) );
                base = sequence;
        }else if ( source->pwpoly != 0 ){
                pwpoly = new GridLocalPolynomial( *(source->pwpoly) );
                base = pwpoly;
        }else if ( source->wavelet != 0 ){
                wavelet = new GridWavelet( *(source->wavelet) );
                base = wavelet;
        }
        if ( source->domain_transform_a != 0 ){
                setDomainTransform( source->domain_transform_a, source->domain_transform_b );
        }
}

void TasmanianSparseGrid::updateGlobalGrid( int depth, TypeDepth type, const int *anisotropic_weights ){
        if ( global != 0 ){
                global->updateGrid( depth, type, anisotropic_weights );
        }else{
                cerr << "ERROR: updateGlobalGrid called, but the grid is not global" << endl;
        }
}
void TasmanianSparseGrid::updateSequenceGrid( int depth, TypeDepth type, const int *anisotropic_weights ){
        if ( sequence != 0 ){
                sequence->updateGrid( depth, type, anisotropic_weights );
        }else{
                cerr << "ERROR: updateSequenceGrid called, but the grid is not sequence" << endl;
        }
}

double TasmanianSparseGrid::getAlpha() const{
        return ( global != 0 ) ? global->getAlpha() : 0.0;
}
double TasmanianSparseGrid::getBeta() const{
        return ( global != 0 ) ? global->getBeta() : 0.0;
}
int TasmanianSparseGrid::getOrder() const{
        return ( pwpoly != 0 ) ? pwpoly->getOrder() : ( ( wavelet != 0 ) ? wavelet->getOrder() : -1 );
}

int TasmanianSparseGrid::getNumDimensions() const{ return ( base == 0 ) ? 0 : base->getNumDimensions(); }
int TasmanianSparseGrid::getNumOutputs() const{ return ( base == 0 ) ? 0 : base->getNumOutputs(); }
TypeOneDRule TasmanianSparseGrid::getRule() const{ return ( base == 0 ) ? rule_none :base->getRule(); }
const char* TasmanianSparseGrid::getCustomRuleDescription() const{ return (global != 0) ? global->getCustomRuleDescription() : ""; }

int TasmanianSparseGrid::getNumLoaded() const{ return ( base == 0 ) ? 0 : base->getNumLoaded(); }
int TasmanianSparseGrid::getNumNeeded() const{ return ( base == 0 ) ? 0 : base->getNumNeeded(); }
int TasmanianSparseGrid::getNumPoints() const{ return ( base == 0 ) ? 0 : base->getNumPoints(); }

double* TasmanianSparseGrid::getLoadedPoints() const{
        if ( domain_transform_a == 0 ){
                return base->getLoadedPoints();
        }else{
                double *x = base->getLoadedPoints();
                mapCanonicalToTransformed( base->getNumDimensions(), base->getNumLoaded(), base->getRule(), x );
                return x;
        }
}
double* TasmanianSparseGrid::getNeededPoints() const{
        if ( domain_transform_a == 0 ){
                return base->getNeededPoints();
        }else{
                double *x = base->getNeededPoints();
                mapCanonicalToTransformed( base->getNumDimensions(), base->getNumNeeded(), base->getRule(), x );
                return x;
        }
}
double* TasmanianSparseGrid::getPoints() const{
        if ( domain_transform_a == 0 ){
                return base->getPoints();
        }else{
                double *x = base->getPoints();
                mapCanonicalToTransformed( base->getNumDimensions(), base->getNumPoints(), base->getRule(), x );
                return x;
        }
}

double* TasmanianSparseGrid::getQuadratureWeights() const{
        if ( domain_transform_a == 0 ){
                return base->getQuadratureWeights();
        }else{
                double *w = base->getQuadratureWeights();
                double scale = getQuadratureScale( base->getNumDimensions(), base->getRule() );
                #pragma omp parallel for schedule(static)
                for( int i=0; i<getNumPoints(); i++ ) w[i] *= scale;
                return w;
        }
}
double* TasmanianSparseGrid::getInterpolationWeights( const double x[] ) const{
        if ( domain_transform_a == 0 ){
                return base->getInterpolationWeights( x );
        }else{
                int num_dimensions = base->getNumDimensions();
                double *x_canonical = new double[num_dimensions];  std::copy( x, x + num_dimensions, x_canonical );
                mapTransformedToCanonical( num_dimensions, base->getRule(), x_canonical );
                double *w = base->getInterpolationWeights( x_canonical );
                delete[] x_canonical;
                return w;
        }
}

void TasmanianSparseGrid::loadNeededPoints( const double *vals ){ base->loadNeededPoints( vals ); }

void TasmanianSparseGrid::evaluate( const double x[], double y[] ) const{
        if ( domain_transform_a == 0 ){
                base->evaluate( x, y );
        }else{
                int num_dimensions = base->getNumDimensions();
                double *x_canonical = new double[num_dimensions];  std::copy( x, x + num_dimensions, x_canonical );
                mapTransformedToCanonical( num_dimensions, base->getRule(), x_canonical );
                base->evaluate( x_canonical, y );
                delete[] x_canonical;
        }
}
void TasmanianSparseGrid::integrate( double q[] ) const{
        base->integrate( q );
        if ( domain_transform_a != 0 ){
                double scale = getQuadratureScale( base->getNumDimensions(), base->getRule() );
                for( int k=0; k<getNumOutputs(); k++ ) q[k] *= scale;
        }
}

bool TasmanianSparseGrid::isGlobal() const{
        return (global != 0);
}
bool TasmanianSparseGrid::isSequence() const{
        return (sequence != 0);
}
bool TasmanianSparseGrid::isLocalPolynomial() const{
        return (pwpoly != 0);
}
bool TasmanianSparseGrid::isWavelet() const{
        return (wavelet != 0);
}

void TasmanianSparseGrid::setDomainTransform( const double a[], const double b[] ){
        if ( (base == 0) || (base->getNumDimensions() == 0) ){
                cerr << "ERROR: cannot call setDomainTransform on uninitialized grid!" << endl;
                return;
        }
        clearDomainTransform();
        int num_dimensions = base->getNumDimensions();
        domain_transform_a = new double[num_dimensions];  std::copy( a, a + num_dimensions, domain_transform_a );
        domain_transform_b = new double[num_dimensions];  std::copy( b, b + num_dimensions, domain_transform_b );
}
bool TasmanianSparseGrid::isSetDomainTransfrom() const{
        return ( domain_transform_a != 0 );
}
void TasmanianSparseGrid::clearDomainTransform(){
        if ( domain_transform_a != 0 ){ delete[] domain_transform_a; domain_transform_a = 0; }
        if ( domain_transform_b != 0 ){ delete[] domain_transform_b; domain_transform_b = 0; }
}
void TasmanianSparseGrid::getDomainTransform( double a[], double b[] ) const{
        if ( (base == 0) || (base->getNumDimensions() == 0) || (domain_transform_a == 0) ){
                cerr << "ERROR: cannot call getDomainTransform on uninitialized grid or if no transform has been set!" << endl;
                return;
        }
        int num_dimensions = base->getNumDimensions();
        std::copy( domain_transform_a, domain_transform_a + num_dimensions, a );
        std::copy( domain_transform_b, domain_transform_b + num_dimensions, b );
}

void TasmanianSparseGrid::mapCanonicalToTransformed( int num_dimensions, int num_points, TypeOneDRule rule, double x[] ) const{
        if ( rule == rule_gausslaguerre ){
                for( int i=0; i<num_points; i++ ){
                        for( int j=0; j<num_dimensions; j++ ){
                                x[i*num_dimensions+j] /= domain_transform_b[j];
                                x[i*num_dimensions+j] += domain_transform_a[j];
                        }
                }
        }else if ( (rule == rule_gausshermite) || (rule == rule_gausshermiteodd) ){
                for( int i=0; i<num_points; i++ ){
                        for( int j=0; j<num_dimensions; j++ ){
                                x[i*num_dimensions+j] /= sqrt( domain_transform_b[j] );
                                x[i*num_dimensions+j] += domain_transform_a[j];
                        }
                }
        }else{
                for( int i=0; i<num_points; i++ ){
                        for( int j=0; j<num_dimensions; j++ ){
                                x[i*num_dimensions+j] *= 0.5* ( domain_transform_b[j] - domain_transform_a[j] );
                                x[i*num_dimensions+j] += 0.5* ( domain_transform_b[j] + domain_transform_a[j] );
                        }
                }
        }
}
void TasmanianSparseGrid::mapTransformedToCanonical( int num_dimensions, TypeOneDRule rule, double x[] ) const{
        if ( rule == rule_gausslaguerre ){
                for( int j=0; j<num_dimensions; j++ ){
                        x[j] -= domain_transform_a[j];
                        x[j] *= domain_transform_b[j];
                }
        }else if ( (rule == rule_gausshermite) || (rule == rule_gausshermiteodd) ){
                for( int j=0; j<num_dimensions; j++ ){
                        x[j] -= domain_transform_a[j];
                        x[j] *= sqrt( domain_transform_b[j] );
                }
        }else{
                for( int j=0; j<num_dimensions; j++ ){
                        x[j] *= 2.0;
                        x[j] -= ( domain_transform_b[j] + domain_transform_a[j] );
                        x[j] /= ( domain_transform_b[j] - domain_transform_a[j] );
                }
        }
}
double TasmanianSparseGrid::getQuadratureScale( int num_dimensions, TypeOneDRule rule ) const{
        double scale = 1.0;
        if ( (rule == rule_gausschebyshev1) || (rule == rule_gausschebyshev2) || (rule == rule_gaussgegenbauer) || (rule == rule_gaussjacobi) ){
                double alpha = ( rule == rule_gausschebyshev1 ) ? -0.5 : ( rule == rule_gausschebyshev2 ) ? 0.5 : global->getAlpha();
                double beta = ( rule == rule_gausschebyshev1 ) ? -0.5 : ( rule == rule_gausschebyshev2 ) ? 0.5 : ( rule == rule_gaussgegenbauer ) ? global->getAlpha() : global->getBeta();
                for( int j=0; j<num_dimensions; j++ ) scale *= pow( 0.5*( domain_transform_b[j] - domain_transform_a[j] ), alpha + beta + 1.0 );
        }else if ( rule == rule_gausslaguerre ){
                for( int j=0; j<num_dimensions; j++ ) scale *= pow( domain_transform_b[j], global->getAlpha() + 1.0 );
        }else if ( (rule == rule_gausshermite) || (rule == rule_gausshermiteodd) ){
                double power = -0.5 * ( 1.0 + global->getAlpha() );
                for( int j=0; j<num_dimensions; j++ ) scale *= pow( domain_transform_b[j], power );
        }else{
                for( int j=0; j<num_dimensions; j++ ) scale *= (domain_transform_b[j] - domain_transform_a[j] ) / 2.0;
        }
        return scale;
}

void TasmanianSparseGrid::setAnisotropicRefinement( TypeDepth type, int min_growth, int output ){
        if ( sequence != 0 ){
                sequence->setAnisotropicRefinement( type, min_growth );
        }else if ( global != 0 ){
                if ( OneDimensionalMeta::isNonNested( global->getRule() ) ){
                        cerr << "ERROR: setAnisotropicRefinement called for a global grid with non-nested rule" << endl;
                }else{
                        global->setAnisotropicRefinement( type, min_growth, output );
                }
        }else{
                cerr << "ERROR: setAnisotropicRefinement called for grid that is neither sequence nor Global with sequence rule" << endl;
        }
}
int* TasmanianSparseGrid::estimateAnisotropicCoefficients( TypeDepth type, int output ){
        if ( sequence != 0 ){
                return sequence->estimateAnisotropicCoefficients( type, output );
        }else if ( global != 0 ){
                if ( OneDimensionalMeta::isNonNested( global->getRule() ) ){
                        cerr << "ERROR: estimateAnisotropicCoefficients called for a global grid with non-nested rule" << endl;
                }else{
                        return global->estimateAnisotropicCoefficients( type, output );
                }
        }else{
                cerr << "ERROR: estimateAnisotropicCoefficients called for grid that is neither sequence nor Global with sequence rule" << endl;
        }
        return 0;
}
void TasmanianSparseGrid::setSurplusRefinement( double tolerance, int output ){
        if ( sequence != 0 ){
                sequence->setSurplusRefinement( tolerance, output );
        }else if ( global != 0 ){
                if ( OneDimensionalMeta::isSequence( global->getRule() ) ){
                        global->setSurplusRefinement( tolerance, output );
                }else{
                        cerr << "ERROR: setSurplusRefinement called for a global grid with non-sequence rule" << endl;
                }
        }else{
                cerr << "ERROR: setSurplusRefinement( double, int ) called for grid that is neither sequence nor Global with sequence rule" << endl;
        }
}
void TasmanianSparseGrid::setSurplusRefinement( double tolerance, TypeRefinement criteria, int output ){
        if ( pwpoly != 0 ){
                pwpoly->setSurplusRefinement( tolerance, criteria, output );
        }else if ( wavelet != 0 ){
                wavelet->setSurplusRefinement( tolerance, criteria, output );
        }else{
                cerr << "ERROR: setSurplusRefinement( double, TypeRefinement ) called for grid that is neither local polynomial nor wavelet" << endl;
        }
}
void TasmanianSparseGrid::clearRefinement(){
        base->clearRefinement();
}
void TasmanianSparseGrid::removePointsBySurplus( double tolerance, int output ){
        if ( pwpoly == 0 ){
                cout << "ERROR: removePointsBySurplus() called for a grid that is not local polynomial." << endl;
        }else{
                if ( pwpoly->removePointsBySurplus( tolerance, output ) == 0 ){
                        clear();
                }
        }
}

double* TasmanianSparseGrid::evalHierarchicalFunctions( const double x[] ) const{
        return base->evalHierarchicalFunctions( x );
}
void TasmanianSparseGrid::setHierarchicalCoefficients( const double c[] ){
        base->setHierarchicalCoefficients( c );
}

int* TasmanianSparseGrid::getGlobalPolynomialSpace( bool interpolation, int &num_indexes ) const{
        if ( global != 0 ){
                return global->getPolynomialSpace( interpolation, num_indexes );
        }else if ( sequence != 0 ){
                return sequence->getPolynomialSpace( interpolation, num_indexes );
        }else{
                cout << "ERROR: getGlobalPolynomialSpace() called for a grid that is neither Global nor Sequence." << endl;
                num_indexes = 0;
                return 0;
        }
}
const double* TasmanianSparseGrid::getSurpluses() const{
        if ( pwpoly != 0 ){
                return pwpoly->getSurpluses();
        }else if ( wavelet != 0 ){
                return wavelet->getSurpluses();
        }else if ( sequence != 0 ){
                return sequence->getSurpluses();
        }else{
                cout << "ERROR: getSurplusses() called for a grid that is neither local polynomial nor wavelet nor sequence." << endl;
                return 0;
        }
}
const int* TasmanianSparseGrid::getPointsIndexes() const{
        if ( pwpoly != 0 ){
                //return ( (pwpoly->getNumNeeded()>0) ? pwpoly->getNeededIndexes() : pwpoly->getPointIndexes() );
                return pwpoly->getPointIndexes();
        }else if ( wavelet != 0 ){
                return wavelet->getPointIndexes();
        }else if ( global != 0 ){
                return global->getPointIndexes();
        }else if ( sequence != 0 ){
                return sequence->getPointIndexes();
        }else{
                cout << "ERROR: getPointIndexes() called for a grid that is neither local polynomial nor wavelet nor sequence." << endl;
                return 0;
        }
}
const int* TasmanianSparseGrid::getNeededIndexes() const{
        if ( pwpoly != 0 ){
                //return ( (pwpoly->getNumNeeded()>0) ? pwpoly->getNeededIndexes() : pwpoly->getPointIndexes() );
                return pwpoly->getNeededIndexes();
        }else{
                cout << "ERROR: getPointIndexes() called for a grid that is neither local polynomial nor wavelet nor sequence." << endl;
                return 0;
        }
}

void TasmanianSparseGrid::printStats() const{
        using std::setw;

        const int L1 = 20, L2 = 10;
        cout << endl;
        cout << setw(L1) << "Grid Type:" << setw(L2) << " " << " ";
        if ( isGlobal() ) cout << "Global";
        if ( isSequence() ) cout << "Sequence";
        if ( isLocalPolynomial() ) cout << "Local Polynomial";
        if ( isWavelet() ) cout << "Wavelets";
        if ( !(isGlobal() || isSequence() || isLocalPolynomial() || isWavelet()) ) cout << "none";
        cout << endl;

        cout << setw(L1) << "Dimensions:" << setw(L2) << " " << getNumDimensions() << endl;
        cout << setw(L1) << "Outputs:" << setw(L2) << " " << getNumOutputs() << endl;
        if ( getNumOutputs() == 0 ){
                cout << setw(L1) << "Nodes:" << setw(L2) << " " << getNumPoints() << endl;
        }else{
                cout << setw(L1) << "Loaded nodes:" << setw(L2) << " " << getNumLoaded() << endl;
                cout << setw(L1) << "Needed nodes:" << setw(L2) << " " << getNumNeeded() << endl;
        }
        cout << setw(L1) << "Rule:" << setw(L2) << " " << " " << OneDimensionalMeta::getHumanString( getRule() ) << endl;
        if ( getRule() == rule_customtabulated ){
                cout << setw(L1) << "Description:" << setw(L2) << " " << " " << getCustomRuleDescription() << endl;
        }
        if ( isSetDomainTransfrom() ){
                cout << setw(L1) << "Domain:" << setw(L2) << "Custom" << endl;
        }else{
                cout << setw(L1) << "Domain:" << setw(L2) << "Canonical" << endl;
        }

        if ( isGlobal() ){
                TypeOneDRule rr = getRule();
                if ( (rr == rule_gaussgegenbauer) || (rr == rule_gausslaguerre) || (rr == rule_gausshermite) || (rr == rule_gaussgegenbauerodd) || (rr == rule_gausshermiteodd)  ){
                        cout << setw(L1) << "Alpha:" << setw(L2) << " " << getAlpha() << endl;
                }
                if ( rr == rule_gaussjacobi ){
                        cout << setw(L1) << "Alpha:" << setw(L2) << " " << getAlpha() << endl;
                        cout << setw(L1) << "Beta:" << setw(L2) << " " << getBeta() << endl;
                }
        }else if ( isSequence() ){
                // sequence rules are simple, nothing to specify here
        }else if ( isLocalPolynomial() ){
                int o = getOrder();
                cout << setw(L1) << "Order:" << setw(L2) << " " << o << endl;
        }else if ( isWavelet() ){
                int o = getOrder();
                cout << setw(L1) << "Order:" << setw(L2) << " " << o << endl;
        }else{
                // cerr << endl << "ERROR: unknown grid type!" << endl;
                // show nothing, just like the sequence grid
        }

        cout << endl;
}


// ------------ C Interface for use with Python ctypes and potentially other C codes -------------- //
extern "C" void* tsgConstructTasmanianSparseGrid(){  return (void*) new TasmanianSparseGrid();  }
extern "C" void tsgDestructTasmanianSparseGrid( void * grid ){  delete ((TasmanianSparseGrid*) grid);  }

extern "C" void tsgCopyGrid( void * destination, void * source ){  ((TasmanianSparseGrid*) destination)->copyGrid( ((TasmanianSparseGrid*) source) );  }

extern "C" const char* tsgGetVersion(){  return TasmanianSparseGrid::getVersion();  }
extern "C" const char* tsgGetLicense(){  return TasmanianSparseGrid::getLicense();  }

extern "C" void tsgWrite( void * grid, const char* filename ){  ((TasmanianSparseGrid*) grid)->write( filename );  }
extern "C" int tsgRead( void * grid, const char* filename ){
        bool result = ((TasmanianSparseGrid*) grid)->read( filename );
        return result ? 0 : 1;
}

extern "C" void tsgMakeGlobalGrid( void * grid, int dimensions, int outputs, int depth, const char * sType, const char * sRule, int * anisotropic_weights, double alpha, double beta, const char* custom_filename ){
        TypeDepth depth_type = OneDimensionalMeta::getIOTypeString( sType );
        TypeOneDRule rule = OneDimensionalMeta::getIORuleString( sRule );
        if ( depth_type == type_none ){  cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl;  }
        if ( rule == rule_none ){  cerr << "WARNING: incorrect rule type: " << sType << ", defaulting to clenshaw-curtis." << endl;  }
        ((TasmanianSparseGrid*) grid)->makeGlobalGrid( dimensions, outputs, depth, depth_type, rule, anisotropic_weights, alpha, beta, custom_filename );
}
extern "C" void tsgMakeSequenceGrid( void * grid, int dimensions, int outputs, int depth, const char * sType, const char * sRule, int * anisotropic_weights ){
        TypeDepth depth_type = OneDimensionalMeta::getIOTypeString( sType );
        TypeOneDRule rule = OneDimensionalMeta::getIORuleString( sRule );
        if ( depth_type == type_none ){  cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl;  depth_type = type_iptotal;  }
        if ( rule == rule_none ){  cerr << "WARNING: incorrect rule type: " << sRule << ", defaulting to clenshaw-curtis." << endl;  rule = rule_clenshawcurtis;  }
        ((TasmanianSparseGrid*) grid)->makeSequenceGrid( dimensions, outputs, depth, depth_type, rule, anisotropic_weights );
}
extern "C" void tsgMakeLocalPolynomialGrid( void * grid, int dimensions, int outputs, int depth, int order, const char * sRule ){
        TypeOneDRule rule = OneDimensionalMeta::getIORuleString( sRule );
        if ( rule == rule_none ){  cerr << "WARNING: incorrect rule type: " << sRule << ", defaulting to localp." << endl;  rule = rule_localp;  }
        ((TasmanianSparseGrid*) grid)->makeLocalPolynomialGrid( dimensions, outputs, depth, order, rule );
}
extern "C" void tsgMakeWaveletGrid( void * grid, int dimensions, int outputs, int depth, int order ){
        ((TasmanianSparseGrid*) grid)->makeWaveletGrid( dimensions, outputs, depth, order );
}

extern "C" void tsgUpdateGlobalGrid( void * grid, int depth, const char * sType, const int *anisotropic_weights ){
        TypeDepth depth_type = OneDimensionalMeta::getIOTypeString( sType );
        if ( depth_type == type_none ){  cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl;  depth_type = type_iptotal;  }
        ((TasmanianSparseGrid*) grid)->updateGlobalGrid( depth, depth_type, anisotropic_weights );
}
extern "C" void tsgUpdateSequenceGrid( void * grid, int depth, const char * sType, const int *anisotropic_weights ){
        TypeDepth depth_type = OneDimensionalMeta::getIOTypeString( sType );
        if ( depth_type == type_none ){  cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl;  depth_type = type_iptotal;  }
        ((TasmanianSparseGrid*) grid)->updateSequenceGrid( depth, depth_type, anisotropic_weights );
}

extern "C" double tsgGetAlpha( void * grid ){  return ((TasmanianSparseGrid*) grid)->getAlpha();  }
extern "C" double tsgGetBeta( void * grid ){  return ((TasmanianSparseGrid*) grid)->getBeta();  }
extern "C" int tsgGetOrder( void * grid ){  return ((TasmanianSparseGrid*) grid)->getOrder();  }
extern "C" int tsgGetNumDimensions( void * grid ){  return ((TasmanianSparseGrid*) grid)->getNumDimensions();  }
extern "C" int tsgGetNumOutputs( void * grid ){  return ((TasmanianSparseGrid*) grid)->getNumOutputs();  }
extern "C" const char * tsgGetRule( void * grid ){  return  OneDimensionalMeta::getIORuleString(  ((TasmanianSparseGrid*) grid)->getRule()  );  }
extern "C" const char * tsgGetCustomRuleDescription( void * grid ){  return ((TasmanianSparseGrid*) grid)->getCustomRuleDescription();  }

extern "C" int tsgGetNumLoaded( void * grid ){  return ((TasmanianSparseGrid*) grid)->getNumLoaded();  }
extern "C" int tsgGetNumNeeded( void * grid ){  return ((TasmanianSparseGrid*) grid)->getNumNeeded();  }
extern "C" int tsgGetNumPoints( void * grid ){  return ((TasmanianSparseGrid*) grid)->getNumPoints();  }

extern "C" double* tsgGetLoadedPoints( void * grid ){  return ((TasmanianSparseGrid*) grid)->getLoadedPoints();  }
extern "C" double* tsgGetNeededPoints( void * grid ){  return ((TasmanianSparseGrid*) grid)->getNeededPoints();  }
extern "C" double* tsgGetPoints( void * grid ){  return ((TasmanianSparseGrid*) grid)->getPoints();  }

extern "C" double* tsgGetQuadratureWeights( void * grid ){  return ((TasmanianSparseGrid*) grid)->getQuadratureWeights();  }
extern "C" double* tsgGetInterpolationWeights( void * grid, const double *x ){  return ((TasmanianSparseGrid*) grid)->getInterpolationWeights( x );  }

extern "C" void tsgLoadNeededPoints( void * grid, const double *vals ){  ((TasmanianSparseGrid*) grid)->loadNeededPoints( vals );  }

extern "C" void tsgEvaluate( void * grid, const double *x, double *y ){  ((TasmanianSparseGrid*) grid)->evaluate( x, y );  }
extern "C" void tsgIntegrate( void * grid, double *q ){  ((TasmanianSparseGrid*) grid)->integrate( q );  }

extern "C" void tsgBatchEvaluate( void * grid, const double *x, int num_x, double *y ){
        TasmanianSparseGrid* tsg = (TasmanianSparseGrid*) grid;
        int iNumDim = tsg->getNumDimensions(), iNumOutputs = tsg->getNumOutputs();
        #pragma omp parallel for
        for( int i=0; i<num_x; i++ ){
                tsg->evaluate( &(x[ i*iNumDim ]), &(y[ i*iNumOutputs ]) );
        }
}

extern "C" double* tsgBatchGetInterpolationWeights( void * grid, const double *x, int num_x ){
        TasmanianSparseGrid* tsg = (TasmanianSparseGrid*) grid;
        int iNumDim = tsg->getNumDimensions(), iNumPoints = tsg->getNumPoints();
        double *weights = new double[ num_x * iNumPoints ];
        #pragma omp parallel for
        for( int i=0; i<num_x; i++ ){
                double *w = tsg->getInterpolationWeights( &(x[ i*iNumDim ]) );
                std::copy( w, w + iNumPoints, &( weights[ i*iNumPoints ] ) );
                delete[] w;
        }
        return weights;
}

extern "C" int tsgIsGlobal( void * grid ){  return ( ((TasmanianSparseGrid*) grid)->isGlobal() ? 0 : 1 );  }
extern "C" int tsgIsSequence( void * grid ){  return ( ((TasmanianSparseGrid*) grid)->isSequence() ? 0 : 1 );  }
extern "C" int tsgIsLocalPolynomial( void * grid ){  return ( ((TasmanianSparseGrid*) grid)->isLocalPolynomial() ? 0 : 1 );  }
extern "C" int tsgIsWavelet( void * grid ){  return ( ((TasmanianSparseGrid*) grid)->isWavelet() ? 0 : 1 );  }

extern "C" void tsgSetDomainTransform( void * grid, const double a[], const double b[] ){  ((TasmanianSparseGrid*) grid)->setDomainTransform( a, b );  }
extern "C" bool tsgIsSetDomainTransfrom( void * grid ){  return ( ((TasmanianSparseGrid*) grid)->isSetDomainTransfrom() ? 0 : 1 );  }
extern "C" void tsgClearDomainTransform( void * grid ){  ((TasmanianSparseGrid*) grid)->clearDomainTransform();  }
extern "C" void tsgGetDomainTransform( void * grid, double a[], double b[] ){  ((TasmanianSparseGrid*) grid)->getDomainTransform( a, b );  }

extern "C" void tsgSetAnisotropicRefinement( void * grid, const char * sType, int min_growth, int output ){
        TypeDepth depth_type = OneDimensionalMeta::getIOTypeString( sType );
        if ( depth_type == type_none ){  cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl;  depth_type = type_iptotal;  }
        ((TasmanianSparseGrid*) grid)->setAnisotropicRefinement( depth_type, min_growth, output );
}
extern "C" int* tsgEstimateAnisotropicCoefficients( void * grid, const char * sType, int output, int *num_coefficients ){
        TypeDepth depth_type = OneDimensionalMeta::getIOTypeString( sType );
        if ( depth_type == type_none ){  cerr << "WARNING: incorrect depth type: " << sType << ", defaulting to type_iptotal." << endl;  depth_type = type_iptotal;  }
        if ( (depth_type == type_curved) || (depth_type == type_ipcurved) || (depth_type == type_qpcurved) ){
                *num_coefficients = 2 * ( ((TasmanianSparseGrid*) grid)->getNumDimensions() );
        }else{
                *num_coefficients = ((TasmanianSparseGrid*) grid)->getNumDimensions();
        }
        return ((TasmanianSparseGrid*) grid)->estimateAnisotropicCoefficients( depth_type, output );
}
extern "C" void tsgSetGlobalSurplusRefinement( void * grid, double tolerance, int output ){
        ((TasmanianSparseGrid*) grid)->setSurplusRefinement( tolerance, output );
}
extern "C" void tsgSetLocalSurplusRefinement( void * grid, double tolerance, const char * sRefinementType, int output ){
        TypeRefinement ref_type = OneDimensionalMeta::getIOTypeRefinementString( sRefinementType );
        if ( ref_type == refine_none ){  cerr << "WARNING: incorrect refinement type: " << sRefinementType << ", defaulting to type_classic." << endl;  ref_type = refine_classic;  }
        ((TasmanianSparseGrid*) grid)->setSurplusRefinement( tolerance, ref_type, output );
}
extern "C" void tsgClearRefinement( void * grid ){
        ((TasmanianSparseGrid*) grid)->clearRefinement();
}
extern "C" void tsgRemovePointsBySurplus( void * grid, double tolerance, int output ){
        ((TasmanianSparseGrid*) grid)->removePointsBySurplus( tolerance, output );
}

extern "C" double* tsgEvalHierarchicalFunctions( void * grid, const double *x ){
        return ((TasmanianSparseGrid*) grid)->evalHierarchicalFunctions( x );
}
extern "C" double* tsgBatchEvalHierarchicalFunctions( void * grid, const double *x, int num_x ){
        TasmanianSparseGrid* tsg = (TasmanianSparseGrid*) grid;
        int iNumDim = tsg->getNumDimensions(), iNumPoints = tsg->getNumPoints();
        double *vals = new double[ num_x * iNumPoints ];
        #pragma omp parallel for
        for( int i=0; i<num_x; i++ ){
                double *v = tsg->evalHierarchicalFunctions( &(x[ i*iNumDim ]) );
                std::copy( v, v + iNumPoints, &( vals[ i*iNumPoints ] ) );
                delete[] v;
        }
        return vals;
}
extern "C" void tsgSetHierarchicalCoefficients( void * grid, const double *c ){
        ((TasmanianSparseGrid*) grid)->setHierarchicalCoefficients( c );
}
extern "C" const double* tsgGetSurpluses( void *grid ){
        return ((TasmanianSparseGrid*) grid)->getSurpluses();
}

extern "C" int* tsgGetGlobalPolynomialSpace( void * grid, int interpolation, int *num_indexes ){
        int ni, *idx = ((TasmanianSparseGrid*) grid)->getGlobalPolynomialSpace( (interpolation == 0), ni );
        *num_indexes = ni;
        return idx;
}


extern "C" void tsgPrintStats( void *grid ){  ((TasmanianSparseGrid*) grid)->printStats();  }

extern "C" void tsgDeleteDoubles( double * p ){  delete[] p;  }
extern "C" void tsgDeleteInts( int * p ){  delete[] p;  }

}

#endif

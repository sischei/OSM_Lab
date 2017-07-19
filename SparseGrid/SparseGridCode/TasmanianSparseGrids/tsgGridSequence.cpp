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

#ifndef __TASMANIAN_SPARSE_GRID_GLOBAL_NESTED_CPP
#define __TASMANIAN_SPARSE_GRID_GLOBAL_NESTED_CPP

#include "tsgGridSequence.hpp"

namespace TasGrid{

GridSequence::GridSequence() : num_dimensions(0), num_outputs(0), points(0), needed(0), parents(0), surpluses(0), nodes(0), coeff(0), values(0), max_levels(0)
{}
GridSequence::GridSequence( const GridSequence &seq ) : num_dimensions(0), num_outputs(0), points(0), needed(0), parents(0), surpluses(0), nodes(0), coeff(0), values(0), max_levels(0)
{  copyGrid( &seq ); }
GridSequence::~GridSequence(){ reset(); }

void GridSequence::write( std::ofstream &ofs ) const{
        ofs << std::scientific; ofs.precision(17);
        ofs << num_dimensions << " " << num_outputs << endl;
        if ( num_dimensions > 0 ){
                ofs << OneDimensionalMeta::getIORuleString( rule ) << endl;
                if ( points == 0 ){
                        ofs << "0" << endl;
                }else{
                        ofs << "1 ";
                        points->write( ofs );
                }
                if ( needed == 0 ){
                        ofs << "0" << endl;
                }else{
                        ofs << "1 ";
                        needed->write( ofs );
                }
                if ( surpluses == 0 ){
                        ofs << "0";
                }else{
                        ofs << "1";
                        for( int i=0; i<num_outputs * points->getNumIndexes(); i++ ) ofs << " " << surpluses[i];
                }
                ofs << endl;
                values->write( ofs );
        }
}
void GridSequence::read( std::ifstream &ifs ){
        reset();
        ifs >> num_dimensions >> num_outputs;
        if ( num_dimensions > 0 ){
                int flag;
                std::string T;
                ifs >> T;
                rule = OneDimensionalMeta::getIORuleString( T.c_str() );
                ifs >> flag;  if ( flag == 1 ){  points = new IndexSet(num_dimensions); points->read( ifs );  }
                ifs >> flag;  if ( flag == 1 ){  needed = new IndexSet(num_dimensions); needed->read( ifs );  }
                IndexManipulator IM(num_dimensions);
                parents = IM.computeDAGup( (( points == 0 ) ? needed : points) );
                ifs >> flag;  if ( flag == 1 ){  int n = num_outputs * points->getNumIndexes(); surpluses = new double[n]; for( int i=0; i<n; i++ ) ifs >> surpluses[i];  }
                values = new StorageSet( 0, 0 ); values->read( ifs );
                int mp = 0, mn = 0, max_level;
                max_levels = new int[num_dimensions];
                if ( needed == 0 ){ // points must be non-zero
                        IM.getMaxLevels( points, max_levels, mp );
                }else if ( points == 0 ){ // only needed, no points (right after creation)
                        IM.getMaxLevels( needed, max_levels, mn );
                }else{ // both points and needed are set
                        IM.getMaxLevels( points, max_levels, mp );
                        IM.getMaxLevels( needed, 0, mn );
                }
                max_level = ( mp > mn ) ? mp : mn;
                prepareSequence( max_level + 1 );
        }
}

void GridSequence::reset(){
        if ( points != 0 ){ delete points; points = 0; }
        if ( needed != 0 ){ delete needed; needed = 0; }
        if ( parents != 0 ){ delete[] parents; parents = 0; }
        if ( surpluses != 0 ){ delete[] surpluses; surpluses = 0; }
        if ( nodes != 0 ){ delete[] nodes; nodes = 0; }
        if ( coeff != 0 ){ delete[] coeff; coeff = 0; }
        if ( values != 0 ){ delete values; values = 0; }
        if ( max_levels != 0 ){ delete[] max_levels; max_levels = 0; }
}
void GridSequence::clearRefinement(){
        if ( needed != 0 ){ delete needed; needed = 0; }
}

void GridSequence::makeGrid( int cnum_dimensions, int cnum_outputs, int depth, TypeDepth type, TypeOneDRule crule, const int *anisotropic_weights ){
        IndexManipulator IM(cnum_dimensions);

        IndexSet *pset = IM.selectTensors( depth, type, anisotropic_weights, crule );
        setPoints( pset, cnum_outputs, crule );
}
void GridSequence::copyGrid( const GridSequence *seq ){
        IndexSet *work = ( seq->points == 0 ) ? seq->needed : seq->points;
        IndexSet *pset = new IndexSet( work );
        setPoints( pset, seq->num_outputs, seq->rule );
        if ( (num_outputs > 0) && (seq->points != 0) ){ // if there are values inside the source object
                loadNeededPoints( seq->values->getValues( 0 ) );
        }

        if ( (seq->points != 0) && (seq->needed != 0) ){ // there is a refinement
                IndexManipulator IM(num_dimensions);
                int m1, m2;
                IM.getMaxLevels( seq->points, 0, m1 );
                IM.getMaxLevels( seq->needed, 0, m2 );
                int max_level = (m1 > m2) ? m1 : m2;
                delete[] nodes; delete[] coeff;
                prepareSequence( max_level+1 );

                needed = new IndexSet( seq->needed );
        }
}

void GridSequence::setPoints( IndexSet* &pset, int cnum_outputs, TypeOneDRule crule ){
        reset();
        num_dimensions = pset->getNumDimensions();
        num_outputs = cnum_outputs;
        rule = crule;

        IndexManipulator IM(num_dimensions);

        needed = pset;
        pset = 0;

        parents = IM.computeDAGup( needed );

        max_levels = new int[num_dimensions];
        int max_level; IM.getMaxLevels( needed, max_levels, max_level );

        prepareSequence( max_level + 1 );

        if ( num_outputs == 0 ){
                points = needed;
                needed = 0;
        }else{
                values = new StorageSet( num_outputs, needed->getNumIndexes() );
        }
}

void GridSequence::updateGrid( int depth, TypeDepth type, const int *anisotropic_weights ){
        IndexManipulator IM(num_dimensions);

        IndexSet *pset = IM.selectTensors( depth, type, anisotropic_weights, rule );
        updateGrid( pset );
}

void GridSequence::updateGrid( IndexSet* &update ){
        clearRefinement();
        if ( (num_outputs == 0) || (points == 0) ){
                setPoints( update, num_outputs, rule );
        }else{
                update->addIndexSet( points );
                needed = update->diffSets( points );

                if ( (needed != 0) && (needed->getNumIndexes() > 0) ){
                        OneDimensionalMeta meta;
                        IndexManipulator IM(num_dimensions);
                        int max_level; IM.getMaxLevels( update, 0, max_level );
                        delete[] nodes; delete[] coeff;
                        prepareSequence( max_level+1 );
                }else if ( needed->getNumIndexes() == 0 ){
                        delete needed;
                }

                delete update;
        }
}

int GridSequence::getNumDimensions() const{ return num_dimensions; }
int GridSequence::getNumOutputs() const{ return num_outputs; }
TypeOneDRule GridSequence::getRule() const{ return rule; }

int GridSequence::getNumLoaded() const{ return ((points == 0) ? 0 : points->getNumIndexes()); }
int GridSequence::getNumNeeded() const{ return ((needed == 0) ? 0 : needed->getNumIndexes()); }
int GridSequence::getNumPoints() const{ return ((points == 0) ? getNumNeeded() : getNumLoaded()); }

double* GridSequence::getLoadedPoints() const{
        if ( points == 0 ) return 0;
        int num_points = points->getNumIndexes();
        if ( num_points == 0 ) return 0;
        double *x = new double[ num_dimensions * num_points ];
        #pragma omp parallel for
        for( int i=0; i<num_points; i++ ){
                const int *p = points->getIndex( i );
                for( int j=0; j<num_dimensions; j++ ){
                        x[i*num_dimensions + j] = nodes[ p[j] ];
                }
        }
        return x;
}
double* GridSequence::getNeededPoints() const{
        if ( needed == 0 ) return 0;
        int num_points = needed->getNumIndexes();
        if ( num_points == 0 ) return 0;
        double *x = new double[ num_dimensions * num_points ];
        #pragma omp parallel for
        for( int i=0; i<num_points; i++ ){
                const int *p = needed->getIndex( i );
                for( int j=0; j<num_dimensions; j++ ){
                        x[i*num_dimensions + j] = nodes[ p[j] ];
                }
        }
        return x;
}

double* GridSequence::getPoints() const{
        return (( points == 0 ) ? getNeededPoints() : getLoadedPoints());
}

double* GridSequence::getQuadratureWeights() const{
        IndexSet *work = ( points == 0 ) ? needed : points;
        double *integ = cacheBasisIntegrals();

        int n = work->getNumIndexes();
        double *weights = new double[n];

        for( int i=0; i<n; i++ ){
                const int* p = work->getIndex(i);
                weights[i] = integ[ p[0] ];
                for( int j=1; j<num_dimensions; j++ ){
                        weights[i] *= integ[ p[j] ];
                }
        }
        delete[] integ;

        applyTransformationTransposed( weights );

        return weights;
}

double* GridSequence::getInterpolationWeights( const double x[] ) const{
        double **cache = cacheBasisValues( x );

        IndexSet *work = ( points == 0 ) ? needed : points;
        int n = work->getNumIndexes();

        double *weights = new double[n];

        weights[0] = 1.0;
        for( int i=1; i<n; i++ ){
                const int* p = work->getIndex(i);
                weights[i] = cache[0][ p[0] ];
                for( int j=1; j<num_dimensions; j++ ){
                        weights[i] *= cache[j][ p[j] ];
                }
        }

        for( int j=0; j<num_dimensions; j++ ) delete[] cache[j];
        delete[] cache;

        applyTransformationTransposed( weights );

        return weights;
}

void GridSequence::loadNeededPoints( const double *vals ){
        if ( points == 0 ){
                values->setValues( vals );
                points = needed;
                needed = 0;
        }else if ( needed == 0 ){
                values->setValues( vals );
        }else{
                values->addValues( points, needed, vals );
                points->addIndexSet( needed );
                delete needed; needed = 0;
                delete[] parents;
                IndexManipulator IM(num_dimensions);
                parents = IM.computeDAGup( points );
                int m; IM.getMaxLevels( points, max_levels, m );
        }
        recomputeSurpluses();
}

void GridSequence::evaluate( const double x[], double y[] ) const{
        double **cache = cacheBasisValues( x );

        std::fill( y, y + num_outputs, 0.0 );

        int n = points->getNumIndexes();

        for( int i=0; i<n; i++ ){
                const int* p = points->getIndex(i);
                double basis_value = cache[0][ p[0] ];
                for( int j=1; j<num_dimensions; j++ ){
                        basis_value *= cache[j][ p[j] ];
                }

                for( int k=0; k<num_outputs; k++ ){
                        y[k] += basis_value * surpluses[ i*num_outputs + k ];
                }
        }

        for( int j=0; j<num_dimensions; j++ ) delete[] cache[j];
        delete[] cache;
}

void GridSequence::integrate( double q[] ) const{
        double *integ = cacheBasisIntegrals();

        std::fill( q, q + num_outputs, 0.0 );

        int n = points->getNumIndexes();

        for( int i=0; i<n; i++ ){
                const int* p = points->getIndex(i);
                double w = integ[ p[0] ];
                for( int j=1; j<num_dimensions; j++ ){
                        w *= integ[ p[j] ];
                }

                for( int k=0; k<num_outputs; k++ ){
                        q[k] += w * surpluses[ i*num_outputs + k ];
                }
        }
        delete[] integ;
}

double* GridSequence::evalHierarchicalFunctions( const double x[] ) const{
        IndexSet *work = ( points == 0 ) ? needed : points;
        int num_points = work->getNumIndexes();
        double *vals = new double[ num_points ];

        double **cache = cacheBasisValues( x );

        for( int i=0; i<num_points; i++ ){
                const int* p = work->getIndex( i );
                vals[i] = cache[0][ p[0] ];
                for( int j=1; j<num_dimensions; j++ ){
                        vals[i] *= cache[j][ p[j] ];
                }
        }

        for( int j=0; j<num_dimensions; j++ ) delete[] cache[j];
        delete[] cache;

        return vals;
}
void GridSequence::setHierarchicalCoefficients( const double c[] ){
        int num_ponits = points->getNumIndexes();
        if ( surpluses != 0 ) delete[] surpluses;
        surpluses = new double[ num_ponits * num_outputs ];
        std::copy( c, c + num_ponits * num_outputs, surpluses );
}

int* GridSequence::estimateAnisotropicCoefficients( TypeDepth type, int output ) const{
        double tol = 1000 * TSG_NUM_TOL;
        int num_points = points->getNumIndexes();
        double *max_surp = new double[num_points];

        if ( output == -1 ){
                double *norm = new double[num_outputs]; std::fill( norm, norm + num_outputs, 0.0 );
                for( int i=0; i<num_points; i++ ){
                        const double *val = values->getValues( i );
                        for( int k=0; k<num_outputs; k++ ){
                                double v = fabs( val[k] );
                                if ( norm[k] < v ) norm[k] = v;
                        }
                }
                #pragma omp parallel for
                for( int i=0; i<num_points; i++ ){
                        double smax = fabs( surpluses[i*num_outputs] ) / norm[0];
                        for( int k=1; k<num_outputs; k++ ){
                                double v = fabs( surpluses[i*num_outputs + k] ) / norm[k];
                                if ( smax < v ) smax = v;
                        }
                        max_surp[i] = smax;
                }
                delete[] norm;
        }else{
                for( int i=0; i<num_points; i++ ){
                        max_surp[i] = surpluses[i*num_outputs + output];
                }
        }



        int n = 0, m;
        for( int i=0; i<num_points; i++ ){
                n += ( max_surp[i] > tol ) ? 1 : 0;
        }

        double *A, *b;

        if ( (type == type_curved) || (type == type_ipcurved) || (type == type_qpcurved) ){
                m = 2*num_dimensions + 1;
                A = new double[ n * m ];
                b = new double[ n ];

                int count = 0;
                for( int c=0; c<num_points; c++ ){
                        const int *indx = points->getIndex(c);
                        if ( max_surp[c] > tol ){
                                for( int j=0; j<num_dimensions; j++ ){
                                        A[ j*n + count ] = ( (double) indx[j] );
                                }
                                for( int j=0; j<num_dimensions; j++ ){
                                        A[ (num_dimensions + j)*n + count ] = log( (double) (indx[j] + 1 ) );
                                }
                                A[ 2*num_dimensions*n + count ] = 1.0;
                                b[count++] = -log( max_surp[c] );
                        }
                }
        }else{
                m = num_dimensions + 1;
                A = new double[ n * m ];
                b = new double[ n ];

                int count = 0;
                for( int c=0; c<num_points; c++ ){
                        const int *indx = points->getIndex(c);
                        if ( max_surp[c] > tol ){
                                for( int j=0; j<num_dimensions; j++ ){
                                        A[ j*n + count ] = - ( (double) indx[j] );
                                }
                                A[ num_dimensions*n + count ] = 1.0;
                                b[count++] = log( max_surp[c] );
                        }
                }
        }
        delete[] max_surp;

        double *x = new double[m];
        TasmanianDenseSolver::solveLeastSquares( n, m, A, b, 1.E-5, x );

        int *weights = new int[m--];
        for( int j=0; j<m; j++ ){
                weights[j] = (int)( x[j] * 1000.0 + 0.5 );
        }

        int min_weight = weights[0];  for( int j=1; j<num_dimensions; j++ ){  if ( min_weight < weights[j] ) min_weight = weights[j];  } // start by finding the largest weight
        if ( min_weight < 0 ){ // all directions are diverging, default to isotropic total degree
                for( int j=0; j<num_dimensions; j++ )  weights[j] = 1;
                if ( m == 2*num_dimensions ) for( int j=num_dimensions; j<2*num_dimensions; j++ )  weights[j] = 0;
        }else{
                for( int j=0; j<num_dimensions; j++ )  if( (weights[j]>0) && (weights[j]<min_weight) ) min_weight = weights[j];  // find the smallest positive weight
                for( int j=0; j<num_dimensions; j++ ){
                        if ( weights[j] <= 0 ){
                                weights[j] = min_weight;
                                if ( m == 2*num_dimensions ){
                                        if ( std::abs( weights[num_dimensions + j] ) > weights[j] ){
                                                weights[num_dimensions + j] = ( weights[num_dimensions + j] > 0.0 ) ? weights[j] : -weights[j];
                                        }
                                }
                        }
                }
        }

        delete[] A;
        delete[] b;
        delete[] x;

        return weights;
}

void GridSequence::setAnisotropicRefinement( TypeDepth type, int min_growth, int output ){
        clearRefinement();

        int *weights = estimateAnisotropicCoefficients( type, output );

        IndexManipulator IM(num_dimensions);
        int level = IM.getMinChildLevel( points, type, weights, rule );

        IndexSet* total = IM.selectTensors( level, type, weights, rule );
        needed = total->diffSets( points ); // this exploits the 1-1 correspondence between points and tensors

        while( (needed == 0) || (needed->getNumIndexes() < min_growth) ){
                delete total;
                if ( needed != 0 ) delete needed;
                total = IM.selectTensors( ++level, type, weights, rule );
                needed = total->diffSets( points );
        }
        total->addIndexSet( points );
        delete[] weights;

        //OneDimensionalMeta meta;
        int max_level; IM.getMaxLevels( total, 0, max_level );
        delete[] nodes; delete[] coeff;
        prepareSequence( max_level+1 );
        delete total;
}
void GridSequence::setSurplusRefinement( double tolerance, int output ){
        clearRefinement();

        int num_points = points->getNumIndexes();
        bool *flagged = new bool[num_points];

        double *norm = new double[num_outputs]; std::fill( norm, norm + num_outputs, 0.0 );
        for( int i=0; i<num_points; i++ ){
                const double *val = values->getValues( i );
                for( int k=0; k<num_outputs; k++ ){
                        double v = fabs( val[k] );
                        if ( norm[k] < v ) norm[k] = v;
                }
        }

        if ( output == -1 ){
                #pragma omp parallel for
                for( int i=0; i<num_points; i++ ){
                        double smax = fabs( surpluses[i*num_outputs] ) / norm[0];
                        for( int k=1; k<num_outputs; k++ ){
                                double v = fabs( surpluses[i*num_outputs + k] ) / norm[k];
                                if ( smax < v ) smax = v;
                        }
                        flagged[i] = ( smax > tolerance );
                }
        }else{
                for( int i=0; i<num_points; i++ ){
                        flagged[i] = ( (fabs( surpluses[i*num_outputs + output] ) / norm[output]) > tolerance );
                }
        }
        delete[] norm;

        IndexManipulator IM(num_dimensions);
        IndexSet *kids = IM.selectFlaggedChildren( points, flagged );
        if ( (kids != 0) && (kids->getNumIndexes() > 0) ){
                kids->addIndexSet( points );

                IndexSet *total = IM.getLowerCompletion( kids );
                if ( total == 0 ){
                        total = kids;
                        kids = 0;
                }else{
                        total->addIndexSet( kids );
                        delete kids;
                }

                OneDimensionalMeta meta;
                int max_level; IM.getMaxLevels( total, 0, max_level );
                delete[] nodes; delete[] coeff;
                prepareSequence( max_level+1 );

                needed = total->diffSets( points );
                delete total;

        }
        delete[] flagged;
}

int* GridSequence::getPolynomialSpace( bool interpolation, int &n ) const{
        IndexSet *work = ( points == 0 ) ? needed : points;
        int *poly;
        if ( interpolation ){
                n = work->getNumIndexes();

                poly = new int[ n * num_dimensions ];
                const int* p = work->getIndex( 0 );

                std::copy( p, p + n * num_dimensions, poly );
        }else{
                IndexManipulator IM( num_dimensions );
                IndexSet* set = IM.getPolynomialSpace( work, rule, interpolation );

                n = set->getNumIndexes();

                poly = new int[ n * num_dimensions ];
                const int* p = set->getIndex( 0 );

                std::copy( p, p + n * num_dimensions, poly );

                delete set;

                return poly;
        }
        return poly;
}
const double* GridSequence::getSurpluses() const{
        return surpluses;
}
const int* GridSequence::getPointIndexes() const{
        return (( points == 0 ) ? needed->getIndex( 0 ) : points->getIndex( 0 ));
}

void GridSequence::prepareSequence( int n ){
        GreedySequences greedy;
        if ( rule == rule_leja ){
                nodes = greedy.getLejaNodes( n );
        }else if ( rule == rule_maxlebesgue ){
                nodes = greedy.getMaxLebesgueNodes( n );
        }else if ( rule == rule_minlebesgue ){
                nodes = greedy.getMinLebesgueNodes( n );
        }else if ( rule == rule_mindelta ){
                nodes = greedy.getMinDeltaNodes( n );
        }else if ( rule == rule_rleja ){
                OneDimensionalNodes core;
                nodes = core.getRLeja( n );
        }
        coeff = new double[ n ];
        coeff[0] = 1.0;
        for( int i=1; i<n; i++ ){
                coeff[i] = 1.0;
                for( int j=0; j<i; j++ ) coeff[i] *= ( nodes[i] - nodes[j] );
        }
}

double** GridSequence::cacheBasisValues( const double x[] ) const{
        double **cache = new double*[num_dimensions];
        for( int j=0; j<num_dimensions; j++ ){
                cache[j] = new double[ max_levels[j] + 1 ];
                cache[j][0] = 1.0;
                for( int i=0; i<max_levels[j]; i++ ){
                        cache[j][i+1] = cache[j][i] * ( x[j] - nodes[i] );
                }
                for( int i=1; i<=max_levels[j]; i++ ){
                        cache[j][i] /= coeff[i];
                }
        }
        return cache;
}
double* GridSequence::cacheBasisIntegrals() const{
        int max_level = max_levels[0];  for( int j=1; j<num_dimensions; j++ ) if ( max_level < max_levels[j] ) max_level = max_levels[j];

        double *integ = new double[++max_level]; // integrals of basis functions
        std::fill( integ, integ + max_level, 0.0 );

        OneDimensionalNodes core;
        int n = 1 + max_level / 2; // number of Gauss-Legendre points needed to integrate the basis functions
        double *lag_x = 0, *lag_w = 0;
        core.getGaussLegendre( n, lag_w, lag_x );

        for( int i=0; i<n; i++ ){
                double v = 1.0;
                for( int j=1; j<max_level; j++ ){
                        v *= ( lag_x[i] - nodes[j-1] );
                        integ[j] += lag_w[i] * v / coeff[j];
                }
        }
        integ[0] = 2.0;

        delete[] lag_w;
        delete[] lag_x;

        return integ;
}

double GridSequence::evalBasis( const int f[], const int p[] ) const{
        double v = 1.0;
        for( int j=0; j<num_dimensions; j++ ){
                double x = nodes[ p[j] ];
                double w = 1.0;
                for( int i=0; i<f[j]; i++ ){
                        w *= ( x - nodes[i] );
                }
                v *= w / coeff[ f[j] ];
        }
        return v;
}

void GridSequence::recomputeSurpluses(){
        int n = points->getNumIndexes();
        if ( surpluses != 0 ) delete[] surpluses;
        surpluses = new double[ n * num_outputs ];

        if ( num_outputs > 2 ){
                #pragma omp parallel for schedule(static)
                for( int i=0; i<n; i++ ){
                        const double* v = values->getValues( i );
                        std::copy( v, v + num_outputs, &( surpluses[i*num_outputs] ) );
                }
        }else{
                const double* v = values->getValues( 0 );
                std::copy( v, v + n * num_outputs, surpluses );
        }

        IndexManipulator IM(num_dimensions);
        int *level = IM.computeLevels( points );
        int top_level = level[0];  for( int i=1; i<n; i++ ){  if ( top_level < level[i] ) top_level = level[i];  }

        for( int l=1; l<=top_level; l++ ){
                #pragma omp parallel for schedule(dynamic)
                for( int i=0; i<n; i++ ){
                        if ( level[i] == l ){
                                const int* p = points->getIndex( i );

                                int *monkey_count = new int[top_level + 1 ];
                                int *monkey_tail = new int[top_level + 1 ];
                                bool *used = new bool[n];
                                std::fill( used, used + n, false );

                                int current = 0;

                                monkey_count[0] = 0;
                                monkey_tail[0] = i;

                                while( monkey_count[0] < num_dimensions ){
                                        if ( monkey_count[current] < num_dimensions ){
                                                int branch = parents[ monkey_tail[current] * num_dimensions + monkey_count[current] ];
                                                if ( (branch == -1) || (used[branch]) ){
                                                        monkey_count[current]++;
                                                }else{
                                                        double basis_value = evalBasis( points->getIndex(branch), p );
                                                        for( int k=0; k<num_outputs; k++ ){
                                                                surpluses[ i*num_outputs + k ] -= basis_value * surpluses[ branch * num_outputs + k ];
                                                        }
                                                        used[branch] = true;

                                                        monkey_count[++current] = 0;
                                                        monkey_tail[current] = branch;
                                                }
                                        }else{
                                                monkey_count[--current]++;
                                        }
                                }

                                delete[] used;
                                delete[] monkey_tail;
                                delete[] monkey_count;
                        }
                }
        }

        delete[] level;
}

void GridSequence::applyTransformationTransposed( double weights[] ) const{
        IndexSet *work = ( points == 0 ) ? needed : points;
        int n = work->getNumIndexes();

        IndexManipulator IM(num_dimensions);
        int *level = IM.computeLevels( work );

        int top_level = level[0];  for( int i=1; i<n; i++ ){  if ( top_level < level[i] ) top_level = level[i];  }

        int *monkey_count = new int[top_level + 1];
        int *monkey_tail = new int[top_level + 1];
        bool *used = new bool[n];

        for( int l=top_level; l>0; l-- ){
                for( int i=0; i<n; i++ ){
                        if ( level[i] == l ){
                                const int* p = work->getIndex( i );
                                int current = 0;

                                monkey_count[0] = 0;
                                monkey_tail[0] = i;
                                std::fill( used, used + n, false );

                                while( monkey_count[0] < num_dimensions ){
                                        if ( monkey_count[current] < num_dimensions ){
                                                int branch = parents[  monkey_tail[current] * num_dimensions + monkey_count[current]  ];
                                                if ( (branch == -1) || used[branch] ){
                                                        monkey_count[current]++;
                                                }else{
                                                        weights[branch] -= weights[i] * evalBasis( work->getIndex( branch ), p );
                                                        used[branch] = true;

                                                        monkey_count[++current] = 0;
                                                        monkey_tail[current] = branch;
                                                }
                                        }else{
                                                monkey_count[--current]++;
                                        }
                                }
                        }
                }
        }

        delete[] level;
        delete[] used;
        delete[] monkey_tail;
        delete[] monkey_count;
}

}

#endif

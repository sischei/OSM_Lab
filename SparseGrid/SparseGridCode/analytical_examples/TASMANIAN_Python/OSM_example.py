#!/usr/bin/python

##############################################################################################################################################################################
# Copyright (c) 2017, Miroslav Stoyanov
#
# This file is part of
# Toolkit for Adaptive Stochastic Modeling And Non-Intrusive ApproximatioN: TASMANIAN
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
#    and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
#    or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
# OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
# OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# UT-BATTELLE, LLC AND THE UNITED STATES GOVERNMENT MAKE NO REPRESENTATIONS AND DISCLAIM ALL WARRANTIES, BOTH EXPRESSED AND IMPLIED.
# THERE ARE NO EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT,
# COPYRIGHT, TRADEMARK, OR OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL ACCOMPLISH THE INTENDED RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE.
# THE USER ASSUMES RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS, CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING FROM OR ARISING OUT OF,
# IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE SOFTWARE.
##############################################################################################################################################################################

# necessary import for every use of TASMANIAN
#
import TasmanianSG
import numpy as np

# imports specifically needed by the examples
import math
from random import uniform
from datetime import datetime

print("TasmanianSG version: {0:s}".format(TasmanianSG.__version__))
print("TasmanianSG license: {0:s}".format(TasmanianSG.__license__))

grid  = TasmanianSG.TasmanianSparseGrid()
grid1 = TasmanianSG.TasmanianSparseGrid()
grid2 = TasmanianSG.TasmanianSparseGrid()

#############################################################################

# EXAMPLE 1 for OSM:
# interpolate: f(x,y) = cos(0.5 * pi * x) * cos(0.5 * pi * y)
# using localp and semilocalp grids

# 1000 2-dimensional sample points 
aPnts = np.empty([1000, 2])  
for iI in range(1000):
    for iJ in range(2):
        aPnts[iI][iJ] = uniform(-1.0, 1.0)
#aPnts = aPnts[:,:2]

aTres = np.empty([1000,])
for iI in range(1000):
    aTres[iI] = math.cos(0.5 * math.pi * aPnts[iI][0]) * math.cos(0.5 * math.pi * aPnts[iI][1])

iDim = 2
iOut = 1
iDepth = 7

print("\n-------------------------------------------------------------------------------------------------")
print("Example 1 for OSM: interpolate f(x,y) = cos(0.5 * pi * x) * cos(0.5 * pi * y)")
print("       using localp and localp-zero rules with depth {0:1d}".format(iDepth))
print("       the error is estimated as the maximum from 1000 random points\n")

grid.makeLocalPolynomialGrid(iDim, iOut, iDepth, 2, "localp")

aPoints = grid.getPoints()
iNumP1 = aPoints.shape[0]
aVals = np.empty([aPoints.shape[0], 1])
for iI in range(aPoints.shape[0]):
    aVals[iI] = math.cos(0.5 * math.pi * aPoints[iI][0]) * math.cos(0.5 * math.pi * aPoints[iI][1])
grid.loadNeededPoints(aVals)

aRes = grid.evaluateBatch(aPnts)
fError1 = max(np.fabs(aRes[:,0] - aTres))

grid.makeLocalPolynomialGrid(iDim, iOut, iDepth-1, 2, "localp-zero")

aPoints = grid.getPoints()
iNumP2 = aPoints.shape[0]
aVals = np.empty([aPoints.shape[0], 1])
for iI in range(aPoints.shape[0]):
    aVals[iI] = math.cos(0.5 * math.pi * aPoints[iI][0]) * math.cos(0.5 * math.pi * aPoints[iI][1])
grid.loadNeededPoints(aVals)

aRes = grid.evaluateBatch(aPnts)
fError2 = max(np.fabs(aRes[:,0] - aTres))

print(" For localp    Number of points: {0:1d}   Error: {1:1.16e}".format(iNumP1, fError1))
print(" For localp-zero   Number of points: {0:1d}   Error: {1:1.16e}".format(iNumP2, fError2))
print(" Note: localp-zero wins this competition because the function is zero at the boundary")

#############################################################################


# EXAMPLE 2 for OSM:
# interpolate: f(x,y) = exp(-x) / (1 + 100 * exp(-10 * y))
# using different refinement schemes

aTres = np.empty([1000,])
for iI in range(1000):
    aTres[iI] = math.exp(-aPnts[iI][0]) / (1.0 + 100.0 * math.exp(-10.0 * aPnts[iI][1]))

iDim = 2
iOut = 1
iDepth = 2
fTol = 1.E-5

grid1.makeLocalPolynomialGrid(iDim, iOut, iDepth, -1, "localp")

aPoints = grid1.getPoints()
aVals = np.empty([aPoints.shape[0], 1])
for iI in range(aPoints.shape[0]):
    aVals[iI] = math.exp(-aPoints[iI][0]) / (1.0 + 100.0 * math.exp(-10.0 * aPoints[iI][1]))
grid1.loadNeededPoints(aVals)

grid2.makeLocalPolynomialGrid(iDim, iOut, iDepth, -1, "localp")
grid2.loadNeededPoints(aVals)

print("\n-------------------------------------------------------------------------------------------------")
print("Example 2: interpolate f(x,y) = exp(-x) / (1 + 100 * exp(-10 * y))")
print("   the error is estimated as the maximum from 1000 random points")
print("   tolerance is set at 1.E-5 and maximal order polynomials are used\n")

print("               Classic         FDS")
print(" precision    points     error    points     error")

for iK in range(7):
    grid1.setSurplusRefinement(fTol, -1, "classic")
    aPoints = grid1.getNeededPoints()
    aVals = np.empty([aPoints.shape[0], 1])
    for iI in range(aPoints.shape[0]):
        aVals[iI] = math.exp(-aPoints[iI][0]) / (1.0 + 100.0 * math.exp(-10.0 * aPoints[iI][1]))
    grid1.loadNeededPoints(aVals)

    aRes = grid1.evaluateBatch(aPnts)
    fError1 = max(np.fabs(aRes[:,0] - aTres))

    grid2.setSurplusRefinement(fTol, -1, "fds")
    aPoints = grid2.getNeededPoints()
    aVals = np.empty([aPoints.shape[0], 1])
    for iI in range(aPoints.shape[0]):
        aVals[iI] = math.exp(-aPoints[iI][0]) / (1.0 + 100.0 * math.exp(-10.0 * aPoints[iI][1]))
    grid2.loadNeededPoints(aVals)

    aRes = grid2.evaluateBatch(aPnts)
    fError2 = max(np.fabs(aRes[:,0] - aTres))

    print(" {0:9d} {1:9d}  {2:1.2e} {3:9d}  {4:1.2e}".format(iK+1, grid1.getNumPoints(), fError1, grid2.getNumPoints(), fError2))
  

#############################################################################
#
# EXAMPLE 3 for OSM:
# interpolate: f(x,y) = exp(-x) / (1 + 100 * exp(-10 * y))
# using local polynomails and wavelets

grid1.makeLocalPolynomialGrid(iDim, iOut, iDepth=3, iOrder=1, sRule="localp")
aPoints = grid1.getPoints()
aVals = np.empty([aPoints.shape[0], 1])
for iI in range(aPoints.shape[0]):
    aVals[iI] = math.exp(-aPoints[iI][0]) / (1.0 + 100.0 * math.exp(-10.0 * aPoints[iI][1]))
grid1.loadNeededPoints(aVals)

grid2.makeWaveletGrid(iDim, iOut, iDepth=1, iOrder=1)
aPoints = grid2.getPoints()
aVals = np.empty([aPoints.shape[0], 1])
for iI in range(aPoints.shape[0]):
    aVals[iI] = math.exp(-aPoints[iI][0]) / (1.0 + 100.0 * math.exp(-10.0 * aPoints[iI][1]))
grid2.loadNeededPoints(aVals)

print("\n-------------------------------------------------------------------------------------------------")
print("Example 3: interpolate f(x,y) = exp(-x) / (1 + 100 * exp(-10 * y))")
print("   the error is estimated as the maximum from 1000 random points")
print("   using local polynomials and wavelets\n")

print("           Polynomials        Wavelets")
print(" precision    points     error    points     error")

for iK in range(8):
    grid1.setSurplusRefinement(fTol, -1, "fds")
    aPoints = grid1.getNeededPoints()
    aVals = np.empty([aPoints.shape[0], 1])
    for iI in range(aPoints.shape[0]):
        aVals[iI] = math.exp(-aPoints[iI][0]) / (1.0 + 100.0 * math.exp(-10.0 * aPoints[iI][1]))
    grid1.loadNeededPoints(aVals)

    aRes = grid1.evaluateBatch(aPnts)
    fError1 = max(np.fabs(aRes[:,0] - aTres))

    grid2.setSurplusRefinement(fTol, -1, "fds")
    aPoints = grid2.getNeededPoints()
    aVals = np.empty([aPoints.shape[0], 1])
    for iI in range(aPoints.shape[0]):
        aVals[iI] = math.exp(-aPoints[iI][0]) / (1.0 + 100.0 * math.exp(-10.0 * aPoints[iI][1]))
    grid2.loadNeededPoints(aVals)

    aRes = grid2.evaluateBatch(aPnts)
    fError2 = max(np.fabs(aRes[:,0] - aTres))

    print(" {0:9d} {1:9d}  {2:1.2e} {3:9d}  {4:1.2e}".format(iK+1, grid1.getNumPoints(), fError1, grid2.getNumPoints(), fError2))

print("\n-------------------------------------------------------------------------------------------------\n")    


    
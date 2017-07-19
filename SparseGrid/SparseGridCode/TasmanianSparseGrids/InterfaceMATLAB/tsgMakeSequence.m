function  [ lGrid, points ] = tsgMakeSequence( sGridName, iDim, iOut, s1D, sType, iDepth, mTransformAB, vAnisotropy )
%
% [ lGrid, points ] = tsgMakeSequence( sGridName, iDim, iOut, s1D, sType, iDepth, mTransformAB, vAnisotropy )
%
% creates a new sparse grid using a sequence rule
%
% INPUT:
%
% sGridName: the name of the grid, give it a string name, 
%            i.e. 'myGrid' or '1' or 'pi314'
%            DO NOT LEAVE THIS EMPTY
%
% iDim: (integer, positive)
%       the number of inputs
%
% iOut: (integer, non-negative)
%       the number of outputs
%
% s1D: (string for the underlying 1-D rule that induces the grid)
%
%      'leja'           'rleja'          'rleja-shifted'   
%      'max-lebesgue'   'min-lebesgue'   'min-delta'
%
% sType: (string giving the tensor selection strategy)
%       'level'       'curved'         'tensor'         'iptensor'
%       'iptotal'     'ipcurved'       'qptotal'        'qpcurved'
%       'hyperbolic'  'iphyperbolic'   'qphyperbolic'
%
% iDepth: (integer non-negative)
%       controls the density of the grid, i.e., the offset for the tensor
%       selection, the meaning of iDepth depends on sType
%       Example 1: sType == 'iptotal' will give a grid that interpolates
%              exactly all polynomials of degree up to and including iDepth
%       Example 2: sType == 'qptotal' will give a grid that integrates
%              exactly all polynomials of degree up to and including iDepth
%
% vAnisotropy: (optional vector of positive integers, length iDim or 2*iDim)
%       the anisotropic weights associated with sType
%
% mTransformAB: (optional matrix of size iDim x 2)
%               for all but gauss-laguerre and gauss-hermite grids, the
%               transform specifies the lower and upper bound of the domain
%               in each direction. For gauss-laguerre and gauss-hermite
%               grids, the transform gives the a and b parameters that
%               change the weight to 
%               exp( -b ( x - a ) )  and  exp( -b ( x - a )^2 )
%
% OUTPUT:
%
% lGrid: list containing information about the sparse grid, can be used to
%        call other functions
%
% points: (optional) the points of the grid in an array of dimension [ num_poits, dim ]
%
% [ lGrid, points ] = tsgMakeSequence( sGridName, iDim, iOut, s1D, sType, iDepth, mTransformAB, vAnisotropy )
%

% create lGrid object
lGrid.sName = sGridName;
lGrid.iDim  = iDim;
lGrid.iOut  =  iOut;
lGrid.sType = 'sequence';

% check for conflict with tsgMakeQuadrature
if ( strcmp( sGridName, '' ) )
    error('sGridName cannot be empty');
end

% generate filenames
[ sFiles, sTasGrid ] = tsgGetPaths();
[ sFileG, sFileX, sFileV, sFileO, sFileW, sFileC ] = tsgMakeFilenames( lGrid.sName );

sCommand = [sTasGrid,' -makesequence'];

sCommand = [ sCommand, ' -gridfile ',   sFileG];
sCommand = [ sCommand, ' -dimensions ', num2str(lGrid.iDim)];
sCommand = [ sCommand, ' -outputs ',    num2str(lGrid.iOut)];
sCommand = [ sCommand, ' -onedim ',     s1D];
sCommand = [ sCommand, ' -depth ',      num2str(iDepth)];
sCommand = [ sCommand, ' -type ',       sType];

% set the domain transformation
if ( exist('mTransformAB') && (max( size( mTransformAB )) ~= 0) )
    if ( size( mTransformAB, 2 ) ~= 2 )
        error(' mTransformAB must be a matrix with 2 columns');
    end
    if ( size( mTransformAB, 1 ) ~= lGrid.iDim )
        error(' mTransformAB must be a matrix with iDim number of rows');
    end
    tsgWriteMatrix( sFileV, mTransformAB );
    lClean.sFileV = 1;
    sCommand = [ sCommand, ' -tf ',sFileV];
end

% set anisotropy
if ( exist('vAnisotropy') && (max( size( vAnisotropy )) ~= 0) )
    if ( min( size( vAlphaBeta ) ) ~= 1 )
        error(' vAnisotropy must be a vector, i.e., one row or one column');
    end
    if ( max( size( vAlphaBeta ) ) ~= lGrid.iDim )
        error(' vAnisotropy must be a vector of size iDim');
    end
    if ( size( vAnisotropy, 1 ) > size( vAnisotropy, 2 ) )
        tsgWriteMatrix( sFileW, vAnisotropy );
    else
        tsgWriteMatrix( sFileW, vAnisotropy' );
    end
    lClean.sFileW = 1;
    sCommand = [ sCommand, ' -anisotropyfile ',sFileW];
end

% read the points for the grid
if ( nargout > 1 )
    sCommand = [ sCommand, ' -of ',sFileO];
    lClean.sFileO = 1;
end

[status, cmdout] = system(sCommand);

if ( max( size( findstr( 'ERROR', cmdout ) ) ) ~= 0 )
    disp(cmdout);
    error('The tasgrid execurable returned an error, see above');
    return;
else
    if ( ~isempty(cmdout) )
        fprintf(1,['WARNING: Command had non-empty output:\n']);
        disp(cmdout);
    end
    if ( nargout > 1 )
        points = tsgReadMatrix( sFileO );
    end
end

if ( exist( 'lClean' ) )
    tsgCleanTempFiles( lGrid, lClean );
end

end

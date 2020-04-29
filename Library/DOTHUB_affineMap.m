function [A,B] = DOTHUB_affineMap(pointsFrom,pointsTo)

%  This function calculates an affine transform (A matrix and B vector) to
%  from the points defined by pointsFrom to those of pointsTO

% ####################### INPUTS ##########################################

% pointsFrom    :   [N x 3] matrix of coordinates in original space

% pointsTo      :   [N x 3] matrix of coordinates in target space

% ####################### OUTPUTS #########################################

% A             :   [3 x 3] Affine transformation matrix
% B             :   [3 x 1] Translation vector.

% preproFileName :  The full path of the resulting .prepro file.

% ####################### Dependencies ####################################

% #########################################################################
% This code was copied almost verbatim from Qianqian Fang <fangq at nmr.mgh.harvard.edu>
% The original function is part of "metch" toobox, see COPYING for license
% at https://github.com/fangq/metch/
% RJC, UCL, April 2020
%
% ############################# Updates ###################################
% #########################################################################


% parameters: 
%      pfrom: nx3 matrix, each row is a 3d point in original space
%      pto: nx3 matrix, each row is a 3d point in the mapped space
%
% outputs:
%      A: 3x3 matrix, the calculated affine A matrix
%      b: 3x1 vector, the calculated affine b vector
%
% the solution will satisfy the following equation: A*pfrom'+b=pto
%  author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
%  date: 12/12/2008

% Please find more information at http://iso2mesh.sf.net/cgi-bin/index.cgi?metch
%
% this function is part of "metch" toobox, see COPYING for license

% Please find more information at http://iso2mesh.sf.net/cgi-bin/index.cgi?metch
%
% this function is part of "metch" toobox, see COPYING for license

bsubmat=eye(3);
nPoints=size(pointsFrom,1);
if size(pointsTo,1)~=nPoints
    error('Both inputs should have the same dimensions [N x 3]');
end

amat = zeros(nPoints*3,9);

for i = 1:nPoints
    amat(i*3-2:i*3,:) = kron(bsubmat, pointsFrom(i,:));
end
amat = [amat,repmat(bsubmat,nPoints,1)];

bvec = pointsTo';
bvec = bvec(:);

tmp = amat\bvec;
A = reshape(tmp(1:9),3,3)';
B = tmp(end-2:end);

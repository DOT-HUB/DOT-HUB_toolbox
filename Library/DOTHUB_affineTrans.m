function outputPoints = DOTHUB_affineTrans(inputPoints,A,B)

%Performs 3D transformation of nx3 'input' using affine transformation matrix
%formed from A (3X3) and B (3x1)

%RJC UCL, Jan 2020 #########################################################

size_in = size(inputPoints);
if size_in(2) ~= 3
    error('Input points must be of dimensions nx3...');
    return
end

outputPoints = ((A*inputPoints') + repmat(B,1,size_in(1)))';


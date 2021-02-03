function range = DOTHUB_range(x,dim)
% Function to calculate the range of an array

% x: n dimensional array
% dim: optional input specifying array dimension

% ZK, Gowerlabs, February 2021

if nargin == 1 
    range = max(x) - min(x);
else
    range = max(x,[],dim) - min(x,[],dim);
end

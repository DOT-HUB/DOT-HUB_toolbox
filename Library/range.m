function [y] = range(A)

% This function is a toolbox-independent replacement for range from
% Statistics and Machine Learning toolbox
%
% ########## INPUTS ##########
% A             : array (vector or multidimensional, see doc min)
%
% ########## OUTPUTS ##########
% y             : a scalar of largest element minus smallet
%                 (for vector inputs) or a vector (for  
%                 multidimensional input)

y = max(A) - min(A);
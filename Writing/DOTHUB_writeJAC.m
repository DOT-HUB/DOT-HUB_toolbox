function [jac, jacFileName] = DOTHUB_writeJAC(jacFileName,logData,J,basis)

% This script creates a .jac file, which organizes and stores Jacobians for
% DOT reconstruction.

% ####################### INPUTS ##########################################

% jacFileName       :  The desired path &/ filename for the .jac file.
%                      This can be anything, but we recommend this variable be defined with the
%                      following code snippet, where: rmapFileName = full path and name
%                      of rmap on which jacobian calculated; transportPackage = 'toast++' or
%                      equivalent. opticalPropertiesByTissue = matrix of optical properties used.
%                      This snippet also provides recommended input variable 'logData'. 

                        % ds = datestr(now,'yyyymmDDHHMMSS');
                        % [pathstr, name, ~] = fileparts(rmapFileName);
                        % jacFileName = fullfile(pathstr,[name '.jac']);
                        % transportPackage = 'toast';
                        % logData(1,:) = {'Created on: ', ds};
                        % logData(2,:) = {'Derived from rmap file: ', rmapFileName};
                        % logData(3,:) = {'Calculated using: ', transportPackage};
                        % logData(4,:) = {'Optical properties (tissueInd, wavelength,[mua musPrime refInd]): ', opticalPropertiesByTissue};

% logData           :  (Optional). logData is a cell array of strings containing useful
%                      info as per snippet above. Parse empty to ignore.
                        
% J                 :   A cell structure of length nWavs, containing .vol and .gm 
%                       of dimensions #channels x #nodes
%                       Units are mm:
%                       d(ln(Intensity_active/Intensity_baseline))/d(absorbtion coefficient (mm-1))
%                       Note that the spatial dimension in the .vol variable can either be #nodes of
%                       the corresponding volume mesh in rmapFileName or the
%                       number of elements in the basis. The J{i}.gm is for 
%                       visualization, masking and cortex-constrained recon purposes.
%                       It contains the GM extraction of the volume Jacobian. 
%                       Dimensions of channel x #gm_nodes

% basis             :   (Optional) 1x3 vector specifying the basis in which J is defined.
%                       If basis is not parsed, jac.basis will be saved empty = [] 
%                       and J{i}.vol is assumed to be in the volume mesh space

% ####################### OUTPUTS #########################################

% jac                      :  Structure containing all data inputs

% jacFileName              :  The full path of the resulting .jac file

% jacFileName.jac file     :  File containing all data inputs

% ####################### Dependencies ####################################
% #########################################################################
% RJC, UCL, April 2020
%
% ############################# Updates ###################################
% #########################################################################

% MANAGE VARIABLES
% #########################################################################

jac.J = J;

if isempty(logData)
    logData = {};
    warning('logData is empty: this might make it harder to keep track of your data...');
end
jac.logData = logData;

if ~exist('basis','var')
    jac.basis = [];
else
    jac.basis = basis;
end

%Create filename ##########################################################
[pathstr, name, ext] = fileparts(jacFileName);
if isempty(ext) || ~strcmpi(ext,'.jac')
    ext = '.jac';
end
if isempty(pathstr)
    pathstr = pwd;
end
jacFileName = fullfile(pathstr,[name ext]);
jac.fileName = jacFileName; %including the fileName within the structure is very useful 
%for tracking and naming things derived further downstream.

if exist(jacFileName,'file')
    warning([name ext ' will be overwritten...']);
end

%Save .jac file ###########################################################
save(jacFileName,'-struct','jac');
fprintf('##################### Writing .jac file #########################\n');
fprintf(['.jac data file saved as ' jacFileName '\n']);
fprintf('\n');

function [invjac, invjacFileName] = DOTHUB_writeINVJAC(invjacFileName,logData,invJ,basis,saveFlag)

% This script creates a .invjac file, which organizes and stores inverted Jacobian data for
% DOT reconstruction.

% ####################### INPUTS ##########################################

% invjacFileName       :  The desired path &/ filename for the .invjac file.
%                      This can be anything, but we recommend this variable be defined with the
%                      following code snippet.
                        
%                      [pathstr, name, ~] = fileparts(jac.fileName);
%                      invjacFileName = fullfile(pathstr,[name '.invjac']);
%                      invjac.fileName = invjacFileName;
%                      ds = datestr(now,'yyyymmDDHHMMSS');
%                      logData(1,:) = {'Created on: ', ds};
%                      logData(2,:) = {'Derived from jac file: ', jac.fileName};
%                      logData(3,:) = {'reconMethod: ', varInputs.reconMethod};
%                      logData(4,:) = {'regMethod: ', varInputs.regMethod};
%                      logData(5,:) = {'hyperParameter: ', varInputs.regMethod}; 

% logData           :  (Optional). logData is a cell array of strings containing useful
%                      info as per snippet above. Parse empty to ignore.
                        
% invJ              :   A cell array, with each cell containing a
%                       wavelength-specific jacobian of dimensions channel x node. 
%                       If the reconMethod was multispectral, invJ only has one entry 
%                       - the inverted multispectral jacobian.

% basis             :   (Optional) 1x3 vector specifying the basis in which J is defined.
%                       If no basis is used, parse empty = [] and 
%                       space is assumed to be in the volume mesh space

% saveFlag          :   (Optional) 1 or 0 to indicate if invjac should be saved to
%                       disk. Defaults to

% ####################### OUTPUTS #########################################

% invjac             :  Structure containing all data inputs

% invjacFileName     :  The full path of the resulting .jac file

% invjacFileName.invjac file  :  File containing all data inputs (if saved)

% ####################### Dependencies ####################################
% #########################################################################
% RJC, UCL, April 2020
%
% ############################# Updates ###################################
% #########################################################################

% MANAGE VARIABLES
% #########################################################################

invjac.invJ = invJ;
invjac.basis = basis;
if isempty(logData)
    warning('logData is empty: this might make it harder to keep track of your data...');
end
invjac.logData = logData;


%Create filename ##########################################################
[pathstr, name, ext] = fileparts(invjacFileName);
if isempty(ext) || ~strcmpi(ext,'.invjac')
    ext = '.invjac';
end
if isempty(pathstr)
    pathstr = pwd;
end
invjacFileName = fullfile(pathstr,[name ext]);
invjac.fileName = invjacFileName; %including the fileName within the structure is very useful 
%for tracking and naming things derived further downstream.

if saveFlag
    %Save .invjac file ###########################################################
    save(invjacFileName,'-struct','invjac');
    fprintf(['.invjac data file saved as ' invjacFileName '\n']);
end

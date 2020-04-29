function [prepro, preproFileName] = DOTHUB_writePREPRO(preproFileName,logData,SD3D,tDOD,dod,tHRF,dcAvg,dcAvgStd)

% This script creates a .prepro file, containing pre-processed HD-DOT data. 
% The .pre file is primarily designed to support image reconstruction, 
% but in cases where the HRF is calculated in channel space, it also 
% provides a space to save HRF-related variables that might not be needed
% for reconstruction but for channel plotting (if appropriate). The dod data
% in a .pre file will be used for reconstruction. It can be the averaged HRF 
% data, or a full time course of data.  Structures mirror those of Homer2. 
% The .pre file is designed to be the output of the project/study specific 
% pre-processing pipeline, the name of which is specified by streamMname.

% ####################### INPUTS ##########################################

% preproFileName    :  The desired path &/ filename for the .prepro file.
%                      Can be anything. Can be parsed empty to prevent saving 
%                      but we recommend this variable be defined with the
%                      following snippet, where 'nirsFileName' is the full
%                      path and name of the .nirs file being pre-processed.
%                      This snippet also applies for input variable
%                      'logData'.
                       
                        %[pathstr, name, ~] = fileparts(nirsFileName);
                        %ds = datestr(now,'yyyymmDDHHMMSS');
                        %preproFileName = fullfile(pathstr,[name '.prepro']);
                        %logData(1,:) = {'Created on: '; ds};
                        %logData(2,:) = {'Derived from data: ', nirsFileName};
                        %logData(3,:) = {'Pre-processed using:', mfilename('fullpath')};
                         
% logData           :   (Optional). logData is a cell array of strings containing useful
%                       info as per snippet above. Parse empty to ignore.

% SD3D             :   The source-detector structure (Homer2 style) with 3D
%                       optode positions

% tDOD              :   The time vector in seconds corresponding to the first
%                       dimension of dod.

% dod               :   Change in absorbance (natural log definition).  As per
%                       Homer2 definition (time x channel) if full timecourse,
%                       or time x channel x condition for multi-condition HRF

% tHRF              :   (Optional) The time vector in seconds corresponding to 
%                       the HRF (this can be the same as tDOD if only saving
%                       HRF data or can be different if saving both HRF and
%                       full time course.

% dcAvg             :   (Optional) Change in concentration HRF data.  As per Homer2
%                       definition (time x chromophore x channel) for one condition,
%                       or time x channel x condition for multi-condition HRF

% dcAvgStd          :   (Optional) Std of change in concentration HRF.  As per Homer2
%                       definition (time x chromophore x channel) for one condition,
%                       or time x channel x condition for multi-condition HRF

% ####################### OUTPUTS #########################################

% prepro            :  The prepro structure

% preproFileName    :  The full path of the resulting .prepro file. Just
%                      returns subjectName if savePath is not saved

% .prepro           :  A file containing a structure of all the data inputs parsed

% ####################### Dependencies ####################################
% #########################################################################
% RJC, UCL, April 2020
%
% ############################# Updates ###################################
% #########################################################################

% MANAGE VARIABLES
% #########################################################################

prepro.SD3D = SD3D; %These are must-haves
prepro.tDOD = tDOD;
prepro.dod = dod;

if isempty(logData)
    warning('logData is empty: this might make it harder to keep track of your data. See comments for more info.');
end
prepro.logData = logData; %This might be empty

%Add optional inputs:
if exist('tHRF','var')
    if ~isempty(tHRF)
        prepro.tHRF = tHRF;
    end
end

if exist('dcAvg','var')
    if ~isempty(dcAvg)
        prepro.dcAvg = dcAvg;
    end
end

if exist('dcStd','var')
    if ~isempty(dcAvgStd)
         prepro.dcStd = dcAvgStd;
    end
end

% Skip naming and saving if preproFileName parsed empty (streaming?)
if isempty(preproFileName); return; end

%Create filename ##########################################################
[pathstr, name, ext] = fileparts(preproFileName);
if isempty(ext) || ~strcmpi(ext,'.prepro')
    ext = '.prepro';
end
if isempty(pathstr)
    pathstr = pwd;
end
preproFileName = fullfile(pathstr,[name ext]);
prepro.fileName = preproFileName; %including the fileName within the structure is very useful 
%for tracking and naming things derived further downstream.

if exist(preproFileName,'file')
    warning([name ext ' will be overwritten...']);
end
%Save .pre file ###########################################################
save(preproFileName,'-struct','prepro');
fprintf('##################### Writing .prepro file #######################\n');
fprintf(['.prepro data file saved as ' preproFileName '\n']);
fprintf('\n');
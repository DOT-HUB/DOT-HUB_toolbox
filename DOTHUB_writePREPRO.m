function DOTHUB_LUMOwritepreproc(nirsFilename,streamMname,SD_3D,tDOD,dod,tHRF,dcAvg,dcStd)

% This script creates a .pre file, containing pre-processed HD-DOT data. 
% The .pre file is primarily designed to support image reconstruction, 
% but in cases where the HRF is calculated in channel space, it also 
% provides a space to save HRF-related variables that might not be needed
% for reconstruction but for channel plotting (if appropriate). The dod data
% in a .pre file will be used for reconstruction. It can be the averaged HRF 
% data, or a full time course of data.  Structures mirror those of Homer2. 
% The .pre file is designed to be the output of the project/study specific 
% pre-processing pipeline, the name of which is specified by streamMname.

% ####################### INPUTS ##########################################

% nirsFilename  :   The full path of the .nirs data file from which the
%                   pre-processing data was derived

% streamMname   :   Please parse mfilename('fullpath'). This gives the filename 
%                   of the pro-processing script used to create this data.  
%                   This is just for book-keeping purposes. 

% SD_3D         :   The source-detector structure (Homer2 style) with 3D
%                   optode positions

% tDOD          :   The time vector in seconds corresponding to the first
%                   dimension of dod.

% dod           :   Change in absorbance (natural log definition).  As per
%                   Homer2 definition (time x channel) if full timecourse,
%                   or time x channel x condition for multi-condition HRF

% tHRF           :  (Optional) The time vector in seconds corresponding to 
%                   the HRF (this can be the same as tDOD if only saving
%                   HRF data)

% dcAvg         :   (Optional) Change in concentration HRF data.  As per Homer2
%                   definition (time x chromophore x channel) for one condition,
%                   or time x channel x condition for multi-condition HRF

% dcAvgStd      :   (Optional) Std of change in concentration HRF.  As per Homer2
%                   definition (time x chromophore x channel) for one condition,
%                   or time x channel x condition for multi-condition HRF

% ####################### OUTPUTS #########################################

% .pre file containing all the data inputs parsed plus a logdata cell array 
% specifying the .nirs path, date and the processing stream name (streamMname)

% ####################### Dependencies ####################################

% #########################################################################
% RJC, UCL, April 2020
%
% ############################# Updates ###################################
% #########################################################################

% MANAGE VARIABLES
% #########################################################################

savelist = {'SD_3D','tDOD','dod','logdata'}; %These are required

if exist('tHRF','var')
    if ~isempty(tHRF)
        savelist{end+1} = 'tHRF';
    end
end

if exist('dcAvg','var')
    if ~isempty(dcAvg)
        savelist{end+1} = 'dcAvg';
    end
end

if exist('dcStd','var')
    if ~isempty(dcStd)
        savelist{end+1} = 'dcStd';
    end
end
    
%Create logdata ###########################################################
ds = datestr(now,'yyyymmDDHHMMSS');
logdata{1,1} = ['Created on: ' ds];
logdata{2,1} = ['Derived from file: ' nirsFilename];
logdata{3,1} = ['Pre-processed using: ' streamMname '.m'];

%Save .pre file ###########################################################
[pathstr, name, ext] = fileparts(nirsFilename);
outname = fullfile(pathstr,[name '_' ds '.pre']);
save(outname,savelist{:});
fprintf(['.pre data file saved as ' outname '\n']);


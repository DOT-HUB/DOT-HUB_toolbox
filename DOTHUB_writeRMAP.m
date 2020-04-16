function DOTHUB_LUMOwritermap(nirsFilename,HeadVolumeMesh_reg,GMSurfaceMesh_reg,ScalpSurfaceMesh_reg,vol2gm,srcPos,detPos)

% This script creates a registered mesh and positions (.rmap) file, 

% ####################### INPUTS ##########################################

% mapFilename   :   The full path of the .map data file from which the
%                   pre-processing data was derived

% ####################### OUTPUTS #########################################

% .jac file containing all

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


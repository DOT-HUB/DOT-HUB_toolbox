% DOT-HUB toolbox script to run five-wavelength data

% based off RJC 2020 UCL
% modified by GCVL 2022 Cambridge

%% Specify paths of pre-defined elements (.LUMO, atlas .mshs, Homer2 preprocessing .cfg file).

[filepath,~,~] = fileparts(mfilename('fullpath'));
if ~exist('LUMODirName','var')
    disp('Select .LUMO directory...');
    LUMODirName = uigetdir(pwd,'Select .LUMO directory');
elseif isempty(LUMODirName)
    disp('Select .LUMO directory...');
    LUMODirName = uigetdir(pwd,'Select .LUMO directory');
elseif ~exist(LUMODirName,'dir')
    disp('Specified directory not found, please select .LUMO directory...');
    LUMODirName = uigetdir(pwd,'Select .LUMO directory');
end
origMeshFileName = [filepath '/ExampleMeshes/AdultMNI152.mshs'];
%cfgFileName = [filepath '/ExampleData/Example1/preproPipelineExample1.cfg'];

%%
addpath('LUMO_toolbox');
%% Covert .LUMO to .nirs
[nirs, nirsFileName, SD3DFileName] = DOTHUB_LUMO2nirs_5wav(LUMODirName);
%s[nirs, nirsFileName, SD3DFileName] = DOTHUB_LUMO2nirs (LUMODirName);

%SD3D.DetPos(9:16) = [];
%% Run data quality checks - this produces multiple figures, so comment out for speed.
DOTHUB_dataQualityCheck(nirsFileName);
disp('Examine data quality figures, press any key to continue');
pause 



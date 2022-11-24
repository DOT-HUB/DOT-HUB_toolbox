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

%% Run data quality checks - this produces multiple figures, so comment out for speed.
DOTHUB_dataQualityCheck(nirsFileName);
disp('Examine data quality figures, press any key to continue');
pause 

%% Run Homer2 pre-processing pipeline using .cfg file. Alternatively you can run line by line (as per commented below).

%[prepro, preproFileName] = DOTHUB_runHomerPrepro(nirsFileName,cfgFileName);

%%%%%Equivalent line-by-line Homer2 calls and prepro write:
 dod = hmrIntensity2OD(nirs.d);
 SD3D = enPruneChannels(nirs.d,nirs.SD3D,ones(size(nirs.t)),[0 1e11],12,[0 100],0); 

 %Force MeasListAct to be the same across wavelengths
 SD3D = DOTHUB_balanceMeasListAct(SD3D);

%Set SD2D
 SD2D = nirs.SD; 
 SD2D.MeasListAct = SD3D.MeasListAct;
 
%Bandpass filter and convert to Concentration
 dod = hmrBandpassFilt(dod,nirs.t,0,0.5);
 dc = hmrOD2Conc(dod,SD3D,[6 6]);
 dc = dc*1e6; %Homer works in Molar by default, we use uMolar.

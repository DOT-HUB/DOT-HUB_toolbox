% DOT-HUB toolbox Pipeline Example 1.
%
% What follows is an example of a wrapper script that employs the main 
% steps of the toolbox. Most steps output variables into the work
% space and also writes them out as the key file types, so you can comment 
% out steps as you work through them and pick up where you left off, rather
% than re-running every step. The whole script runs in ~12 minutes on a 2018
% MacBook Pro with 16Gb RAM.
%
% Example 1 is the simplest application of the toolbox to LUMO data. It is
% assumed that we have no subject-specific information about the position
% of the optodes (and therefore 3D positioning information is derived from
% the default 3D values in  the LUMOcap .json file) nor do we have
% subject-specific structural MRI; we use an adult atlas.
%
% The dataset is an adult visual eccentricity experiment equivalent to that
% described in Vidal-Rosas et. al. 2021 Neurophotonics (in review):
% "Evaluating a new generation of wearable high-density diffuse optical
% tomography technology via retinotopic mapping of the adult visual cortex"
%
% RJC, UCL, Dec 2020.

%% Specify paths of pre-defined elements (.LUMO, atlas .mshs, Homer2 preprocessing .cfg file).
[filepath,~,~] = fileparts(mfilename('fullpath'));
LUMODirName = [filepath '/ExampleData/Example1/Example1_VisualCortexEccentricity.LUMO'];
origMeshFileName = [filepath '/ExampleMeshes/AdultMNI152.mshs'];
cfgFileName = [filepath '/ExampleData/Example1/preproPipelineExample1.cfg'];

%% Covert .LUMO to .nirs
[nirs, nirsFileName, SD3DFileName] = DOTHUB_LUMO2nirs(LUMODirName);

%% Run data quality checks - this produces multiple figures, so comment out for speed.
%DOTHUB_dataQualityCheck(nirsFileName);
%disp('Examine data quality figures, press any key to continue');
%pause 

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

%Regress short channels
dc = DOTHUB_hmrSSRegressionByChannel(dc,SD3D,12,1); %This is a custom SS regression script. 

%Block avg
[dcAvg,dcAvgStd,tHRF] = hmrBlockAvg(dc,nirs.s,nirs.t,[-5 25]);

%Convert back to dod for reconstruction
dodRecon = DOTHUB_hmrConc2OD(dcAvg/1e6,SD3D,[6 6]); %Note converting back to Molar units here for Homer function
tDOD = tHRF;

% Use code snippet from DOTHUB_writePREPRO to define contents of logs:
[pathstr, name, ~] = fileparts(nirsFileName);
ds = datestr(now,'yyyymmDDHHMMSS');
preproFileName = fullfile(pathstr,[name '.prepro']);
logData(1,:) = {'Created on: '; ds};
logData(2,:) = {'Derived from data: ', nirsFileName};
logData(3,:) = {'Pre-processed using:', mfilename('fullpath')};
[prepro, preproFileName] = DOTHUB_writePREPRO(preproFileName,logData,dodRecon,tDOD,SD3D,nirs.s,dcAvg,dcAvgStd,tHRF,nirs.CondNames,SD2D);

%% Plot prepro HRF results as array map if desired. Make sure you parse the 2D version of the array.
conditionToPlot = 1;
y = squeeze(prepro.dcAvg(:,:,:,conditionToPlot)); %Crop out chosen condition to plot
figure;
DOTHUB_LUMOplotArray(y,prepro.tHRF,prepro.SD2D);

%% Register chosen mesh to subject SD3D and create rmap
[rmap, rmapFileName] = DOTHUB_meshRegistration(nirsFileName,origMeshFileName);
figure;
DOTHUB_plotRMAP(rmap)

%% Calculate Jacobian 
basis = [30 30 30];
[jac, jacFileName] = DOTHUB_makeToastJacobian(rmapFileName,basis);

%% Invert Jacobian
%Note that you can either separately calculate the inverse, or just run
%DOTHUB_reconstruction, which will then call the inversion itself
[invjac, invjacFileName] = DOTHUB_invertJacobian(jacFileName,preproFileName,'saveFlag',true,'reconMethod','multispectral','hyperParameter',0.01);

%% Reconstruction
%Reconstruct
[dotimg, dotimgFileName] = DOTHUB_reconstruction(preproFileName,[],invjacFileName,rmapFileName,'saveVolumeImages',true);

%% Display peak response results on atlas surface and in volume
timeRange = [10 15]; %seconds post-onset
fs = length(dotimg.tImg)./DOTHUB_range(dotimg.tImg);
frameRange = round((timeRange + abs(min(dotimg.tImg))).*fs);
frames = frameRange(1):frameRange(2);
DOTHUB_plotSurfaceDOTIMG(dotimg,origMeshFileName,frames,'condition',3,'view',[0 20]);
DOTHUB_plotVolumeDOTIMG(dotimg,origMeshFileName,frames);




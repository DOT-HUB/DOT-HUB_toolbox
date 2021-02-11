% DOT-HUB toolbox Pipeline Example 2.
%
% What follows is an example of a wrapper script that employs the main 
% steps of the toolbox. Most steps output variables into the work
% space and also writes them out as the key file types, so you can comment 
% out steps as you work through them and pick up where you left off, rather
% than re-running every step. The whole script runs in ~12 minutes on a 2018
% MacBook Pro with 16Gb RAM.
%
% Example 2 shows the application of the toolbox to LUMO data when
% subject-specific information about the position of the optodes is
% available (in the form of a CSV) but we don't have subject-specific 
% structural MRI, so we use an adult atlas.
%
% The dataset is an adult visual eccentricity experiment equivalent to that
% described in Vidal-Rosas et. al. 2021(?) Neurophotonics (in review):
% "Evaluating a new generation of wearable high-density diffuse optical
% tomography technology via retinotopic mapping of the adult visual cortex"
%
% RJC, UCL, Dec 2020.
% ZK, Gowerlabs, Dec 2020.

%% Specify paths of pre-defined elements (.LUMO, atlas .mshs, Homer2 preprocessing .cfg file).
[filepath,~,~] = fileparts(mfilename('fullpath'));
LUMODirName = [filepath '/ExampleData/Example2/Example2_VisualCortexEccentricity.LUMO'];
origMeshFileName = [filepath '/ExampleMeshes/AdultMNI152.mshs'];
cfgFileName = [filepath '/ExampleData/Example2/preproPipelineExample2.cfg'];
posCSVFileName = [filepath '/ExampleData/Example2/Example2_Polhemus.csv'];

%% Covert .LUMO to .nirs - here is the only real difference with example pipeline 1 - we call the polhemus data.
[nirs, nirsFileName, SD3DFileName] = DOTHUB_LUMO2nirs(LUMODirName, [], posCSVFileName);

%% Run data quality checks - this produces multiple figures, so comment out for speed.
%DOTHUB_dataQualityCheck(nirsFileName);
%disp('Examine data quality figures, press any key to continue');
%pause 

%% Run Homer2 pre-processing pipeline 
[prepro, preproFileName] = DOTHUB_runHomerPrepro(nirsFileName,cfgFileName);

%% Plot prepro HRF results as array map if desired. Make sure you parse the 2D version of the array.
conditionToPlot = 1;
y = squeeze(prepro.dcAvg(:,:,:,conditionToPlot)); %Crop out chosen condition to plot
figure;
DOTHUB_LUMOplotArray(y,prepro.tHRF,prepro.SD2D);

%% Register chosen mesh to subject SD3D and create rmap
[rmap, rmapFileName] = DOTHUB_meshRegistration(nirsFileName,origMeshFileName);
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
DOTHUB_plotSurfaceDOTIMG(dotimg,rmap,frames,'condition',3,'view',[0 20]);
DOTHUB_plotVolumeDOTIMG(dotimg,rmap,frames,'slicePos',[-12 40 48]);




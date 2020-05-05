% DOTHUB Recon Pipeline Example 2.
%
% Individual optode positioning information, with Atlas head model
%
% This example uses a real HD-DOT dataset for a single subject from the micro-NTS (uNTS)
% dataset obtained in 2016 (Chitnis et al, 2016 https://doi.org/10.1364/BOE.7.004275).
% The pre-exisiting inputs can be found in ExampleData/uNTS_fingerTapping, and 
% consist of the raw .nirs file for our subject and tje separate .SD3D file. 
% The atlas mesh we will use for reconstruction (AdultMNI152.mshs) is in the
% meshes folder. In this example, each file associated with the pipeline 
% is saved to disk, and the subsequent function is then parsed the filename.  
% This is slower than parsing the structure directly, but means you can
% pause/comment out any step in turn without having to re-run everything
% (although the whole pipeline runs in only ~6 minutes on my Macbook.
%
% RJC, UCL, April 2020.

%Paths of pre-defined elements (.nirs file, atlas mesh, SD3D file, Homer2
%cfg file)
nirsFileName = 'ExampleData/Example2/uNTS_FingerTap_Subj01.nirs';
origMeshFileName = 'Meshes/AdultMNI152.mshs';
SD3DFileName = 'ExampleData/Example2/uNTS_FingerTap_Subj01.SD3D';
cfgFileName = 'ExampleData/Example2/preproPipelineExample2.cfg';

%Run data quality checks
DOTHUB_dataQualityCheck(nirsFileName,1)

%Run bespoke pre-processing script (simplest possible example included below)
[prepro, preproFileName] = DOTHUB_runHomerPrepro(nirsFileName,cfgFileName);

%Register chosen mesh to subject SD3D and create rmap
[rmap, rmapFileName] = DOTHUB_meshRegistration(nirsFileName,origMeshFileName);

%Calculate Jacobian (use small basis for speed of example)
[jac, jacFileName] = DOTHUB_makeToastJacobian(rmapFileName,[10 10 10]);

%You can either separately calculate the inverse, or just run
%DOTHUB_reconstruction, which will then call the inversion.
[invjac, invjacFileName] = DOTHUB_invertJacobian(jacFileName,preproFileName,'saveFlag',true,'reconMethod','multispectral','hyperParameter',0.01);

%Reconstruct
[dotimg, dotimgFileName] = DOTHUB_reconstruction(preproFileName,[],invjacFileName,rmapFileName,'saveVolumeImages',true);

%Display peak response results on surface and in volume
frames = 55:75;
DOTHUB_plotSurfaceDOTIMG(dotimg,rmap,frames,'view',[-130 15])
DOTHUB_plotVolumeDOTIMG(dotimg,rmap,frames);




% DOTHUB Recon Pipeline Example 1.
%
% This example uses a real HD-DOT dataset for a single subject from the micro-NTS (uNTS)
% dataset obtained in 2016 (Chitnis et al, 2016 https://doi.org/10.1364/BOE.7.004275).
% The pre-exisiting inputs can be found in ExampleData/uNTS_fingerTapping, and 
% consist of theraw .nirs file for our subject (already containing the SD3D).
% A separate .SD3D file, and the atlas mesh we will use for reconstruction
% (AdultMNI152.mshs). In this example, each file associated with the pipeline 
% is saved to disk, and the subsequent function is given the filename. This 
% is slower than parsing the structure directly, but means you can
% pause/comment out each step in turn without having to re-run everything
% (although the whole pipeline runs in ~6 minutes on my Macbook.
%
% RJC, UCL, April 2020.

nirsFileName = 'ExampleData/uNTS_fingerTapping/uNTS_FingerTap_Subj01.nirs';
origMeshFullFileName = 'ExampleData/uNTS_fingerTapping/AdultMNI152.mshs';
SD3DFileName = 'ExampleData/uNTS_fingerTapping/uNTS_FingerTap_Subj01.SD3D';

%Run bespoke pre-processing script (simplest possible example included below)
[prepro, preproFileName] = examplePreProcessingScript(nirsFileName);

%Register chosen mesh to subject SD3D and create rmap
[rmap, rmapFileName] = DOTHUB_meshRegistration(nirsFileName,origMeshFullFileName);

%Calculate Jacobian (use small basis for speed of example)
[jac, jacFileName] = DOTHUB_makeToastJacobian(rmap,[10 10 10]);

%You can either separately calculate the inverse, or just run
%DOTHUB_reconstruction, which will then call the inversion.
[invjac, invjacFileName] = DOTHUB_invertJacobian(jac,prepro,'saveFlag',true,'reconMethod','multispectral','hyperParameter',0.02);%,'regMethod','covariance');

%Reconstruct
[dotimg, dotimgFileName] = DOTHUB_reconstruction(preproFileName,[],invjacFileName,rmapFileName);

%Display results
rmap = load(rmapFileName,'-mat');
DOTHUB_displayOnMesh(rmap.gmSurfaceMesh,dotimg.hbo.gm(60,:))

%jacFileName = 'ExampleData/uNTS_fingerTapping/uNTS_FingerTap_Subj01.jac';
%preproFileName = 'ExampleData/uNTS_fingerTapping/uNTS_FingerTap_Subj01.prepro';
%rmapFileName = 'ExampleData/uNTS_fingerTapping/uNTS_FingerTap_Subj01.rmap';






function [prepro, preproFileName] = examplePreProcessingScript(nirsFileName)

load(nirsFileName,'-mat');

od = hmrIntensity2OD(d);
od = hmrBandpassFilt(od,1/mean(diff(t)),0.01,0.5);
dc = hmrOD2Conc(od,SD3D,[6 6]);
[dcAvg, dcStd, tHRF] = hmrBlockAvg(dc,s,t,[-5 25]);
dodAvg = DOTHUB_hmrHRFConc2OD(dcAvg, SD3D,[6 6]);

%USE CODE SNIPPET FROM DOTHUB_writePREPRO to define filename and logData
[pathstr, name, ~] = fileparts(nirsFileName);
ds = datestr(now,'yyyymmDDHHMMSS');
preproFileName = fullfile(pathstr,[name '.prepro']);
logData(1,:) = {'Created on: '; ds};
logData(2,:) = {'Derived from data: ', nirsFileName};
logData(3,:) = {'Pre-processed using:', mfilename('fullpath')};

%DOTHUB_writePREPRO(preproFileName,logData,SD3D,tDOD,dod,tHRF,dcAvg,dcStd)
[prepro, preproFileName] = DOTHUB_writePREPRO(preproFileName,logData,SD3D,tHRF,dodAvg,tHRF,dcAvg,dcStd);

end
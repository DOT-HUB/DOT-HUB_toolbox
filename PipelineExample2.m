%Run data quality assessment
%DOTHUB_dataQualityCheck(nirsFileName); Write this function? Just turn off plot labels.

%nirsFileName = 'ExampleData/uNTS_fingerTapping/uNTS_FingerTap_Subj01.nirs';

%Run bespoke pre-processing script (simplest possible example below)
%[prepro, preproFileName] = examplePreProcessingScript('ExampleData/uNTS_fingerTapping/uNTS_FingerTap_Subj01.nirs');

%Register chosen mesh to subject SD3D
%origMeshFullFileName = 'ExampleData/uNTS_fingerTapping/AdultMNI152.mshs';
%SD3DFileName = 'ExampleData/uNTS_fingerTapping/uNTS_FingerTap_Subj01.SD3D';
%[rmap, rmapFileName] = DOTHUB_meshRegistration(nirsFileName,origMeshFullFileName);

%Calculate Jacobian
%[jac, jacFileName] = DOTHUB_makeToastJacobian(rmap,[10 10 10]);

%You can either separately calculate the inverse, or just run DOTHUB_reconstruction. 
%[invjac, invjacFileName] = DOTHUB_invertJacobian(jacFileName,preproFileName,'saveFlag',true);

preproFileName = 'ExampleData/uNTS_fingerTapping/uNTS_FingerTap_Subj01.prepro';
jacFileName = 'ExampleData/uNTS_fingerTapping/uNTS_FingerTap_Subj01.jac';
rmapFileName = 'ExampleData/uNTS_fingerTapping/uNTS_FingerTap_Subj01.rmap';

%Reconstruct
[dot, dotFileName] = DOTHUB_reconstruction(preproFileName,jacFileName,[],rmapFileName);









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
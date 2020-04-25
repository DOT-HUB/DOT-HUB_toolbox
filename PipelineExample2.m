
%Run data quality assessment
%DOTHUB_dataQualityCheck(nirsFileName); Write this function? Just turn off plot labels.

%Run bespoke pre-processing script (simplest possible example below)
%[prepro, preproFileName] = examplePreProcessingScript('ExampleData/uNTS_fingerTapping/uNTS_FingerTap_Subj01.nirs');

%Register chosen mesh to subject SD_3D
%origMeshFullFileName = 'ExampleData/uNTS_fingerTapping/AdultMNI152.mshs';
%SD_3DFullFileName = 'ExampleData/uNTS_fingerTapping/uNTS_FingerTap_Subj01_3D.SD';
%[rmap, rmapFileName] = DOTHUB_meshRegistration(SD_3DFullFileName,origMeshFullFileName);

%Calculate Jacobian
rmapFileName = 'ExampleData/uNTS_fingerTapping/AdultMNI152_Reg2_uNTS_FingerTap_Subj01_3D.rmap';
[jac, jacFileName] = DOTHUB_makeToastJacobian(rmapFileName,[10 10 10]);

%JUST Reconstruction to go!
[dot, dotFileName] = DOTHUB_reconstruction();


function [prepro, preproFileName] = examplePreProcessingScript(nirsFileName)

load(nirsFileName,'-mat');

od = hmrIntensity2OD(d);
od = hmrBandpassFilt(od,1/mean(diff(t)),0.01,0.5);
dc = hmrOD2Conc(od,SD_3D,[6 6]);
[dcAvg, dcStd, tHRF] = hmrBlockAvg(dc,s,t,[-5 25]);
dodAvg = DOTHUB_hmrHRFConc2OD(dcAvg, SD_3D,[6 6]);

%USE CODE SNIPPET FROM DOTHUB_writePREPRO to define filename and logData
[pathstr, name, ~] = fileparts(nirsFileName);
ds = datestr(now,'yyyymmDDHHMMSS');
preproFileName = fullfile(pathstr,[name '_' ds '.prepro']);
logData(1,:) = {'Created on: '; ds};
logData(2,:) = {'Derived from data: ', nirsFileName};
logData(3,:) = {'Pre-processed using:', mfilename('fullpath')};

%DOTHUB_writePREPRO(preproFileName,logData,SD_3D,tDOD,dod,tHRF,dcAvg,dcStd)
[prepro, preproFileName] = DOTHUB_writePREPRO(preproFileName,logData,SD_3D,tHRF,dodAvg,tHRF,dcAvg,dcStd);

end
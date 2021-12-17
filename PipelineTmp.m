% DOTHUB Recon Pipeline Example 3.
%
% Individual optode positioning information, subject-specific .mshs file.
% LUMO DATA








%% function [prepro, preproFileName] = examplePreProcessingScript(nirsFileName)
% 
% load(nirsFileName,'-mat');
% 
% od = hmrIntensity2OD(d);
% od = hmrBandpassFilt(od,1/mean(diff(t)),0.01,0.5);
% dc = hmrOD2Conc(od,SD3D,[6 6]);
% [dcAvg, dcStd, tHRF] = hmrBlockAvg(dc,s,t,[-5 25]);
% dodAvg = DOTHUB_hmrHRFConc2OD(dcAvg, SD3D,[6 6]);
% 
% %USE CODE SNIPPET FROM DOTHUB_writePREPRO to define filename and logData
% [pathstr, name, ~] = fileparts(nirsFileName);
% ds = datestr(now,'yyyymmDDHHMMSS');
% preproFileName = fullfile(pathstr,[name '.prepro']);
% logData(1,:) = {'Created on: '; ds};
% logData(2,:) = {'Derived from data: ', nirsFileName};
% logData(3,:) = {'Pre-processed using:', mfilename('fullpath')};
% 
% %DOTHUB_writePREPRO(preproFileName,logData,SD3D,tDOD,dod,tHRF,dcAvg,dcStd)
% [prepro, preproFileName] = DOTHUB_writePREPRO(preproFileName,logData,SD3D,tHRF,dodAvg,tHRF,dcAvg,dcStd);
% 
% end

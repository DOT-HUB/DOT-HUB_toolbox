function [nirs, outname] = DOTHUB_cropNIRSfile(nirsFileName,timeStart,timeStop,outname)

%This function loads the specified nirs file and crops the contents between
%the specified timeStart and timeStop, then saves the result as a new .nirs
%file. Useful if a dataset contains more than one experiment that need
%separate processing pipelines.

% #########################################################################
% INPUTS ##################################################################
% nirsFileName      =   The .nirs filename to crop from

% timeStart         =   The time point (in seconds) at which to begin the new dataset.

% timeStop          =   The time point (in seconds) at which to end the new dataset.

% outputFilename    =   (Optional) specified output filename.

% OUTPUTS #################################################################

% nirs              = the nirs structure of the cropped dataset

% outname    = full filename of new saved .nirs file

% #########################################################################
% #########################################################################
% Version 0.1
% RJC, University College London, Nov 2021
% #########################################################################
% #########################################################################

% Parse inputs
if ~exist('nirsFileName','var')
    [filename, pathname, ~] = uigetfile('*.nirs','Select .nirs file');
    nirsFileName = [pathname '/' filename];
else
    [pathname,filename] = fileparts(nirsFileName);
    if isempty(pathname)
        pathname = pwd;
    end
    nirsFileName = fullfile(pathname,[filename '.nirs']);
end

if ~exist('outname','var')
    outname = [nirsFileName(1:end-5) '_crop_' num2str(timeStart) '-' num2str(timeStop) '.nirs'];
else
    [pathname,filename] = fileparts(outname);
    if isempty(pathname)
        pathname = pwd;
    end
    outname = fullfile(pathname,[filename '.nirs']);
end

%Load data
nirsOrig = load(nirsFileName,'-mat');

%Find start and stop indices;
[~,indStart] = min(abs(nirsOrig.t - timeStart));
[~,indStop] = min(abs(nirsOrig.t - timeStop));

%Crop d,t,s
nirs.d = nirsOrig.d(indStart:indStop,:);
nirs.t = nirsOrig.t(indStart:indStop);
nirs.t = nirs.t - nirs.t(1);
nirs.aux = zeros(length(nirs.t),8);

tmp = nirsOrig.s(indStart:indStop,:);
activeStim = find(any(tmp));
if ~isempty(activeStim)
    nirs.s = tmp(:,activeStim);
    nirs.CondNames = nirsOrig.CondNames{activeStim};
else
    nirs.s = zeros(size(nirs.t));
    nirs.CondNames = '';
end
nirs.SD = nirsOrig.SD;
if isfield(nirsOrig,'SD3D')
    nirs.SD3D = nirsOrig.SD3D;
end

fprintf(['Saving file to ' outname ' ...\n']);
save(outname,'-struct','nirs','-v7.3');






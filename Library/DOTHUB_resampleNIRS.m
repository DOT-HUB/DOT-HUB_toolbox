function [nirsResampled,outname] = DOTHUB_resampleNIRS(nirsFileName, fout)

% This function amends the sample rate of a .nirs tile to the frequency
% fout, then saves a new .nirs file at that new rate.
%
% INPUTS
% nirsFileName  =        The filepath of an existing .nirs file
% fout          =        The desired sample rate (defaults to 1Hz)
%
% OUTPUTS
% nirsResampled =        A nirs structure containing all components of the
%                        resampled .nirs file
% outname       =        The name to which the resampled .nirs file is
%                        saved.

%############################################################
if ~exist('nirsFileName','var')
    disp('Select .nirs file...');
    [file,path] = uigetfile('*.nirs','Select .nirs file');
    nirsFileName = fullfile(path,file);
end

if ~exist('fout','var')
    fout = 1;
    disp('Resampling to 1Hz by default');
end

%Load Data ############################################################
nirs = load(nirsFileName,'-mat');
fin = length(nirs.t)/max(nirs.t);

%Resample
ss = sum(nirs.s,1);
stim_tims = zeros(max(ss),size(nirs.s,2));
for i = 1:size(nirs.s,2)
    tmp = nirs.t(find(nirs.s(:,i)));
    stim_tims(1:length(tmp),i) = tmp ;
end

%Zero mean prior to resample
dtmp = nirs.d;
mnD = mean()

[d_res,t_res] = resample(nirs.d,nirs.t',fout);
%Remove zeros and negatives.
d_res(d_res<=0) = min(d_res(d_res>0));
[aux_res,t_res] = resample(nirs.aux,nirs.t',fout);
if size(t_res,1)~=1
    t_res = t_res';
end

s_res = zeros(length(t_res),size(nirs.s,2));
for i = 1:size(s_res,2)
    stim = stim_tims(stim_tims(:,i)>0,i);
    for j = 1:length(stim)
        tmp = stim(j);
        [~,k] = min(abs(t_res-tmp));
        s_res(k,i) = 1;
    end
end

nirsResampled = nirs;
nirsResampled.d = d_res;
nirsResampled.t = t_res;
nirsResampled.s = s_res;
nirsResampled.aux = aux_res;

[path,filename,~] = fileparts(nirsFileName);
outname = fullfile(path,[filename '_resampled_' num2str(fout) 'Hz.nirs']);
save(outname,'-struct','nirsResampled');




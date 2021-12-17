function DOTHUB_adultMRISegmentSPM(T1FileName,T2FileName,savePathName)

% This script uses SPM12 to segment the MRIs in 6 tissue types. Available
% sequences are T1-w, T2-w and FLAIR. Tissue types are fat, skull, air,
% CFS, WM and GM
%
%############################### INPUTS ###################################
% T1FileName    =       Full path of T1.nii file
% T2FileName    =       Full path of T2.nii file
% savePathName  =       (Optional) Full path of destination for /TPM tissue
%                       probability maps. Defaults to T1 directory
%
%############################# Dependencies ###############################
% SPM12 Matlab toolbox in path
%
% #########################################################################
% RJC, UCL and S Brigadoi, Uni Padova, May 2020
%
% ############################## Updates ##################################
% #########################################################################
%
% ############################### TO DO ###################################
% #########################################################################

T1FileName = '/Users/RCooper/Dropbox/Repositories/DOT-HUB_toolbox/ExampleData/Example1/T1.nii';
T2FileName = '/Users/RCooper/Dropbox/Repositories/DOT-HUB_toolbox/ExampleData/Example1/T2.nii';

% ########################## Input handling ###############################
% Ensure full paths
[T1Path,T1Name,ext] = fileparts(T1FileName);
if ~strcmpi(ext,'.nii')
    error('Expecting .nii files');
end
if isempty(T1Path)
    T1FileName = fullfile(pwd,T1Name);
    [T1Path,T1Name,ext] = fileparts(T1FileName);
end

if ~exist('T2FileName','var')
    T2Flag = 0;
elseif isempty(T2FileName)
    T2Flag = 0;
else
    T2Flag = 1;
    [T2Path,T2ame,ext] = fileparts(T2FileName);
    if ~strcmpi(ext,'.nii')
        error('Expecting .nii files');
    end
    if isempty(T2Path)
        T2FileName = fullfile(pwd,T2Name);
        [T2Path,T2Name,ext] = fileparts(T2FileName);
    end
end

if ~exist('savePathName','var')
    pathnameSubjSave = T1Path;
else
    pathnameSubjSave = savePathName;
end
if ~exist(pathnameSubjSave,'dir')
    mkdir(pathnameSubjSave);
end

% ########################## Set up SPM batch #############################
% Matlab job preparation
clear matlabbatch
spmTPMFileName = which('TPM.nii'); %SPM12 must already be on path for this to work
if isempty(spmTPMFileName)
    error('Please ensure SPM12 is in your matlab path.');
end

if ~T2Flag
    matlabbatch{1}.spm.spatial.preproc.channel(1).vols = {[T1FileName ',1']};
    matlabbatch{1}.spm.spatial.preproc.channel(1).biasreg = 0.1;%0.1
    matlabbatch{1}.spm.spatial.preproc.channel(1).biasfwhm = 20;%60
    matlabbatch{1}.spm.spatial.preproc.channel(1).write = [0 0];
    matlabbatch{1}.spm.spatial.preproc.channel(2).vols = {[T2FileName ',1']};
    matlabbatch{1}.spm.spatial.preproc.channel(2).biasreg = 0.1;
    matlabbatch{1}.spm.spatial.preproc.channel(2).biasfwhm = 20;%60
    matlabbatch{1}.spm.spatial.preproc.channel(2).write = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[spmTPMFileName ',1']};
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[spmTPMFileName ',2']};
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[spmTPMFileName ',3']};
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[spmTPMFileName ',4']};
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[spmTPMFileName ',5']};
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[spmTPMFileName ',6']};
    matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
    
    spm_jobman('initcfg');
    spm('defaults', 'FMRI');
    spm_jobman('run',matlabbatch);
else %T1 and T2
    % Matlab job preparation
    matlabbatch{1}.spm.spatial.preproc.channel(1).vols = {[T1FileName ',1']};
    matlabbatch{1}.spm.spatial.preproc.channel(1).biasreg = 0.1;%0.1
    matlabbatch{1}.spm.spatial.preproc.channel(1).biasfwhm = 60;%60
    matlabbatch{1}.spm.spatial.preproc.channel(1).write = [0 0];
    matlabbatch{1}.spm.spatial.preproc.channel(2).vols = {[T2FileName ',1']};
    matlabbatch{1}.spm.spatial.preproc.channel(2).biasreg = 0.1;
    matlabbatch{1}.spm.spatial.preproc.channel(2).biasfwhm = 20;%60
    matlabbatch{1}.spm.spatial.preproc.channel(2).write = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[spmTPMFileName ',1']};
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[spmTPMFileName ',2']};
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[spmTPMFileName ',3']};
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[spmTPMFileName ',4']};
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[spmTPMFileName ',5']};
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[spmTPMFileName ',6']};
    matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;%1
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;%1
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
    
    spm_jobman('initcfg');
    spm('defaults', 'FMRI');
    spm_jobman('run',matlabbatch);
end

% This function creates the mask. The original MRI (T1 and T2)
% was segmented with SPM12 which outputs 6 tissue segmentation mask: c1 for
% the GM, c2 for the WM, c3 for the CSF, c4 for the skull, c5 for the skin
% and c6 for air. 
%
% Input:  pathname         pathname to the directory of the subject, which should
%                          contain the tissue probability maps from SPM12.
%         subj             subject name
%
% Output: maskOriginal     segmented mask 
%         maskCut          segmented mask with 0 where neck and lower part
%                          of the head are present
%
% Written by S. Brigadoi 01/2017

% Load probability maps for the 6 tissue types and normalize to the maximum
% to get maps with values from 0 to 1.
pathnameSPM = fullfile(pathname,'SPM_segm','T1_T2_N4_ANTS','BiasReg01_FWHM20');
GM_prob  = load_untouch_nii(fullfile(pathnameSPM,'c10004-Sag_MPRAGE_1mm_ND_N4.nii'));
GM_prob = GM_prob.img/max(GM_prob.img(:)); 
WM_prob  = load_untouch_nii(fullfile(pathnameSPM,'c20004-Sag_MPRAGE_1mm_ND_N4.nii'));
WM_prob = WM_prob.img/max(WM_prob.img(:));
CSF_prob = load_untouch_nii(fullfile(pathnameSPM,'c30004-Sag_MPRAGE_1mm_ND_N4.nii'));
CSF_prob = CSF_prob.img/max(CSF_prob.img(:));
skull_prob = load_untouch_nii(fullfile(pathnameSPM,'c40004-Sag_MPRAGE_1mm_ND_N4.nii'));
skull_prob = skull_prob.img/max(skull_prob.img(:));
scalp_prob = load_untouch_nii(fullfile(pathnameSPM,'c50004-Sag_MPRAGE_1mm_ND_N4.nii'));
scalp_prob = scalp_prob.img/max(scalp_prob.img(:));
air_prob = load_untouch_nii(fullfile(pathnameSPM,'c60004-Sag_MPRAGE_1mm_ND_N4.nii'));
air_prob = air_prob.img/max(air_prob.img(:));

% Plot to check number assigned to tissue by SPM is correct
figure;
subplot(321)
imagesc(rot90(GM_prob(:,:,round(end/2)+40)))
axis square
title('GM')
subplot(322)
imagesc(rot90(WM_prob(:,:,round(end/2)+40)))
axis square
title('WM')
subplot(323)
imagesc(rot90(CSF_prob(:,:,round(end/2)+40)))
axis square
title('CSF')
subplot(324)
imagesc(rot90(skull_prob(:,:,round(end/2)+40)))
axis square
title('Skull')
subplot(325)
imagesc(rot90(scalp_prob(:,:,round(end/2)+40)))
axis square
title('Soft tissue')
subplot(326)
imagesc(rot90(air_prob(:,:,round(end/2)+40)))
axis square
title('Air')
set(gcf,'PaperPositionMode','auto','Position',[560 85 1085 863])

% Mask dimensions
dim = size(GM_prob);

% Compute segmentation mask and threshold it to get rid of external noise,
% thus creating a scalp mask. Air is not contained because the external
% space is segmented as air as well. 
mask_tot = WM_prob + CSF_prob + GM_prob + scalp_prob + skull_prob;
scalpMask = zeros(dim);
scalpMask(mask_tot > .99) = 1; 
scalpMask = imfill(scalpMask,'holes');

scalpMaskNew = zeros(dim);
for iSlice = 1:dim(3)
    CC = bwconncomp(scalpMask(:,:,iSlice));
    tmp = scalpMask(:,:,iSlice);
    if CC.NumObjects ~= 1
        nOfVoxels = zeros(CC.NumObjects,1);
        for iCC = 1:CC.NumObjects
            nOfVoxels(iCC,1) = length(CC.PixelIdxList{iCC});
        end
        [~,idx] = max(nOfVoxels);
        for iCC = 1:CC.NumObjects
            if iCC ~= idx
                tmp(CC.PixelIdxList{iCC}) = 0;
            end
        end
    end
    scalpMaskNew(:,:,iSlice) = reshape(tmp,dim(1),dim(2));
end

% Mask probability masks with scalp mask
GM_prob_brainMasked = zeros(dim);
GM_prob_brainMasked(scalpMaskNew==1) = GM_prob(scalpMaskNew == 1);

WM_prob_brainMasked = zeros(dim);
WM_prob_brainMasked(scalpMaskNew==1) = WM_prob(scalpMaskNew == 1);

CSF_prob_brainMasked = zeros(dim);
CSF_prob_brainMasked(scalpMaskNew==1) = CSF_prob(scalpMaskNew == 1);

scalp_prob_brainMasked = zeros(dim);
scalp_prob_brainMasked(scalpMaskNew==1) = scalp_prob(scalpMaskNew == 1);

skull_prob_brainMasked = zeros(dim);
skull_prob_brainMasked(scalpMaskNew==1) = skull_prob(scalpMaskNew == 1);

air_prob_brainMasked = zeros(dim);
air_prob_brainMasked(scalpMaskNew==1) = air_prob(scalpMaskNew == 1);

% Assign each voxel to the tissue with highest probability in that voxel
brain_tissues_t(:,:,:,1) = scalp_prob_brainMasked;
brain_tissues_t(:,:,:,2) = skull_prob_brainMasked;
brain_tissues_t(:,:,:,3) = CSF_prob_brainMasked;
brain_tissues_t(:,:,:,4) = GM_prob_brainMasked;
brain_tissues_t(:,:,:,5) = WM_prob_brainMasked;
brain_tissues_t(:,:,:,6) = air_prob_brainMasked;

[~,segm_prob] = max(brain_tissues_t,[],4);
segm_prob = segm_prob .*scalpMaskNew;
maskOriginal = segm_prob;

tiss_type = {'soft tissue','skull','CSF','GM','WM','air'};

% Plot to check segmentation mask
figure;
subplot(131)
imagesc(rot90(maskOriginal(:,:,round(end/2))));
axis equal tight
subplot(132)
imagesc(rot90(squeeze(maskOriginal(:,round(end/2),:))));
axis equal tight
subplot(133)
imagesc((rot90(squeeze(maskOriginal(round(end/2),:,:)))));
axis equal tight
set(gcf,'PaperPositionMode','auto','Position',[560 85 1085 863])

% Write the mask volume to the .mat file and the tissue type to
% the tissue type file
save(fullfile(pathname,['maskOriginal_' subj '.mat']),'maskOriginal');
save_tiss_type(fullfile(pathname,'mask_tiss_type.txt'), tiss_type);

% Remove voxels of neck and lower part of the head (set to 0)
maskCut = maskOriginal;
maskCut(:,:,1:50) = 0;

% Plot to check segmentation mask
figure;
subplot(131)
imagesc(rot90(maskCut(:,:,round(end/2))));
axis equal tight
subplot(132)
imagesc(rot90(squeeze(maskCut(:,round(end/2),:))));
axis equal tight
subplot(133)
imagesc((rot90(squeeze(maskCut(round(end/2),:,:)))));
axis equal tight
set(gcf,'PaperPositionMode','auto','Position',[560 85 1085 863])

% Write the mask volume to the .mat file 
save(fullfile(pathname,['maskCut_' subj '.mat']),'maskCut');

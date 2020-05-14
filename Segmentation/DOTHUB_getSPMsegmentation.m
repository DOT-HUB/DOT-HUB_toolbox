% This script uses SPM12 to segment the MRIs in 6 tissue types. Available
% sequences are T1-w, T2-w and FLAIR. Tissue types are fat, skull, air,
% CFS, WM and GM

addpath('C:\Users\brigadoi\Documents\MATLAB\spm12')
pathnameSPM = 'C:\Users\brigadoi\Documents\MATLAB\spm12\';

dataset = 'C:\Users\brigadoi\Dropbox\uNTS_MRI\AdultMRIs\';
subjs = dir(fullfile(dataset,'subj*'));

for s = 4%1:length(subjs)
    
    pathnameSubj = fullfile(dataset,subjs(s).name);
    pathnameSubjSave = fullfile(dataset,subjs(s).name,'SPM_segm');
    
    if ~exist(pathnameSubjSave,'dir')
        mkdir(pathnameSubjSave);
    end
    
    % Load MRIs
    T1 = fullfile(pathnameSubj,'0004-Sag_MPRAGE_1mm_ND_N4.nii');
    T2 = fullfile(pathnameSubj,'0002-t2-to-t1.nii.gz');
    %FLAIR = fullfile(pathnameSubj,'0006-t2_spc_FLAIR_fs_sag_p2_1600TI_DIS3D_NDToT1.nii.gz');
    gunzip(T2);
    %gunzip(FLAIR);
    
%     matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {'C:\Users\brigadoi\Dropbox\uNTS_MRI\Adult MRIs\subj05\0004-Sag_MPRAGE_1mm_ND.nii,1'};
%     matlabbatch{1}.spm.spatial.coreg.estwrite.source = {'C:\Users\brigadoi\Dropbox\uNTS_MRI\Adult MRIs\subj05\0002-t2_spc_ns_sag_p2_1mm_ND.nii,1'};
%     matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
%     matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
%     matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
%     matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
%     matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
%     matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
%     matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
%     matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
%     matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
% 
%     T2reg = fullfile(pathnameSubj,'r0002-t2_spc_ns_sag_p2_1mm_ND.nii');
    
    clear matlabbatch
    % Matlab job preparation
    matlabbatch{1}.spm.spatial.preproc.channel(1).vols = {[T1 ',1']};
    matlabbatch{1}.spm.spatial.preproc.channel(1).biasreg = 0.1;%0.1
    matlabbatch{1}.spm.spatial.preproc.channel(1).biasfwhm = 20;%60
    matlabbatch{1}.spm.spatial.preproc.channel(1).write = [0 0];
    matlabbatch{1}.spm.spatial.preproc.channel(2).vols = {[T2(1:end-3) ',1']};
    matlabbatch{1}.spm.spatial.preproc.channel(2).biasreg = 0.1;
    matlabbatch{1}.spm.spatial.preproc.channel(2).biasfwhm = 20;%60
    matlabbatch{1}.spm.spatial.preproc.channel(2).write = [0 0];
    %matlabbatch{1}.spm.spatial.preproc.channel(3).vols = {[FLAIR(1:end-3) ',1']};
    %matlabbatch{1}.spm.spatial.preproc.channel(3).biasreg = 0.001;
    %matlabbatch{1}.spm.spatial.preproc.channel(3).biasfwhm = 60;
    %matlabbatch{1}.spm.spatial.preproc.channel(3).write = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[pathnameSPM '\tpm\TPM.nii,1']};
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[pathnameSPM '\tpm\TPM.nii,2']};
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[pathnameSPM '\tpm\TPM.nii,3']};
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[pathnameSPM '\tpm\TPM.nii,4']};
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[pathnameSPM '\tpm\TPM.nii,5']};
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[pathnameSPM '\tpm\TPM.nii,6']};
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
    
    clear matlabbatch
    % Matlab job preparation
    matlabbatch{1}.spm.spatial.preproc.channel(1).vols = {[T1 ',1']};
    matlabbatch{1}.spm.spatial.preproc.channel(1).biasreg = 0.1;%0.1
    matlabbatch{1}.spm.spatial.preproc.channel(1).biasfwhm = 60;%60
    matlabbatch{1}.spm.spatial.preproc.channel(1).write = [0 0];
    %matlabbatch{1}.spm.spatial.preproc.channel(2).vols = {[T2(1:end-3) ',1']};
    %matlabbatch{1}.spm.spatial.preproc.channel(2).biasreg = 0.1;
    %matlabbatch{1}.spm.spatial.preproc.channel(2).biasfwhm = 20;%60
    %matlabbatch{1}.spm.spatial.preproc.channel(2).write = [0 0];
    %matlabbatch{1}.spm.spatial.preproc.channel(3).vols = {[FLAIR(1:end-3) ',1']};
    %matlabbatch{1}.spm.spatial.preproc.channel(3).biasreg = 0.001;
    %matlabbatch{1}.spm.spatial.preproc.channel(3).biasfwhm = 60;
    %matlabbatch{1}.spm.spatial.preproc.channel(3).write = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[pathnameSPM '\tpm\TPM.nii,1']};
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[pathnameSPM '\tpm\TPM.nii,2']};
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[pathnameSPM '\tpm\TPM.nii,3']};
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[pathnameSPM '\tpm\TPM.nii,4']};
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[pathnameSPM '\tpm\TPM.nii,5']};
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[pathnameSPM '\tpm\TPM.nii,6']};
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


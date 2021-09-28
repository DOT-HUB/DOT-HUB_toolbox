function [GMSensitivity, GMmask] = DOTHUB_MakeGMSensitivityMap(jac,SD,Thresh)

%EM Frijia's function to create a normalized GM sensitivity map, and a
%thresholded mask.

% #########################################################################
% INPUTS ##################################################################
% jac  =            The .jac variable includes the GM and volume jacobians

% J_GM =            The GM jacobian values taken from DOTHUB_Make_ToastJacobian, size
%                   nchan x nnodes GM

% J_vol =           The volume jacobian values taken from DOTHUB_Make_ToastJacobian, size
%                   nchan x nnodes volume

% SD =              The SD file for the array used in the study.  Note that
%                   the srcpos and detpos form the SD file are ignored in
%                   favour of the above defined inputs

% Thresh =          The threshold at which the sensitivity should be cut to
%                   create the binary mask


% OUTPUTS #################################################################

% GMSensitivity =   The normalized sum of the jacobian values of good channels on the GM

% GMMask =          The binary mask created by thresholding GMSensitivity
%                   at the normalized value Thresh.

% #########################################################################
% #########################################################################
% Version 0.2
% EMF, University College London, December 2020 (ho ho hbo)
% #########################################################################
% #########################################################################

%Extract the J_GM and J_vol from the .jac file for the 2 wavelengths
%GM values
J_GM_wav1=J{1}.gm;
J_GM_wav2=J{2}.gm;
%Volume values
J_vol_wav1=J{1}.vol;
J_vol_wav2=J{2}.vol;

%First, find the maximum value of the sum of J_vol for all good channels
%(maxJvol). (BE CAREFUL AS THESE ARE ALL NEGATIVE!)

J_vol_wav1_cropped = J_vol_wav1(SD.MeasListAct(1:end/2)==1,:);
J_vol_wav2_cropped = J_vol_wav2(SD.MeasListAct(end/2+1:end)==1,:);

MaxJvol_wav1_cropped = max(abs(J_vol_wav1_cropped),[],2);
MaxJvol_wav2_cropped = max(abs(J_vol_wav2_cropped),[],2);

%Then divide J_GM by maxJvol on a channel wise basis. This 
%produce GMSensitivity

J_GM_wav1_cropped = abs(J_GM_wav1(SD.MeasListAct(1:end/2)==1,:));
J_GM_wav2_cropped = abs(J_GM_wav2(SD.MeasListAct(end/2+1:end)==1,:));

nnodes = size(J_GM_wav1_cropped,2);
J_GM_wav1_norm = J_GM_wav1_cropped./repmat(MaxJvol_wav1_cropped,1,nnodes);
J_GM_wav2_norm = J_GM_wav2_cropped./repmat(MaxJvol_wav2_cropped,1,nnodes);

%as we are now normalizing within channel, the GMSensitivity is now of a
%confusing scale.
GMSensitivity_wav1 = sum(J_GM_wav1_norm);
GMSensitivity_wav2 = sum(J_GM_wav2_norm);
GMSensitivity=(GMSensitivity_wav1+GMSensitivity_wav2)./2;

GMmask1 = any(J_GM_wav1_norm>Thresh,1);%If any channel is sensitive at both wavs, we are sensitive
GMmask2 = any(J_GM_wav2_norm>Thresh,1);
GMmask= GMmask1 & GMmask2;


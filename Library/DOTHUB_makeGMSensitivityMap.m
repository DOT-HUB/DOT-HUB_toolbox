function [GMSensitivity, GMmask] = DOTHUB_makeGMSensitivityMap(jac,rmap,SD)

% Function to create a normalized GM sensitivity map, and a
% thresholded mask based on the 'covered' definition from Array Designer.
% See Brigadoi et al., Array Designer: automated optimized array design for
% functional near-infrared spectroscopy. Neurophotonics (2018).
%
% #########################################################################
% INPUTS ##################################################################
% jac  =            The jac structure or pathname. .jac includes the GM and
%                   volume jacobians
%
% rmap =            The rmap associated with the .jac - either as a
%                   structure or pathname
%
% SD   =            The SD file for the array used in the study containing
%                   the appropriate MeasListAct.
%
% OUTPUTS #################################################################
%
% GMSensitivity =   log10(sum(abs(J))./max(sum(abs(J))): The normalized total 
%                   sensitivity on a log scale
%
% GMMask =          The binary mask created by thresholding GMSensitivity
%                   on the basis defined by Brigadoi et al.
%
% #########################################################################
% #########################################################################
% Version 0.2
% RJC and EMF, University College London, December 2020 (ho ho hbo), Jan 2021 
%
% TO DO:
% Expand to N wavelengths
% #########################################################################
% #########################################################################
%
% MANAGE VARIABLES
% #########################################################################
if ischar(jac)
    jacFileName = jac;
    [jacPath, jacName,~] = fileparts(jacFileName);
    if isempty(jacPath) %Ensure full path.
        jacPath = pwd;
        jacFileName = fullfile(jacPath,[jacName '.jac']);
    end
    jac = load(jacFileName,'-mat');
end

if ischar(rmap)
    rmapFileName = rmap;
    [rmapPath, rmapName,~] = fileparts(rmapFileName);
    if isempty(rmapPath) %Ensure full path.
        rmapPath = pwd;
        rmapFileName = fullfile(rmapPath,[rmapName '.rmap']);
    end
    rmap = load(rmapFileName,'-mat'); %head, gm, scalp, SD3Dmesh, vol2gm, logData
end

if ischar(SD)
    SDFileName = SD;
    [SDPath, SDName, ext] = fileparts(SDFileName);
    if isempty(SDPath) %Ensure full path.
        SDPath = pwd;
        SDFileName = fullfile(SDPath,[SDName ext]);
    end
    tmp = load(SDFileName,'-mat'); %Could be SD or SD3D, so need to account for both
    if isfield(tmp,'SD')
        SD = tmp.SD;
    elseif isfield(tmp,'SD3D')
        SD = tmp.SD3D;
    else
        error('Cannot parse SD variable - neither SD or SD3D structure found')
    end
end

%% 

%Determine J_vol
if isempty(jac.J{1}.vol) %If using a basis
    %Need to replicate toast mesh to allow transform from basis to mesh
    fprintf('Rebuilding TOAST mesh...\n');
    eltp = ones(length(rmap.headVolumeMesh.elem),1)*3;
    hMesh = toastMesh(rmap.headVolumeMesh.node(:,1:3),rmap.headVolumeMesh.elem(:,1:4),eltp);
    hBasis = toastBasis(hMesh,jac.basis,jac.basis*2);
    BB = hMesh.BoundingBox;
    vVol = prod((BB(2,:)-BB(1,:))./jac.basis);
    %It's too slow to convert to a volume jac, so calculate in basis
    %Calculate sum over volume J val by channel
    volSum_wav1 = sum(abs(jac.J{1}.basis),2);
    volSum_wav2 = sum(abs(jac.J{2}.basis),2);
    %for i = 1:size(jac.J{1}.basis,1)
    %    jac.J{1}.vol(i,:) = hBasis.Map('S->M',jac.J{1}.basis(i,:));
    %    jac.J{2}.vol(i,:) = hBasis.Map('S->M',jac.J{2}.basis(i,:));
    %end
else
    volSum_wav1 = sum(abs(jac.J{1}.vol),2);
    volSum_wav2 = sum(abs(jac.J{2}.vol),2);
    vVol = mean(elemvolume(rmap.headVolumeMesh.node(:,1:3),rmap.headVolumeMesh.elem(:,1:4)));
end

%Extract GM values and normalise
J_GM_wav1 = abs(jac.J{1}.gm)./volSum_wav1;
J_GM_wav2 = abs(jac.J{2}.gm)./volSum_wav2;

%% Calculate threshold we define as 'sensitive'
pthresh = 1; % % change in measured signal to be considered 'measured'
actVol = 1000; %mm3 volume of a theoretical activation (1 cm3)
deltaMua = 0.001; %mm-1 change in mua of theoretical activation (10%)
tmp = vVol./actVol;
SenseThresh = log((100+pthresh)/100)*tmp;

J_GM_wav1_cropped = J_GM_wav1(SD.MeasListAct(1:end/2)==1,:);
J_GM_wav2_cropped = J_GM_wav2(SD.MeasListAct(end/2+1:end)==1,:);

GMSensitivity_wav1 = sum(J_GM_wav1_cropped);
GMSensitivity_wav2 = sum(J_GM_wav2_cropped);

%Create mask
GMmask1 = any(J_GM_wav1_cropped > SenseThresh);
GMmask2 = any(J_GM_wav2_cropped > SenseThresh);
GMmask = GMmask1 & GMmask2;

%Take mean sensitivity across two wavelengths and threshold by mask
GMSensitivity = (GMSensitivity_wav1+GMSensitivity_wav2)./2;
GMSensitivity = log10(GMSensitivity./max(GMSensitivity));

figure
subplot(1,2,1)
set(gcf,'color','w');
DOTHUB_plotSurfaceImage(rmap.gmSurfaceMesh,double(GMmask))
cmap = [0.7 0.7 0.7;0.6 0 0];
colormap(cmap);caxis([0 1]);colorbar off;
title('Array Sensitivity Mask');

subplot(1,2,2)
set(gcf,'color','w');
DOTHUB_plotSurfaceImage(rmap.gmSurfaceMesh,GMSensitivity)
load('greyJet.mat');
cmap = greyJet(end/2:end,:);
colormap(cmap);
cb = colorbar;
caxis([-2 0]);
ylabel(cb,'log10(Normalised Sensitivity)');
title('Array Sensitivity');


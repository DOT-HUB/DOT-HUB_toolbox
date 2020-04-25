function [dot, dotFileName] = DOTHUB_reconstruction(prepro,jac,invjac,rmap,varargin)

% DOTHUB_reconstruction.m
% Linear reconstruction of concentration changes from prepro data.

% ####################### INPUTS ##########################################
% prepro    =  prepro structure or path to .prepro

% jac       =  jac structure or path to .jac. If parsed empty (i.e. as []) and
%              invjac is parsed, invjac is used. One of the two must be parsed.
%              If jac is parsed and jac.basis is not empty, a toast mesh basis is
%              assumed and rebuilt in order to create the volume and then GM
%              images. If you don't have jac, please use DOTHUB_makeToastJacobian

% invjac    =  invjac structure or path to .invjac. If parsed empty (i.e. as [])
%              invjac is calculated but not saved. If you wish to
%              pre-calculate invjac, please use DOTHUB_invertJacobian.m.
%              built (and both saved?)

% rmap      =  rmap structure or path to .rmap

% varargin  =  optional input pairs:
%              'reconMethod' - 'multispec' or 'standard' (default 'standard');
%                   Specifying whether to construct and invert a multispectral
%                   jacobian or whether to recontruct each wavelength
%                   separately and then combine them
%              'regMethod' - 'tikhonov' or 'covariance' or 'spatial' (default 'tikhonov')
%                   Regularization method. See DOTHUB_invertJacobian for
%                   more details
%              'hyperParameter' - numerical value or vector (for 'spatial') (default 0.01);
%                   Regularization hyperparamter. See DOTHUB_invertJacobian for more details
%              'imageType' - 'haem', 'mua' or 'both' (default 'haem');
%                   Determines whether to output haemoglobin images, mua images or
%                   both. Calls 'mua' and 'both' must be coupled with reconMethod 'standard';
%              'saveVolumeImages' - 'true' or 'false' (default 'true');
%                   Flag whether to output volume images to dot structure in addition GM.
%              'saveFlag' - 'true' or 'false' (default 'true');
%                   Flag whether to save the output images to a .dot file
%                   (default true)

% ######################### OUTPUTS #######################################
%[dot, dotFileName]
% ####################### Dependencies ####################################

% #########################################################################
% RJC, UCL, April 2020

% MANAGE VARIABLES
% #########################################################################
varInputs = inputParser;
varInputs.CaseSensitive = false;
addParameter(varInputs,'reconMethod','standard',@ischar);
addParameter(varInputs,'regMethod','tikhonov',@ischar);
addParameter(varInputs,'hyperParameter',0.01,@isnumeric);
validateImageType = @(x) assert(any(strcmpi({'haem','mua','both'},x)));
addParameter(varInputs,'imageType','haem',validateImageType);
validateFlag = @(x) assert(x==0 || x==1);
addParameter(varInputs,'saveVolumeImages',false,validateFlag);
addParameter(varInputs,'saveFlag',true,validateFlag);
parse(varInputs,varargin{:});

%Basic error handling
varInputs = varInputs.Results;
if strcmpi(varInputs.regMethod,'spatial')
    if length(varInputs.hyperParameter)<2
        error('For spatial regularization, the hyperParameter must be a vector');
    end
end

if (strcmpi(varInputs.imageType,'mua') || strcmpi(varInputs.imageType,'both')) && ~strcmpi(varInputs.reconMethod,'standard')
    error('To call for images of mua requires reconMethod = standard');
end

%Print selected parameters
fnames = fieldnames(varInputs);
fprintf(['Input parameters...\n'])
for i = 1:numel(fnames)
    fprintf([fnames{i} ' = ' num2str(getfield(varInputs,fnames{i})) '\n'])
end

%Load core variables if parsed as paths
if ischar(prepro)
    load(prepro,'-mat');
end

if ischar(jac)
    load(jac,'-mat');
end

if ischar(invjac)
    load(invjac,'-mat');
end

if ischar(rmap)
    load(rmap,'-mat');
end

% #########################################################################
% Calculate Inverted Jacobian if not parsed
if isempty(invjac)
    varargininvjac = {'reconMethod',varInputs.reconMethod,'regMethod',varInputs.regMethod,...
        'hyperParameter',varInputs.hyperParameter,'rmap',rmap,'saveFlag',false};
    [invjac, ~] = DOTHUB_invertJacobian(jac,prepro,varargininvjac);
end

% #########################################################################
% Set data in TOAST format
% Convert data into toast style (toast wants = ln(Intensity_active)-ln(intensity_baseline)
% Parsed data is OD (i.e. data_OD = -ln(intensity_active/mean));
% Also, crop out bad channels;
datarecon = -prepro.dod(:,prepro.SD_3D.MeasListAct==1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create matrices to save recontruction results
nNodeNat = size(rmap.headVolumeMesh.node,1); %The spatial size of the native image space (basis or full volume)
nNodeVol = size(invjac.invJ{1},1)/2; %The node count of the volume mesh
nNodeGM = size(rmap.gmSurfaceMesh.node,1); %The node count of the GM mesh
nFrames = size(datarecon,1);
SD_3D = prepro.SD_3D;
wavelengths = SD_3D.Lambda;
nWavs = length(prepro.SD_3D.Lambda);

% pre-assign large things
if ~(strcmpi(varInputs.imageType,'mua')) %only need haem variables if imageType not 'mua'
    hbo.vol = zeros(nFrames,nNodeVol);
    hbr.vol = zeros(nFrames,nNodeVol);
    hbo.gm = zeros(nFrames,nNodeGM);
    hbr.gm = zeros(nFrames,nNodeGM);
end
if (strcmpi(varInputs.imageType,'mua') || strcmpi(varInputs.imageType,'both')) %If mua images are called for, pre-assign
    muaFlag = 1;
    mua = cell(nWavs,1);
    for wav = 1:nWavs
        mua{wav}.vol = zeros(nFrames,nNodeVol);
        mua{wav}.gm = zeros(nFrames,nNodeGM);
    end
else
    muaFlag = 0;
end

if ~isempty(jac.basis) %If using a basis
    basisFlag = 1;
    %Need to replicate toast mesh to allow transform from basis to mesh
    fprintf('Rebuilding TOAST mesh...\n');
    eltp = ones(length(rmap.headVolumeMesh.elem),1)*3;
    hMesh = toastMesh(rmap.headVolumeMesh.node(:,1:3),rmap.headVolumeMesh.elem(:,1:4),eltp);
    hBasis = toastBasis(hMesh,basis,basis*2);
else
    basisFlag = 0;
end

%###################### reconMethod = multispectral #######################
%##########################################################################
if strcmpi(varInputs.ReconMethod,'multispectral')
    fprintf('Reconstructing images...\n');
    for frame = 1:nFrames
        fprintf('Reconstructing frame %d of %d\n',frame,nFrames);
        
        data = datarecon(frame,:);
        img = invJ * data';
        
        if basisFlag
            hbo_tmp = img(1:nNode);
            hbr_tmp = img(nNode+1:2*nNode);
            hbo.vol(frame,:) = hBasis.Map('S->M',hbo_tmp');
            hbr.vol(frame,:) = hBasis.Map('S->M',hbr_tmp');
        else
            hbo.vol(frame,:) = img(1:nNode);
            hbr.vol(frame,:) = img(nNode+1:2*nNode);
        end
        hbo.gm(frame,:) = (vol2gm*hbo.vol(frame,:)');
        hbr.gm(frame,:) = (vol2gm*hbr.vol(frame,:)');
    end
end
end

%###################### reconMethod = standard ############################
%##########################################################################
if strcmpi(varInputs.ReconMethod,'standard')
    fprintf('Reconstructing images...\n');
    
    if ~strcmpi(varInputs.imageType,'mua') %Need to calculate haem images except if mua called
        Eall = [];
        for i = 1:nWavs
            Etmp = GetExtinctions(wavelengths(i));
            Etmp = Etmp(1:2); %HbO and HbR only
            Eall = [Eall; Etmp./1e7]; %This will be nWavs x 2;
        end
        Eallinv = pinv(Eall); %This will be (n chromophores(2)) x nWavs;
    end
    
    for frame = 1:nFrames
        fprintf('Reconstructing frame %d of %d\n',frame,nFrames);
        
        muaImageAll = zeros(nWav,nNodeNat);
        for wav = 1:nWavs
            data = datarecon(frame,SD_3D.MeasList(:,4)==wav);
            invJtmp = invJ{wav};
            tmp = invJtmp * data(frame,:)';
            muaImageAll(wav,:) = tmp; %This will be nWavs * nNodeNat
        end
        
        if ~strcmpi(varInputs.imageType,'mua') %Need to calculate haem images unless imageType = 'mua'
            
            %##### CHECK THIS ########
            img = Eallinv*muaImage;% Should be (chromophores by nWavs)*(nWavs by nNodeNat) = chromophore x node
            %#########################
            
            if basisFlag %Map first from basis to volume
                hbo_tmp = img(1,:);
                hbr_tmp = img(2,:);
                hbo.vol(frame,:) = hBasis.Map('S->M',hbo_tmp);
                hbr.vol(frame,:) = hBasis.Map('S->M',hbr_tmp);
            else %Already in volume
                hbo.vol(frame,:) = img(1,:);
                hbr.vol(frame,:) = img(2,:);
            end
            hbo.gm(frame,:) = (vol2gm*hbo.vol(:))';
            hbr.gm(frame,:) = (vol2gm*hbr.vol(:))';
        end
        
        if muaFlag %Calculate mua images if imageType = 'both' or 'mua'
            for wav = 1:nWavs
                if basisFlag
                    tmp = hBasis.Map('S->M',muaImageAll(wav,:));
                    mua{wav}.vol(frame,:) = tmp;
                    mua{wav}.gm(frame,:) = (vol2gm*tmp)';
                else
                    tmp = muaImageAll(wav,:);
                    mua{wav}.vol(frame,:) = tmp;
                    mua{wav}.gm(frame,:) = (vol2gm*tmp)';
                end
            end
        end
    end
end

if ~varInputs.saveVolumeImages
    hbo.vol = [];
    hbr.vol = [];
    for wav = 1:nWavs;
        mua{wav}.vol = [];
    end
end

%################ Create dot structure and write .dot #####################
%##########################################################################
[pathstr, name, ~] = fileparts(prepro.fileName);
ds = datestr(now,'yyyymmDDHHMMSS');
dotFileName = fullfile(pathstr,[name '_' ds '.dot']);
logData(1,:) = {'Created on: ', ds};
logData(2,:) = {'Associated prepro file: ', prepro.fileName};
logData(3,:) = {'Associated invjac file: ', invjac.fileName};
logData(4,:) = {'reconMethod: ', varInputs.reconMethod};
logData(5,:) = {'regMethod: ', varInputs.regMethod};
logData(6,:) = {'hyperParameter: ', varInputs.regMethod};

[dot, dotFileName] = DOTHUB_writeDOT(dotFileName,logData,hbo,hbr,mua,timebase,varInputs.saveFlag);



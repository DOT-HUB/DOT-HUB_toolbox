function [dotimg, dotimgFileName] = DOTHUB_reconstruction(prepro,jac,invjac,rmap,varargin)

% DOTHUB_reconstruction.m
% Linear reconstruction of concentration changes from prepro data.
%
% ####################### INPUTS ##########################################
% prepro    =  prepro structure or path to .prepro
%
% jac       =  jac structure or path to .jac. If parsed empty (i.e. as []) and
%              invjac is parsed, invjac is used. One of the two must be parsed.
%              If jac is parsed and jac.basis is not empty, a toast mesh basis is
%              assumed and rebuilt in order to create the volume and then GM
%              images. If you don't have jac, please use DOTHUB_makeToastJacobian
%
% invjac    =  invjac structure or path to .invjac. If parsed empty (i.e. as [])
%              invjac is calculated but not saved. If you wish to
%              pre-calculate invjac, please use DOTHUB_invertJacobian.m.
%              built (and both saved?). Note that if invJac is parsed, any
%              vaargin inputs related to inversion are superceded by those
%              saved within the invjac
%
% rmap      =  rmap structure or path to .rmap
%
% varargin  =  optional input pairs:
%              'reconMethod' - 'multispec' or 'standard' (default 'standard');
%                   Specifying whether to construct and invert a multispectral
%                   jacobian or whether to recontruct each wavelength
%                   separately and then combine them
%              'reconSpace' - 'volume' or 'cortex' (default 'volume');
%              'regMethod' - 'tikhonov' or 'covariance' or 'spatial' (default 'tikhonov')
%                   Regularization method. See DOTHUB_invertJacobian for
%                   more details
%              'hyperParameter' - numerical value or vector (for 'spatial') (default 0.01);
%                   Regularization hyperparamter. See DOTHUB_invertJacobian for more details
%              'imageType' - 'haem', 'mua' or 'both' (default 'haem');
%                   Determines whether to output haemoglobin images, mua images or
%                   both. Calls 'mua' and 'both' must be coupled with reconMethod 'standard';
%              'saveVolumeImages' - 'true' or 'false' (default 'true');
%                   Flag whether to output volume images to dotimg structure in addition GM.
%              'saveFlag' - 'true' or 'false' (default 'true');
%                   Flag whether to save the output images to a .dotimg file
%                   (default true)
%
% ######################### OUTPUTS #######################################
% [dotimg, dotimgFileName]
% ####################### Dependencies ####################################
%
% #########################################################################
% RJC, UCL, April 2020

fprintf('################### Running DOTHUB_reconstruction ###################\n');

% MANAGE VARIABLES
% #########################################################################
varInputs = inputParser;
varInputs.CaseSensitive = false;
validateReconMethod = @(x) assert(any(strcmpi({'standard','multispectral'},x)));
addParameter(varInputs,'reconMethod','standard',validateReconMethod);
validateSpace = @(x) assert(any(strcmpi({'volume','cortex'},x)));
addParameter(varInputs,'reconSpace','volume',validateSpace);
validateRegMethod = @(x) assert(any(strcmpi({'tikhonov','covariance','spatial'},x)));
addParameter(varInputs,'regMethod','tikhonov',validateRegMethod);
addParameter(varInputs,'hyperParameter',0.01,@isnumeric);
validateImageType = @(x) assert(any(strcmpi({'haem','mua','both'},x)));
addParameter(varInputs,'imageType','haem',validateImageType);
validateFlag = @(x) assert(x==0 || x==1);
addParameter(varInputs,'saveVolumeImages',true,validateFlag);
addParameter(varInputs,'saveFlag',true,validateFlag);
parse(varInputs,varargin{:});
varInputs = varInputs.Results;

%Basic error handling
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
    preproFileName = prepro;
    prepro = load(preproFileName,'-mat');
    if ndims(prepro.dod)==3
        nCond = size(prepro.dod,3);
    else
        nCond = 1;
        tmp = prepro.dod;
        prepro = rmfield(prepro,'dod');
        prepro.dod(:,:,1) = tmp; %Easier to assign a dummy condition dimension
    end
else %loaded as structure
    if ndims(prepro.dod)==3
        nCond = size(prepro.dod,3);
    else
        nCond = 1;
        tmp = prepro.dod;
        prepro = rmfield(prepro,'dod');
        prepro.dod(:,:,1) = tmp; %Easier to assign a dummy condition dimension
    end
end


if ischar(jac)
    jacFileName = jac;
    jac = load(jacFileName,'-mat');
end

if ischar(invjac)
    invJacFileName = invjac;
    invjac = load(invJacFileName,'-mat');
end

if ischar(rmap)
    rmapFileName = rmap;
    rmap = load(rmapFileName,'-mat');
end

if (strcmpi(varInputs.imageType,'mua') || strcmpi(varInputs.imageType,'both')) %mua images called for
    muaFlag = 1;
else
    muaFlag = 0;
end

% #########################################################################
% Calculate Inverted Jacobian if not parsed
if isempty(invjac)
    varargininvjac = {'reconMethod',varInputs.reconMethod,'reconSpace',varInputs.reconSpace,'regMethod',varInputs.regMethod,...
        'hyperParameter',varInputs.hyperParameter,'rmap',rmap,'saveFlag',false};
    [invjac, ~] = DOTHUB_invertJacobian(jac,prepro,varargininvjac{:},'saveFlag',false);
else
    %invjac is being loaded directly
    %Overwrite varInputs to match those of specified invjac;
    fprintf('invjac parsed directly, reverted to invjac input parameters...\n');
    fnames = fieldnames(varInputs);
    toOverWrite = {'hyperParameter','reconMethod','regMethod','reconSpace'};
    for i = 1:size(invjac.logData,1)
        lgInd = find(strcmpi(toOverWrite,invjac.logData{i,1}));
        if lgInd
            varInputs = setfield(varInputs,toOverWrite{lgInd},invjac.logData{i,2});
        end
    end
    fprintf(['***INPUT PARAMETERS***\n'])
    fnames = fieldnames(varInputs);
    for i = 1:numel(fnames)
        if strcmpi(fnames{i},'rmap');continue;end
        fprintf([fnames{i} ' = ' num2str(getfield(varInputs,fnames{i})) '\n'])
    end
    fprintf('\n');
end

% #########################################################################
% Set data in TOAST format
% Convert data into toast style (toast wants = ln(Intensity_active)-ln(intensity_baseline)
% Parsed data is OD (i.e. data_OD = -ln(intensity_active/mean));
datarecon = -prepro.dod;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unpack variables, define useful counts
nNodeVol = size(rmap.headVolumeMesh.node,1);  %The node count of the volume mesh
nNodeNat = size(invjac.invJ{1},1);%The spatial size of the native image space (basis or full volume)
if strcmpi(varInputs.reconMethod,'multispectral') %In multispectral case, invJ is double the length.
    nNodeNat = size(invjac.invJ{1},1)/2;
end
nNodeGM = size(rmap.gmSurfaceMesh.node,1); %The node count of the GM mesh
nFrames = size(datarecon,1);
vol2gm = rmap.vol2gm;
SD3D = prepro.SD3D;
wavelengths = SD3D.Lambda;
nWavs = length(prepro.SD3D.Lambda);

% pre-assign large things
hbo.vol = zeros(nFrames,nNodeVol,nCond);
hbr.vol = zeros(nFrames,nNodeVol,nCond);
hbo.gm = zeros(nFrames,nNodeGM,nCond);
hbr.gm = zeros(nFrames,nNodeGM,nCond);
mua = cell(nWavs,1);
for wav = 1:nWavs
    mua{wav}.vol = zeros(nFrames,nNodeVol,nCond);
    mua{wav}.gm = zeros(nFrames,nNodeGM,nCond);
end

if ~isempty(invjac.basis) %If using a basis
    basisFlag = 1;
    %Need to replicate toast mesh to allow transform from basis to mesh
    fprintf('Rebuilding TOAST mesh...\n');
    eltp = ones(length(rmap.headVolumeMesh.elem),1)*3;
    hMesh = toastMesh(rmap.headVolumeMesh.node(:,1:3),rmap.headVolumeMesh.elem(:,1:4),eltp);
    hBasis = toastBasis(hMesh,invjac.basis,invjac.basis*2);
else
    basisFlag = 0;
end

%###################### reconMethod = multispectral #######################
%##########################################################################
if strcmpi(varInputs.reconMethod,'multispectral')
    fprintf('Reconstructing images...\n');
    for cond = 1:nCond
        fprintf(['Reconstructing condition ' num2str(cond) ' of ' num2str(nCond) '...\n']);
        for frame = 1:nFrames
            fprintf('Reconstructing frame %d of %d\n',frame,nFrames);
            
            dataTmp = squeeze(datarecon(frame,SD3D.MeasListAct==1,cond));
            img = invjac.invJ{1} * dataTmp'; %invjac.invJ should only have one entry.
            
            if basisFlag %basis to volume to gm
                hbo_tmp = img(1:end/2);
                hbr_tmp = img(end/2+1:end);
                hbo_tmpVol = hBasis.Map('S->M',hbo_tmp);
                hbr_tmpVol = hBasis.Map('S->M',hbr_tmp);
                hbo_tmpGM = (vol2gm*hbo_tmpVol(:));
                hbr_tmpGM = (vol2gm*hbr_tmpVol(:));
                
                hbo.vol(frame,:,cond) = hbo_tmpVol;
                hbr.vol(frame,:,cond) = hbr_tmpVol;
                hbo.gm(frame,:,cond) = hbo_tmpGM;
                hbr.gm(frame,:,cond) = hbr_tmpGM;
            else
                if strcmpi(varInputs.reconSpace,'cortex') %GM to GM only
                    hbo.gm(frame,:,cond) = img(1:end/2);
                    hbr.gm(frame,:,cond) = img(end/2+1:end);
                else                                      %Vol to GM
                    hbo_tmpVol = img(1:end/2);
                    hbr_tmpVol = img(end/2+1:end);
                    hbo_tmpGM = (vol2gm*hbo_tmpVol(:));
                    hbr_tmpGM = (vol2gm*hbr_tmpVol(:));
                    
                    hbo.vol(frame,:,cond) = hbo_tmpVol;
                    hbr.vol(frame,:,cond) = hbr_tmpVol;
                    hbo.gm(frame,:,cond) = hbo_tmpGM;
                    hbr.gm(frame,:,cond) = hbr_tmpGM;
                end
            end
        end
    end
end

%###################### reconMethod = standard ############################
%##########################################################################
if strcmpi(varInputs.reconMethod,'standard')
    fprintf('Reconstructing images...\n');
    for cond = 1:nCond
        fprintf(['Reconstructing condition ' num2str(cond) ' of ' num2str(nCond) '...\n']);
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
            
            muaImageAll = zeros(nWavs,nNodeNat);
            for wav = 1:nWavs
                dataTmp = squeeze(datarecon(frame,SD3D.MeasList(:,4)==wav & SD3D.MeasListAct==1,cond));
                invJtmp = invjac.invJ{wav};
                tmp = invJtmp * dataTmp';
                muaImageAll(wav,:) = tmp; %This will be nWavs * nNodeNat
            end
            
            if ~strcmpi(varInputs.imageType,'mua') %Need to calculate haem images unless imageType = 'mua'
                
                %##### CHECK THIS ########
                img = Eallinv*muaImageAll;% Should be (chromophores by nWavs)*(nWavs by nNodeNat) = chromophore x node
                %#########################
                
                if basisFlag %In basis, so map from basis to volume, the to GM
                    hbo_tmp = img(1,:)';
                    hbr_tmp = img(2,:)';
                    hbo_tmpVol = hBasis.Map('S->M',hbo_tmp);
                    hbr_tmpVol = hBasis.Map('S->M',hbr_tmp);
                    hbo_tmpGM = (vol2gm*hbo_tmpVol);
                    hbr_tmpGM = (vol2gm*hbr_tmpVol);
                    
                    hbo.vol(frame,:,cond) = hbo_tmpVol;
                    hbr.vol(frame,:,cond) = hbr_tmpVol;
                    hbo.gm(frame,:,cond) = hbo_tmpGM;
                    hbr.gm(frame,:,cond) = hbr_tmpGM;
                else         %Not using basis
                    if strcmpi(varInputs.reconSpace,'cortex') %In GM already
                        hbo.gm(frame,:,cond) = img(1,:);
                        hbr.gm(frame,:,cond) = img(2,:);
                    else                                      %In volume, map to GM
                        hbo_tmpVol = img(1,:);
                        hbr_tmpVol = img(2,:);
                        hbo_tmpGM = (vol2gm*hbo_tmpVol');
                        hbr_tmpGM = (vol2gm*hbr_tmpVol');
                        
                        hbo.vol(frame,:,cond) = hbo_tmpVol;
                        hbr.vol(frame,:,cond) = hbr_tmpVol;
                        hbo.gm(frame,:,cond) = hbo_tmpGM;
                        hbr.gm(frame,:,cond) = hbr_tmpGM;
                    end
                end
            end
            
            if muaFlag %Calculate mua images if imageType = 'both' or 'mua'
                for wav = 1:nWavs
                    if basisFlag
                        tmp = hBasis.Map('S->M',muaImageAll(wav,:));
                        mua{wav}.vol(frame,:,cond) = tmp;
                        mua{wav}.gm(frame,:,cond) = (vol2gm*tmp)';
                    else
                        if strcmpi(varInputs.reconSpace,'cortex') %In GM already
                            tmp = muaImageAll(wav,:);
                            mua{wav}.gm(frame,:,cond) = tmp;
                        else                                      %In volume, map to GM
                            tmp = muaImageAll(wav,:);
                            mua{wav}.vol(frame,:,cond) = tmp;
                            mua{wav}.gm(frame,:,cond) = (vol2gm*tmp)';
                        end
                    end
                end
            end
        end
    end
end

if nCond==1
    hbo.vol = squeeze(hbo.vol);
    hbo.gm = squeeze(hbo.gm);
    hbr.vol = squeeze(hbr.vol);
    hbr.gm = squeeze(hbr.gm);
    if muaFlag
        for wav = 1:nWavs
            mua{wav}.vol = squeeze(mua{wav}.vol);
            mua{wav}.gm = squeeze(mua{wav}.gm);
        end
    end
end

if ~varInputs.saveVolumeImages || strcmpi(varInputs.reconSpace,'cortex') %If not saving volume, populate empty
    hbo.vol = [];
    hbr.vol = [];
    for wav = 1:nWavs
        mua{wav}.vol = [];
    end
elseif strcmpi(varInputs.imageType,'haem') %saving volume, but not mua
    for wav = 1:nWavs
        mua{wav}.vol = [];
    end
end

%################ Create dotimg structure and write .dotimg #####################
%##########################################################################
if ~isfield(prepro,'fileName')
    dotimgFileName = [];
    ds = datestr(now,'yyyymmDDHHMMSS');
    prepro.fileName = [];
else 
    if isempty(prepro.fileName)
        dotimgFileName = [];
        ds = datestr(now,'yyyymmDDHHMMSS');
        prepro.fileName = [];
    else
        [pathstr, name, ~] = fileparts(prepro.fileName);
        ds = datestr(now,'yyyymmDDHHMMSS');
        dotimgFileName = fullfile(pathstr,[name '.dotimg']);
    end
end

logData(1,:) = {'Created on: ', ds};
logData(2,:) = {'Associated prepro file: ', prepro.fileName};
logData(3,:) = {'Associated invjac file: ', invjac.fileName};
logData(4,:) = {'Associated rmap file: ', rmap.fileName};
logData(5,:) = {'reconMethod: ', varInputs.reconMethod};
logData(6,:) = {'regMethod: ', varInputs.regMethod};
logData(7,:) = {'hyperParameter: ', varInputs.hyperParameter};

[dotimg, dotimgFileName] = DOTHUB_writeDOTIMG(dotimgFileName,logData,hbo,hbr,mua,prepro.tDOD,varInputs.saveFlag);



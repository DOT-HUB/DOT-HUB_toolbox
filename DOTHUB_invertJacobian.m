function [invjac, invjacFileName] = DOTHUB_invertJacobian(jac,prepro,varargin)

% Inverts jacobian in the specified manner
% ####################### INPUTS ##########################################
% jac       =  jac structure or path to .jac. If jac.basis is not empty, a
%              toast mesh basis is assumed and rebuilt in order to create the volume and then GM
%              images. If you don't have jac, please use DOTHUB_makeToastJacobian
%
% prepro    =  prepro structure or path to .prepro. Only needed for
%              covariance regularilization
%
% rmap      =  rmap structure or path to .rmap
%
% varargin  =  optional input pairs:
%              'reconMethod' - 'multispectral' or 'standard' (default 'standard');
%                   Specifying whether to construct and invert a multispectral
%                   jacobian or whether to recontruct each wavelength
%                   separately and then combine them
%              'reconSpace' - 'full' or 'cortex' (default 'full' = volume mesh or basis);
%              'regMethod' - 'tikhonov' or 'covariance' or 'spatial' (default 'tikhonov')
%                   Regularization method.
%                   'tikonov' = standard 0th order
%                   'covariance' = exploits preproc.dod data to normalize by covariance
%                   'spatial' = spatially varying regularization as in:
%                   White, B. (2012). Developing High-Density Diffuse Optical Tomography
%                   for Neuroimaging. Washington University. pg 23, PhD Thesis.
%                   (https://openscholarship.wustl.edu/etd/665/) Accesed on 16.04.2019
%              'hyperParameter' - numerical value or vector (for 'spatial') (default 0.01);
%                   Regularization hyperparamter. See DOTHUB_invertJacobian for more details
%              'rmap' - structure or path to .rmap file.
%                   Necessary for spatially varying regularization.
%              'saveFlag' - 'true' or 'false' (default 'true');
%                   Flag whether to save invjac to disk;
%
% ####################### OUTPUTS #########################################
%
% invjac            : The invjac stucture
%
% invjacFileName    : The full pathname of the target invjac file.
%
% ####################### Dependencies ####################################
% #########################################################################
% EVR & RJC, UCL, April 2020

fprintf('################## Running DOTHUB_invertJacobian ####################\n');

% MANAGE VARIABLES
% #########################################################################
varInputs = inputParser;
varInputs.CaseSensitive = false;
validateReconMethod = @(x) assert(any(strcmpi({'standard','multispectral'},x)));
addParameter(varInputs,'reconMethod','standard',validateReconMethod);
validateSpace = @(x) assert(any(strcmpi({'volume','basis','cortex'},x)));
addParameter(varInputs,'reconSpace','full',validateSpace);
validateRegMethod = @(x) assert(any(strcmpi({'tikhonov','covariance','spatial'},x)));
addParameter(varInputs,'regMethod','tikhonov',validateRegMethod);
addParameter(varInputs,'hyperParameter',0.01,@isnumeric);
addParameter(varInputs,'rmap',[]);
validateFlag = @(x) assert(x==0 || x==1);
addParameter(varInputs,'saveFlag',true,validateFlag);
parse(varInputs,varargin{:});
varInputs = varInputs.Results;
%varInputs.Results
%      hyperParameter: 0.0100
%         reconMethod: 'standard'
%           regMethod: 'tikhonov'
%       saveMuaImages: 0
%    saveVolumeImages: 0

%More parsing and error handling

if ~isempty(varInputs.rmap) %rmap parsed
    if ischar(varInputs.rmap)
        rmapFileName = varInputs.rmap;
        rmap = load(rmapFileName,'-mat');
    else
        rmap = varInputs.rmap;
        rmapFileName = rmap.fileName;
    end
end
fnames = fieldnames(varInputs);
if strcmpi(varInputs.regMethod,'spatial')
    if length(varInputs.hyperParameter)<2
        error('For spatial regularization, the hyperParameter must be a vector');
    end
end
if strcmpi(varInputs.reconSpace,'cortex') && strcmpi(varInputs.regMethod,'spatial')
    warning('You cannot combine cortically-contrained reconstruction and spatial regularization. Reverting to tikhonov...');
    eval('varInputs.regMethod = ''tikhonov'';');
end
if (strcmpi(varInputs.reconMethod,'cortical') || strcmpi(varInputs.regMethod,'spatial')) && ~any(strcmpi(fnames,'rmap'))
    error('Both cortically-contrained reconstruction and spatial regularization require [''rmap'', rmap] as input argument pair');
end

%Display selected parameters
fprintf(['***INPUT PARAMETERS***\n'])
for i = 1:numel(fnames)
    if strcmpi(fnames{i},'rmap');continue;end
    fprintf([fnames{i} ' = ' num2str(getfield(varInputs,fnames{i})) '\n'])
end
fprintf('\n');
    
% #########################################################################
% Load jac and/prepro structures if they are parsed as paths;
if ischar(jac)
    jacFileName = jac; %Might want to force these to be the full path.
    jac = load(jacFileName,'-mat');
else
    jacFileName = jac.fileName;
end
if ~isempty(prepro)
    if ischar(prepro) %Don't necessarily need prepro
        preproFileName = prepro;
        prepro = load(preproFileName,'-mat');
    else
        preproFileName = prepro.fileName;
    end
end

% #########################################################################
% Unpack a few things...
SD3D = prepro.SD3D;
wavelengths = SD3D.Lambda;
nWavs = length(wavelengths);
hyperParameter = varInputs.hyperParameter;

if strcmpi(varInputs.reconSpace,'cortex') %Cortically constrained
    basisFlag = 0;
    invJbasis = [];
    nNodeNat = size(jac.J{1}.gm,2);
    for wav = 1:nWavs
        JNat{wav} = jac.J{wav}.gm;
    end
else  %Full space
    if ~isempty(jac.basis)  %In basis
        basisFlag = 1;
        invJbasis = jac.basis;
        nNodeNat = size(jac.J{1}.basis,2);
        for wav = 1:nWavs
            JNat{wav} = jac.J{wav}.basis;
        end
    else                    %Full mesh volume
        basisFlag = 0;
        nNodeNat = size(jac.J{1}.vol,2);
        for wav = 1:nWavs
            JNat{wav} = jac.J{wav}.vol;
        end
    end
end
    
% #########################################################################
% Perform Channel Rejection prior to inversion.
% Make sure the active channels at all wavelengths are the same!
Eall = [];
for i = 1:nWavs
    tmp = JNat{i};
    JNatCropped{i} = tmp(SD3D.MeasListAct(SD3D.MeasList(:,4)==i)==1,:);
end

% Determine reconMethod
% #########################################################################
% STANDARD
if strcmpi(varInputs.reconMethod,'standard')
    fprintf('Inverting Jacobian...\n');
    %standard reconstruction approach
    %determine regMethod
    switch varInputs.regMethod
        case 'tikhonov'
            fprintf('Running tikhonov regularized inversion...\n');
            % Inversion matrix, wavelength i
            invJ = cell(nWavs,1);
            for i = 1:nWavs
                Jtmp = JNatCropped{i};
                JJT = Jtmp*Jtmp';
                S=svd(JJT);
                invJ{i} = Jtmp'/(JJT + eye(length(JJT))*(varInputs.hyperParameter*max(S)));
            end
        case 'covariance'
            fprintf('Covariance reconstruction (beta!)\n');
            invJ = cell(nWavs,1);
            for i = 1:nWavs
                Jtmp = JNatCropped{i};
                %Extract relevant chunk of prepro.dod to calculate covariance. %Might be better to pass the reference data directly?
                if length(find(prepro.tDOD<0))>1  %Must be an HRF
                    [~,ind] = min(abs(prepro.tDOD));
                    covData = prepro.dod(1:ind,:);
                else %take first 30 seconds, or 10% of data, whichever is less.
                    if DOTHUB_range(tDOD)< 30
                        covData = prepro.dod;
                    else
                        fs = 1/mean(diff(tDOD));
                        covData = prepro.dod(1:round(fs*30),:);
                    end
                end   
                sigma_v = cov(covData); 
                sigma_u = sparse(1:2*nNodeNat,1:2*nNodeNat,1);
                JJT = Jtmp*sigma_u*Jtmp';
                l1 = hyperParameter*trace(JJT)/trace(sigma_v);
                invJ{i} = sigma_u*Jtmp' / ( JJT + sigma_v*l1);
            end
        case 'spatial'
            fprintf('Running spatial regularized inversion (beta!) \n');
            l1 = hyperParameter(1); % Typical Tikhonov regularization parameter
            l2 = hyperParameter(2); % Spatial regularization parameter
            invJ = cell(nWavs,1);
            for i = 1:nWavs
                Jtmp = JNatCropped{i};

                JJT = Jtmp*Jtmp';  % Prepare Jacobian matrix for inversion (i.e. create square matrix)
                L = sqrt(diag(JJT) + l2*max(diag(JJT))); % Apply regularization
                Linv = 1./L; % Invert matrix

                % Find Atild
                Atild = zeros(size(Jtmp));
                for ind = 1:length(Linv)
                    Atild(ind,:) = Jtmp(ind,:)*Linv(ind);
                end
                atildtatild = Atild*Atild';
                [satild] = svd(atildtatild);
                mxsatild = max(satild);

                % Apply spatial regularization
                val2binv = atildtatild;
                for ind = 1:length(satild)
                    val2binv(ind,ind) = atildtatild(ind,ind) + l1*mxsatild;
                end
                clear JJT   % Clear for efficiency

                % Invert matrix
                inva = val2binv\Atild;
                clear Atild val2binv  % Clear huge matrices for efficiency
                invJtmp = zeros(size(inva))';
                for ind = 1:size(inva,1)
                    invJtmp(:,ind) = inva(ind,:)*Linv(ind);
                end
                invJ{i} = invJtmp;
            end
            
    end
    clear JJT   % Clear for efficiency
end

% #########################################################################
% MULTISPEC
if strcmpi(varInputs.reconMethod,'multispectral')
    fprintf('Building Multispectral Jacobian...\n');
    %Use loop to get specific absorption coefficients
    Eall = [];
    Jtiled = [];
    for i = 1:nWavs
        Etmp = GetExtinctions(wavelengths(i));
        Eall = [Eall; Etmp./1e7]; %Combine specific absorption coeffs into matrix (wavelength x chromphore), convert units from cm-1/M to mm-1/uM
        Jtiled = [Jtiled; JNatCropped{i} JNatCropped{i}]; %Tile wavelength-specific jacobians
    end
    
    %Building extinction coefficient dummy
    Ei = ones(size(JNatCropped{1}));
    Etiled = [];
    for c = 1:2 %Chromophore
        El = [];
        for i = 1:length(wavelengths)
            El = [El; Ei.*Eall(i,c)];
        end
        Etiled = [Etiled El];
    end
    
    %Build multispectral Jacobian.
    Jmulti = Jtiled.*Etiled; % This has units of (d(ln(data/reference))/d(mm-1))*(mm-1/micromolar)
    
    fprintf('Inverting Jacobian...\n');
    %determine regMethod
    switch varInputs.regMethod
        case 'tikhonov'
            fprintf('Running tikhonov regularized inversion...\n');
            JJT = Jmulti*Jmulti';
            S=svd(JJT);
            invJ{1} = Jmulti' / (JJT + eye(length(JJT))*(hyperParameter*max(S)));
            
        case 'covariance'
            fprintf('Running covariance regularized inversion (beta!)...\n');
            %Extract relevant chunk of prepro.dod to calculate covariance. %Might be better to pass the reference data directly?
            if length(find(prepro.tDOD<0))>1  %Must be an HRF
                [~,ind] = min(abs(prepro.tDOD));
                covData = prepro.dod(1:ind,:);
            else %take first 30 seconds, or 10% of data, whichever is less.
                if DOTHUB_range(tDOD)< 30
                    covData = prepro.dod;
                else
                    fs = 1/mean(diff(tDOD));
                    covData = prepro.dod(1:round(fs*30),:);
                end
            end   
            sigma_v = cov(covData); 
            sigma_u = sparse(1:2*nNodeNat,1:2*nNodeNat,1);
            JJT = Jmulti*sigma_u*Jmulti';
            l1 = hyperParameter*trace(JJT)/trace(sigma_v);
            invJ{1} = sigma_u*Jmulti' / ( JJT + sigma_v*l1);
            
        case 'spatial'
            fprintf('Running spatial regularized inversion (beta!) \n');
            
            l1 = hyperParameter(1); % Typical Tikhonov regularization parameter
            l2 = hyperParameter(2); % Spatial regularization parameter
            
            JJT = Jmulti*Jmulti';  % Prepare Jacobian matrix for inversion (i.e. create square matrix)
            L = sqrt(diag(JJT) + l2*max(diag(JJT))); % Apply regularization
            Linv = 1./L; % Invert matrix
            
            % Find Atild
            Atild = zeros(size(Jmulti));
            for ind = 1:length(Linv)
                Atild(ind,:) = Jmulti(ind,:)*Linv(ind);
            end
            atildtatild = Atild*Atild';
            [satild] = svd(atildtatild);
            mxsatild = max(satild);
            
            % Apply spatial regularization
            val2binv = atildtatild;
            for ind = 1:length(satild)
                val2binv(ind,ind) = atildtatild(ind,ind) + l1*mxsatild;
            end
            clear JJT;   % Clear for efficiency
            
            % Invert matrix
            inva = val2binv\Atild;
            clear Atild val2binv  % Clear huge matrices for efficiency
            invJtmp = zeros(size(inva))';
            for ind = 1:size(inva,1)
                invJtmp(:,ind) = inva(ind,:)*Linv(ind);
            end
            invJ{1} = invJtmp;
            clear invJtmp;
    end
end
clear JJT;   % Clear for efficiency

% #########################################################################
% Create invjac filename. We want this to be identical to the jacFileName  
[pathstr, name, ~] = fileparts(jac.fileName);
invjacFileName = fullfile(pathstr,[name '.invjac']);
ds = datestr(now,'yyyymmDDHHMMSS');
logData(1,:) = {'Created on: ', ds};
logData(2,:) = {'Derived from jac file: ', jac.fileName};
%Include varInputs
fnames = fieldnames(varInputs);
for i = 1:length(fnames)
    if strcmpi(fnames{i},'rmap');continue; end
    logData(end+1,:) = {fnames{i}, getfield(varInputs,fnames{i})};
end
[invjac, invjacFileName] = DOTHUB_writeINVJAC(invjacFileName,logData,invJ,invJbasis,varInputs.saveFlag);


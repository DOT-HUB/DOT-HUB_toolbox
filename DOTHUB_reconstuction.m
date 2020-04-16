function [hbo_image_GM, hbr_image_GM, hbo_image_vol, hbr_image_vol,hbo_image_basis,hbr_image_basis,...
          J_basis_wav1, J_basis_wav2, J_GM_wav1, J_GM_wav2, J_vol_wav1, J_vol_wav2] = LUMOreconstruction(data_OD,ref_OD,SD,source_pos,detector_pos,headMesh,gmMesh,lambda,varargin)
% LUMOreconstruction.m
% Linear reconstruction of concentration changes

% INPUTS ##################################################################
% data_OD   =   vector of OD values of size (time x channel).  Number of
%               channels should match the jacobians.  This is equivalent to 
%               n time points (frames) of the homer2 .nirs data vector dod.

% ref_OD    =   vector of OD values of size (timeref x channel).  This is the baseline 
%               against which data_OD will be compared.  If you are reconstructing an HRF
%               this should be the vector of average data points prior to
%               stimulus onset. If you are reconstructing a continuous resting state data
%               series, ref_OD should be parsed empty [], and images will
%               be reconstructed relative to the mean of the data_OD (which
%               is usually zero).  If you are reconstructing a whole data
%               series that DOES have stimulations in is, then ref_OD
%               should probabaly be the data from the initial rest period of
%               the whole timecourse. The number of channels should match the 
%               jacobians.

% SD        =   The SD file associated with the parsed data.  Only used for
%               wavelengths here.

% headMesh  =   The tetrahedral head model mesh structure with fields .node 
%               (nnodes x 4) with columns [x,y,z,nodal tissue type flag] and 
%               .elem (nelem x 4)

% gmMesh    =   The triangular GM surface mesh structure with fields .node 
%               (nnodes x 3) with columns [x,y,z] and .face (nface x 3)

% vol2gm    =   The mapping matrix from the volume mesh to the GM mesh, with
%               dimensions NxM where M is the number of tetrahedral mesh nodes 
%               and N is the number of GM surface nodes Can be calculated 
%               using the Vol2GM_Transform function. If not parsed,
%               will be calculated

% lambda    =   The regularization paramter coefficient.  Default is 0.01.

% path_results (Optional) = Define a path to save Jacobian and inverted
%               Jacobian to speed up future reconstructions, provided that
%               the only change is the input data (for example, several 
%               trials with the same subject)

% OUTPUTS (RECONSTRUCTION RESULTS) ########################################

% hbo_image_GM =    the change in hbo, projected to the cortical surface.
%                   Will have dimensions time x number of gmMesh nodes

% hbr_image_GM =    the change in hbr, projected to the cortical surface.
%                   Will have dimensions time x number of gmMesh nodes

% hbo_image_vol =   the change in hbo in the volume mesh space. time x
%                   number of headMesh nodes

% hbr_image_vol =   the change in hbr in the volume mesh space. time x
%                   number of headMesh nodes

% hbo_image_basis = the change in hbo in the basis space. time x
%                   number of valid basis voxels

% hbr_image_basis = the change in hbr in the basis space. time x
%                   number of valid basis voxels

% OUTPUTS (JACOBIAN SECTION) ##############################################

% J_basis_wav1 =    The jacobian for the shorter wavelength (specified in
%                   SD), in basis space

% J_basis_wav2 =    The jacobian for the longer wavelength (specified in
%                   SD), in basis space

% J_GM_wav1 =       The jacobian for the shorter wavelength mapped on to
%                   the GM surface

% J_GM_wav2 =       The jacobian for the shorter wavelength mapped on to
%                   the GM surface

% J_vol_wav1 =      The jacobian for the shorter wavelength in the full
%                   mesh volume space

% J_vol_wav2 =      The jacobian for the longer wavelength in the full
%                   mesh volume space

%% Declare extra arguments (gray-matter surface and results folder)

if exist('vol2gm','var')~=1 % (Necessary to display reconstruction on the gray matter surface)
    vol2gm = []; % If not provided, it will be during the process
end

if exist('path_results','var')~=1  %(Optional)
    path_results = []; % Reconstructions will be saved in the path provided
end

if exist('path_results','var')~=1  %(Optional)
    path_Jacobian = []; % Reconstructions will be saved in the path provided
end

if length(varargin) >= 1
    if isnumeric(varargin{1}) && ~isempty(varargin{1})
        vol2gm = varargin{1};
    else
        disp(' Transformation matrix from Head model to gray matter surface not provided, it will be calculated later')
    end   
end
if length(varargin) >= 2
    if ischar(varargin{2})
        path_results = varargin{2};
        [folder,filename,~] = fileparts(path_results);
        if exist(folder,'dir')~=7
            mkdir(folder)
        end
        % if ~isempty(filename)
        %     if exist([folder,filename,'.mat'],'file') == 7
        %         warning('A file with that name already exists')
        %         opc = input('Are you sure you want to continue (type Y or N)?','s');
        %         if opc == 'N'
        %             error('program aborted')
        %         end
        %     end
        % else
        %     D = dir(folder);
        %     if length(D)>2
        %         f = warndlg('The folder provided is not empty and possible contains reconstruction data');
        %         opc = input('Are you sure you want to continue (type Y or N)?','s');
        %         if opc == 'N'
        %             error('program aborted')
        %         end
        %     end
        % end
        
    else
        warning('Your data won''t be saved but it will be in the workspace as a variable')
    end
end

% Save Jacobian
if length(varargin) >= 3
    if ischar(varargin{3})
        path_Jacobian = varargin{3};
    else
        path_Jacobian = [];
        warning('The Jacobian won''t be saved but it will be in the workspace as a variable')
    end    
end

% Standard of multispectral reconstruction
RecMethodStr = 'multi-spectral';
RecMethod = 'multispec';
if length(varargin) >= 4
    if ischar(varargin{4})
        RecMethod = varargin{4};
        switch RecMethod 
            case 'multispec', RecMethodStr = 'multi-spectral';
            case 'stnd', RecMethodStr = 'standard';
            otherwise, RecMethodStr = 'multi-spectral';
                warning('Unknown reconstruction method')
        end
    else
        warning('Default reconstruction method')
    end    
end 
warning(['Reconstruction method: ',RecMethodStr])
%% Jacobian calculation
if exist([path_Jacobian,'Jacobian.mat'],'file') ~= 2
    [J_basis_wav1, J_basis_wav2, J_GM_wav1, J_GM_wav2, J_vol_wav1, J_vol_wav2,headMesh] = DOTHUB_Make_ToastJacobian_v01(source_pos,detector_pos,SD,headMesh,gmMesh,vol2gm);
    if ~isempty(path_Jacobian)
        [folder,~,~] = fileparts(path_Jacobian);
        if exist(folder,'dir') ~= 7
            warning([folder '  DOES NOT EXIST!, please provide a path'])
            mesh_pathname = uigetdir('*.mat','Select the folder and name to save Jacobian...');
            mesh_filename = 'Jacobian';
            save([mesh_pathname mesh_filename],'J_basis_wav1','J_basis_wav2','J_GM_wav1','J_GM_wav2','J_vol_wav1','J_vol_wav2','headMesh','-v7.3')        
        else
            save([path_Jacobian,'Jacobian'],'J_basis_wav1','J_basis_wav2','J_GM_wav1','J_GM_wav2','J_vol_wav1','J_vol_wav2','headMesh','-v7.3')        
        end
        
    end
else
    load([path_Jacobian,'Jacobian.mat'],'J_basis_wav1','J_basis_wav2','headMesh')
    %load([path_results,'Jacobian.mat'],'headMesh')
end

% Make toast-handle for mesh ##############################################
% It is assumed that the loaded mesh corresponds to the mesh which has been
% corrected if any anomalities were found, such as missing nodes or wrong
% vertex ordering (reflected in negative volumes). Such a mesh has been
% saved together with the Jacobian and uploaded for reconstruction
fprintf('Building TOAST mesh\n');
eltp = ones(length(headMesh.elem),1)*3;
hMesh = toastMesh(headMesh.node(:,1:3),headMesh.elem(:,1:4),eltp);

%Select basis
basis = [50 50 50]; %This must match the basis in which the jacobian was calculated
fineBasis = basis.*2;
% Set basis
fprintf('Building TOAST basis\n');
hBasis = toastBasis(hMesh,basis,fineBasis);
%% Part 1. Account for list of active channels (reject bad channels)
% Part 2. Since, not necessarily all the channels from the Jacobian will be
%         of good quality, the Jacobian prciously calculate will be trimmed
%         for bad channels before inversion

% Make sure the active channels at both wavelenths are the same!
SD.MeasListAct(1:end/2) = SD.MeasListAct(1:end/2) & SD.MeasListAct(end/2+1:end);
SD.MeasListAct(end/2+1:end) = SD.MeasListAct(1:end/2) & SD.MeasListAct(end/2+1:end);

J_basis_wav1_cropped = J_basis_wav1(SD.MeasListAct(1:end/2)==1,:);
J_basis_wav2_cropped = J_basis_wav2(SD.MeasListAct(end/2+1:end)==1,:);
% Only ACTIVE CHANNELS will be used for inversion (Assuming only two wavelentgh)
invJ = DOTHUB_inversion_matrix(hMesh,hBasis,SD,J_basis_wav1_cropped,...
                                               J_basis_wav2_cropped,headMesh,lambda,RecMethod);

%% Set data in TOAST format
% Convert data into toast style (toast wants = ln(Intensity_active)-ln(intensity_baseline)
% Parsed data is OD (i.e. data_OD = -ln(intensity_active/mean), ref_OD = -ln(intensity_baseline/mean)
% So toast wants = -(data_OD - ref_OD);

if isempty(ref_OD)     % Reconstruct relative to mean, use ref_OD
    ref_OD = data_OD;
end
if size(ref_OD,1) == 1 % Ref data is a single frame (note that this won't work for the Covar norm version)
    data_to_recon = -(data_OD-repmat(ref_OD,size(data_OD,1),1));
else
    data_to_recon = -(data_OD-repmat(mean(ref_OD),size(data_OD,1),1));    
end

% Getting rid of bad channels (could be more elegant! try again!)
data_to_recon = data_to_recon(:,SD.MeasListAct==1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create matrices to save recontruction results
nNodesBasis = size(invJ,1)/2;
nNodesGM = size(vol2gm,1);
nFrame = size(data_OD,1);

hbo_image_basis = zeros(nFrame,nNodesBasis);
hbr_image_basis = zeros(nFrame,nNodesBasis);
hbo_image_GM = zeros(nFrame,nNodesGM);
hbr_image_GM = zeros(nFrame,nNodesGM);

disp('Reconstruction ##############################################');
for frame = 1:nFrame
    fprintf('Reconstructing frame %d of %d\n',frame,nFrame);

    data = data_to_recon(frame,:);
    img = invJ * data';

    hbo_image_basis(frame,:) = img(1:nNodesBasis);
    hbr_image_basis(frame,:) = img(nNodesBasis+1:2*nNodesBasis);

    hbo_image_vol = hBasis.Map('S->M',hbo_image_basis(frame,:)');
    hbr_image_vol = hBasis.Map('S->M',hbr_image_basis(frame,:)');

    hbo_image_GM(frame,:) = (vol2gm*hbo_image_vol(:))';
    hbr_image_GM(frame,:) = (vol2gm*hbr_image_vol(:))';

end

disp('Reconstructions complete!');
[folder,filename,~] = fileparts(path_results);
if isempty(filename)
    filename = 'reconstruction';
end
save([folder,filesep,filename],'hbo_image_GM','hbr_image_GM','hbo_image_vol','hbr_image_vol','hbo_image_basis','hbr_image_basis');
function [J_basis_wav1, J_basis_wav2, J_GM_wav1, J_GM_wav2, J_vol_wav1, J_vol_wav2,headMesh] = DOTHUB_Make_ToastJacobian_v01(source_pos,detector_pos,SD,headMesh,gmMesh,varargin)

%Calculates Jacobian using Toast++.

% #########################################################################
% INPUTS ##################################################################

% source_pos =      Positions of source optodes in the headMesh space

% detector_pos =    Positions of detector optodes in the headMesh space

% SD =              The SD file for the array used in the study.  Note that
%                   the srcpos and detpos form the SD file are ignored in
%                   favour of the above defined inputs

% headMesh =        The tetrahedral head model mesh structure with fields .node 
%                   (nnodes x 4) with columns [x,y,z,nodal tissue type flag] and 
%                   .elem (nelem x 4)

% gmMesh =          The triangular GM surface mesh structure with fields .node 
%                   (nnodes x 3) with columns [x,y,z] and .face (nface x 3)

% vol2gm =          The mapping matrix from the volume mesh to the GM mesh, with
%                   dimensions NxM where M is the number of tetrahedral mesh nodes 
%                   and N is the number of GM surface nodes Can be calculated 
%                   using the Vol2GM_Transform function. If not parsed,
%                   will be calculated


% OUTPUTS #################################################################

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

% #########################################################################
% #########################################################################
% Version 0.1
% RJC, University College London, May 2019
% EVR, University College London, Modified on June 2019
% #########################################################################
% #########################################################################



%Flag 1 to run initial single-channel calculation test before running full
%calculation (this is useful because the test should only take a few
%minutes, so if it hangs for longer, you can be certain there is a problem,
%whereas the full calculations can take hours, so you can't be sure it is
%working.
testFlag = 1;

%Select basis
basis = [50 50 50];
fineBasis = basis.*2;
%Wavelengths
wavelengths = SD.Lambda;

% Check for errors in mesh and correct  ###################################
% #########################################################################
% Correct negatives if they exist
if any(headMesh.node(:,4)<0)
    headMesh = CreateNodalTissueSpec_SB_fastversion(headMesh);
end

elem_tmp = headMesh.elem(:,1:4);
node_tmp = headMesh.node(:,1:4);% Original: node_tmp = headMesh.node(:,1:3); % correted on 15th April 2020
included_nodes = ismember(1:length(headMesh.node),elem_tmp(:));
errnodes_ind = find(~included_nodes);
if ~isempty(errnodes_ind)
    fprintf('Erroneous nodes found, correcting mesh...\n');
    %correct node list
    node_corr = node_tmp(included_nodes,:);
    %correct element list
    for i = 1:length(errnodes_ind)
        elem_tmp(elem_tmp > (errnodes_ind(i)-(i-1))) = elem_tmp(elem_tmp > (errnodes_ind(i)-(i-1)))-1;
    end
    
    included_nodes = ismember(1:length(node_corr),elem_tmp(:));
    errnodes_ind = find(~included_nodes,1);
    if isempty(errnodes_ind)
        fprintf('Erroneous nodes removed\n');
        headMesh.node = node_corr;
        headMesh.elem(:,1:4) = elem_tmp;
        clear elem_tmp
        clear node_corr
        clear included_nodes
        clear errnodes_ind
        
        %Save corrected mesh
        [mesh_filename, mesh_pathname] = uigetfile('*.mat','Select mesh file to overwrite with corrected mesh...');
        mesh_tmp = load([mesh_pathname mesh_filename],'-mat');
        mesh_fields = fields(mesh_tmp); clear mesh_tmp
        fprintf('Saving corrected mesh...\n');
        eval([mesh_fields{1} ' = headMesh;']);
        eval('save([mesh_pathname mesh_filename],mesh_fields{1})');
    else
        fprintf('Erroneous nodes remain, ABORTING');
        return
    end
end

% #########################################################################
% Make toast-handle for mesh ##############################################
fprintf('Building TOAST mesh\n');
eltp = ones(length(headMesh.elem),1)*3;
headMesh.elem(:,1:4) = headMesh.elem(:,[4 1:3]);
hMesh = toastMesh(headMesh.node(:,1:3),headMesh.elem(:,1:4),eltp);

%Check mesh is ordered as TOAST prescribes
fprintf('Checking mesh...\n');


if any(hMesh.ElementSize()<0)
    fprintf('Negative volume elements in HD mesh, attempting to reconfigure...\n');
    headMesh.elem(:,1:4) = headMesh.elem(:,[2 3 4 1]);
    hMesh = toastMesh(headMesh.node(:,1:3),headMesh.elem(:,1:4),eltp); 
    
    
    if any(hMesh.ElementSize()<0)
        fprintf('Negative volume elements remain, ABORTING');
        return
    end
end

%Force source_pos and detector_pos to nearestNode in headMesh
for i = size(detector_pos,1):-1:1
    nnode = Nearest_node_RJC(detector_pos(i,:),headMesh.node(:,1:3));
    det_pos(i,:) = nnode;
end
for i = size(source_pos,1):-1:1
    nnode = Nearest_node_RJC(source_pos(i,:),headMesh.node(:,1:3));
    source_pos(i,:) = nnode;
end

%Define linklist
linklist = SD2linklist_RJC(SD);
fprintf('Building QM file...\n');
qmoutname = 'tmp.qm';
writeqm_RJC(qmoutname,det_pos,source_pos,linklist);

% #########################################################################
% Generate Jacobian #######################################################
% Generate optical property vectors on the nodes
% Optical properties
c0 = 0.3;   % speed of light in vacuum
refind = 1.3; % homogeneous refractive index
refindvec = ones(length(headMesh.node),1)*refind;
c_medium = c0./refindvec;

% Set optical properties
[lambda,mua_model,mus_model,tissues] = TissuePropertiesModel;
muavec = zeros(length(wavelengths),length(headMesh.node));
musvec = zeros(length(wavelengths),length(headMesh.node));
headMesh.node(headMesh.node(:,4)>5,4) = 5;
% How many tissues are there?
nTissues = unique(headMesh.node(:,4))';
if max(nTissues) == 5
    %Assume Scalp, Skull, CSF, GM, WM
    for wav_count = 1:length(wavelengths)
        
        [~,wav_ind] = min(abs(lambda - wavelengths(wav_count)));
        
        for i = 1:5 %Tissue types CSF, GM, WM,
            muavec(wav_count,headMesh.node(:,4)==i) = mua_model(wav_ind,i);
            musvec(wav_count,headMesh.node(:,4)==i) = mus_model(wav_ind,i);
        end
    end
else
    if max(nTissues) == 4
        %Assume ECT, CSF, GM, WM
        for wav_count = 1:length(wavelengths)
            
            [~,wav_ind] = min(abs(lambda - wavelengths(wav_count)));
            
            muavec(wav_count,headMesh.node(:,4)==1) = mean(mua_model(wav_ind,1:2),2); %ECT = mean(scalp+skull)
            musvec(wav_count,headMesh.node(:,4)==1) = mean(mus_model(wav_ind,1:2),2);
            for i = 2:4 %Tissue types CSF, GM, WM,
                muavec(wav_count,headMesh.node(:,4)==i) = mua_model(wav_ind,i+1);
                musvec(wav_count,headMesh.node(:,4)==i) = mus_model(wav_ind,i+1);
            end
        end
    else
        error('Nodal tissue indices are not as expect for infant (1:4) or adult (1:5)');
    end
end
muavec = muavec';
musvec = musvec';

%Jacobian inputs
jtype='bicgstab';
bicgstabtol=1e-12;

%Test how long this will take using single channel
% #########################################################################
if testFlag == 1
    disp('Running test calculation...');
    writeqm_RJC('Test.qm',det_pos(1,:),source_pos(1,:),1); % Modif. on 08/April/20: by default linklist(1,1) was used, but crashed if different from 1
    hMesh.ReadQM('Test.qm');
    % Set source locations
    srcPosMesh = hMesh.Qpos;
    % Set position locations
    detPosMesh = hMesh.Mpos;
    hMesh.Display;
    hold on;
    plotmesh(srcPosMesh,'r.','MarkerSize',40);
    plotmesh(detPosMesh,'b.','MarkerSize',40);
    
    qvec = hMesh.Qvec ('Neumann', 'Gaussian', 2);
    mvec = hMesh.Mvec ('Gaussian', 2,refindvec);
    % Set basis
    hBasis = toastBasis(hMesh,basis,fineBasis);
    
    if hMesh.isvalid
        fprintf('Running toastJacobianCW test\n');
        tic;
        Jtest = toastJacobianCW(hMesh, hBasis, qvec, mvec, muavec(:,1), musvec(:,1), refindvec, jtype, bicgstabtol); %#ok<NASGU>
        duration = toc;
    end
    disp(['Time estimate = ' num2str(duration*(size(det_pos,1) + size(source_pos,1))/60) ' mins']);
end

% #########################################################################
% Calculate full Jacobiab
%Re-set parameters
hMesh.ReadQM(qmoutname)
srcPosMesh = hMesh.Qpos;
detPosMesh = hMesh.Mpos;
figure
hMesh.Display;
hold on;
plotmesh(srcPosMesh,'r.','MarkerSize',40);
plotmesh(detPosMesh,'b.','MarkerSize',40);
qvec = hMesh.Qvec ('Neumann', 'Gaussian', 2); 
mvec = hMesh.Mvec ('Gaussian',2,refindvec); 

% Set basis
disp('Setting basis...');    
hBasis = toastBasis(hMesh,basis,fineBasis);
    
disp('Producing Jacobian 1...');
tic
J_basis_1 = toastJacobianCW(hMesh, hBasis, qvec, mvec, muavec(:,1), musvec(:,1), refindvec,jtype,bicgstabtol);
toc
disp('Producing Jacobian 2...');
tic
J_basis_2 = toastJacobianCW(hMesh, hBasis, qvec, mvec, muavec(:,2), musvec(:,2), refindvec,jtype,bicgstabtol);
toc

%Multiple by speed of light
% J_basis_wav1 = J_basis_1.*repmat(hBasis.Map('M->S',c_medium),1,size(J_basis_wav1,1))';
% J_basis_wav2 = J_basis_2.*repmat(hBasis.Map('M->S',c_medium),1,size(J_basis_wav2,1))';
  J_basis_wav1 = J_basis_1.*repmat(hBasis.Map('M->S',c_medium),1,size(J_basis_1,1))';
  J_basis_wav2 = J_basis_2.*repmat(hBasis.Map('M->S',c_medium),1,size(J_basis_2,1))';

disp('Mapping basis Jacobian to volume');
%Map back to volume
J_vol_wav1 = zeros(size(J_basis_wav1,1),length(headMesh.node));
J_vol_wav2 = zeros(size(J_basis_wav1,1),length(headMesh.node));
for i = 1:size(J_basis_wav1,1)
    J_vol_wav1(i,:) = hBasis.Map('S->M',J_basis_wav1(i,:));
    J_vol_wav2(i,:) = hBasis.Map('S->M',J_basis_wav2(i,:));
end

if length(varargin)>=1 % If map from volume to mesh is not provided, then calculate it
    vol2gm = varargin{1};
    if isempty(vol2gm)
        clear vol2gm
    end
end
if ~exist('vol2gm','var')
    disp('vol2gm not provided: building mapping from tetrahedral to surface...\n');
    %default radius of spherical kernal = 3 mm
    radius = 3;
    vol2gm = Vol2GM_Transform(headMesh.node,gmMesh.node,radius,1);
end

disp('Mapping volue Jacobian to GM (for visualization only)');
J_GM_wav1 = (vol2gm*J_vol_wav1')';
J_GM_wav2 = (vol2gm*J_vol_wav2')';

disp('Complete!');




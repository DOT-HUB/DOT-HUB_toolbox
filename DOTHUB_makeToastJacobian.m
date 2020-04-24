function [jac, jacFileName] = DOTHUB_makeToastJacobian(rmapFileName,basis)

%Calculates Jacobian using Toast++.

% #########################################################################
% INPUTS ##################################################################

% rmapFileName  : The full path to the rmap file containing variables:

                    % SD_3Dmesh          :   The SD structure containing registered 3D optode
                    %                        positions on the mesh, and (critically) SD.MeasList

                    % headVolumeMesh     :   The multi-layer volume mesh structure, registered
                    %                        to the relevant individual. Contains fields:
                    %                        node, face, elem, tissue_labels

                    % gmSurfaceMesh      :   The gm surface mesh structure, registered
                    %                           to the relevant individual. Contains fields:
                    %                           node, face.

                    % scalpSurfaceMesh   :   The scalp surface mesh structure, registered
                    %                           to the relevant individual. Contains fields:
                    %                           node, face.

                    % vol2gm                :   The sparse matrix mapping from head volume mesh
                    %                           space to GM surface mesh space

% basis        :   (Optiona) 1x3 vector specifying basis dimensions if desired. 
%                  A basis of [50 50 50] is aassigned by default if
%                  nnodes>200k.
            
% OUTPUTS #################################################################

% jac                      :  Structure containing all data inputs

% jacFileName              :  The full path of the resulting .jac file

% jacFileName.jac file     :  File containing all data inputs

% #########################################################################
% #########################################################################
% RJC, UCL, May 2019
% EVR, UCL, Modified on June 2019
% RJC, UCL, Modified for Github first commit
% #########################################################################

% MANAGE VARIABLES
% #########################################################################

if ~exist(basis,'var')
    if nNode>2e5 %HARD CODE NODE LIMIT AT 200,000
        basisFlag = 1;
        basis = [50 50 50];
        fineBasis = basis.*2;
    else
        basis = []; %No basis
    end
elseif basis==0
    basis = [];
end
   
%testFlag==1 to run initial single-channel calculation test before running full
%calculation. This is useful because the test should only take a few
%minutes, so if it hangs for longer, you can be certain there is a problem,
%whereas the full calculations can take hours, so you can't be sure it is
%working.
testFlag = 1;

% #########################################################################
% Load rmap ###############################################################
[rmapPath, rmapName,~] = fileparts(rmapFileName);
if isempty(rmapPath) %Ensure full path.
    rmapPath = pwd;
    rmapFileName = fullfile(rmapPath,[rmapName '.rmap']);
end
rmap = load(rmapFileName,'-mat'); %head, gm, scalp, SD_3Dmesh, vol2gm, logData

% Unpack key fields #######################################################
headVolumeMesh = rmap.headVolumeMesh;
SD_3Dmesh = rmap.SD_3Dmesh;
wavelengths = SD_3Dmesh.Lambda;

% #########################################################################
% MESH PREP ###############################################################
%Define basis if deemed necessary #########################################
nElem = size(headVolumeMesh.elem,1);
nNode = size(headVolumeMesh.node,1);
nWavs = length(wavelengths);

% Check for errors in mesh and correct if necessary
% Correct negatives if they exist in tissue indices
rewriteRMAP = 0;
if any(headVolumeMesh.node(:,4)<0)
    warning('Nodal tissue indices incorrect, attempting to fix...\n');
    headVolumeMesh = DOTHUB_createNodalTissueInd(headVolumeMesh);
    rewriteRMAP = 1;
end

% Check for erroneous nodes 
elem_tmp = headVolumeMesh.elem(:,1:4);
node_tmp = headVolumeMesh.node(:,1:3);
included_nodes = ismember(1:nNode,elem_tmp(:));
errnodes_ind = find(~included_nodes);
if ~isempty(errnodes_ind)
    warning('Erroneous nodes found, attempting to correct mesh...\n');
    %correct node list
    node_corr = node_tmp(included_nodes,:);
    %correct element list
    for i = 1:length(errnodes_ind)
        elem_tmp(elem_tmp > (errnodes_ind(i)-(i-1))) = elem_tmp(elem_tmp > (errnodes_ind(i)-(i-1)))-1;
    end
    
    included_nodes = ismember(1:length(node_corr),elem_tmp(:));
    errnodes_ind = find(~included_nodes,1);
    if isempty(errnodes_ind)
        fprintf('Erroneous nodes removed...\n');
        headVolumeMesh.node = node_corr;
        headVolumeMesh.elem(:,1:4) = elem_tmp;
        clear elem_tmp node_tmp node_corr included_nodes errnodes_ind
        rewriteRMAP = 1;
    else
        error('Correction failed: erroneous nodes remain in rmap.headVolumeMesh');
    end
end

% Make toast-handle for mesh and check configuration
fprintf('Building TOAST mesh\n');
eltp = ones(length(headVolumeMesh.elem),1)*3;
hMesh = toastMesh(headVolumeMesh.node(:,1:3),headVolumeMesh.elem(:,1:4),eltp);

%Check mesh is ordered as TOAST prescribes
%Would be nbetter if we built toast-friendly meshes from the start (*)
fprintf('Checking TOAST mesh...\n');
if any(hMesh.ElementSize()<0)
    warning('Negative volume elements in HD mesh, attempting to reconfigure...\n');
    headVolumeMesh.elem(:,1:4) = headVolumeMesh.elem(:,[4 1:3]); %(* Add this line to meshing?)
    %headVolumeMesh.elem(:,1:4) = headVolumeMesh.elem(:,[2 3 4 1]);
    hMesh = toastMesh(headVolumeMesh.node(:,1:3),headVolumeMesh.elem(:,1:4),eltp); 
    
    if any(hMesh.ElementSize()<0)
        error('Negative volume elements remain in rmap.headVolumeMesh');
    else
        fprintf('Reconfiguration successful...\n');
    end
    rewriteRMAP = 1;
end

%Save corrected mesh by overwriting rmap
if rewriteRMAP
    fprintf('Updating .rmap file with corrected headVolumeMesh \n');
    rmap.logData(end+1,:) = {'File updated by DOTHUB_makeToastJacobian on ',datestr(now,'yyyymmDDHHMMSS')};
    DOTHUB_writeRMAP(rmapFileName,rmap.logData,rmap.SD_3Dmesh,headVolumeMesh,rmap.gmSurfaceMesh,rmap.scalpSurfaceMesh,rmap.vol2gm);
end

% #########################################################################
% Set optical properties
% First assign
muaVec = zeros(length(wavelengths),length(headVolumeMesh.node));
musPrimeVec = zeros(length(wavelengths),length(headVolumeMesh.node));
refIndVec = zeros(length(wavelengths),length(headVolumeMesh.node));

% Determine tissue optical properties and populate vectors
nTissues = length(headVolumeMesh.tissue_labels);
for i = 1:nTissues
    tmpInd = find(strcmpi({'scalp','skull','ECT','CSF','GM','WM'},headVolumeMesh.tissue_labels{i}), 1);
    if isempty(tmpInd)
        error('Unknown tissue label, please correct rmap.headVolumeMesh.tissue_labels abd try again');
    end
    tissueNodeList = headVolumeMesh.node(:,4)==i;
    for j = 1:nWavs
        [mua, musPrime, refInd] = DOTHUB_getTissueCoeffs(headVolumeMesh.tissue_labels{i},wavelengths(j));%Potentially update this to be age-specific?
        muaVec(j,tissueNodeList) = mua;
        musPrimeVec(j,tissueNodeList) = musPrime;
        refIndVec(j,tissueNodeList) = refInd;
        
        opticalPropertiesByTissue(i,j,1) = mua;
        opticalPropertiesByTissue(i,j,2) = musPrime;
        opticalPropertiesByTissue(i,j,3) = refInd;
    end
end

if any([muaVec(:)==0; musPrimeVec(:)==0; refIndVec(:)==0])   
    error('Zeros remain in optical property vectors, aborting');
end

% #########################################################################
% Define linklist & QM file (remember to delete) ##########################
linklist = DOTHUB_SD2linklist(SD_3Dmesh);
qmfilename = 'tmpfull.qm';
DOTHUB_writeToastQM(qmfilename,SD_3Dmesh.SrcPos,SD_3Dmesh.DetPos,linklist)


% Generate Jacobian #######################################################
% #########################################################################

% Jacobian inputs
c0 = 0.3; % speed of light in vacuum (m/ns ???)
jtype='bicgstab';
bicgstabtol=1e-12;

% Test how long this will take using single channel
% #########################################################################
if testFlag == 1
    disp('Running test calculation...');
    DOTHUB_writeToastQM('tmp.qm',SD_3Dmesh.SrcPos(1,:),SD_3Dmesh.DetPos(1,:),1)
    hMesh.ReadQM('tmp.qm');
    qvec = hMesh.Qvec ('Neumann', 'Gaussian', 2);
    mvec = hMesh.Mvec ('Gaussian', 2,refindvec);
    
    if basisFlag
        hBasis = toastBasis(hMesh,basis,fineBasis);
    else
        hBasis = 0;
    end
    
    if hMesh.isvalid
        fprintf('Running toastJacobianCW test...\n');
        tic;
        Jtest = toastJacobianCW(hMesh, hBasis, qvec, mvec, muaVec(:,1), musPrimeVec(:,1), refindvec, jtype, bicgstabtol); %#ok<NASGU>
        duration = toc;
    end
    %Assume linear with optode number (very approximate)
    fprintf('Test complete...\n');
    fprintf(['Time estimate for full calculation = ' num2str(duration*0.5*nWavs*(size(det_pos,1) + size(source_pos,1))/60) ' mins\n']);
    delete('tmp.qm');
end

% #########################################################################
% Calculate full Jacobian
%Re-set parameters
hMesh.ReadQM(qmfilename)
qvec = hMesh.Qvec ('Neumann', 'Gaussian', 2); 
mvec = hMesh.Mvec ('Gaussian',2,refindvec); 

% Set basis
if basisFlag
    disp('Setting basis...');   
    hBasis = toastBasis(hMesh,basis,fineBasis);
else
    hBasis = 0;
end

for i = 1:nWavs
    fprintf(['Producing Jacobian at wavelength ', num2str(i), '...\n'];
    
    tic
    Jtmp = toastJacobianCW(hMesh, hBasis, qvec, mvec, muaVec(:,i), musPrimeVec(:,i), refindvec(:,i), jtype, bicgstabtol);
    
    c_medium = c0./refindvec(:,i);
    if basisFlag %Multiply by c, then map to volume to GM, delete volume
        Jtmp = Jtmp.*repmat(hBasis.Map('M->S',c_medium),1,size(Jtmp,1))';
        J_voltmp = zeros(size(Jtmp,1),size(headVolumeMesh.node,1)); %Map to volume temporarily
        J_voltmp = hBasis.Map('S->M',Jtmp);
        Jgm(i,:,:) = (vol2gm*J_voltmp')'; 
        J(i,:,:) = Jtmp;
        clear Jtmp J_voltmp
    
    else %In volume nodes
        J(i,:,:) = Jtmp.*c_medium;
        Jgm(i,:,:) = (vol2gm*J_voltmp')'; 
    end
    duration = toc;
    fprintf(['Completed Jacobian at wavelength ', num2str(i), ' in ' num2str(duration/60) ' minutes\n']);
end
delete(qmfilename);

% #########################################################################
%Process complete, write .jac file
% USE CODE SNIPPET FROM DOTHUB_writeJAC.m
ds = datestr(now,'yyyymmDDHHMMSS');
[pathstr, name, ~] = fileparts(rmapFileName);
jacFileName = fullfile(pathstr,[name '_' ds '.jac']);
logData(1,:) = {'Created on: ', ds};
logData(2,:) = {'Derived from rmap file: ', rmapFileName};
logData(3,:) = {'Calculated using: ', transportPackage};
logData(4,:) = {'Optical properties (tissueInd, wavelength,[mua musPrime refInd]): ', opticalPropertiesByTissue};

% Write .jac file ########################################################
[jac, jacFileName] = DOTHUB_writeJAC(jacFileName,logData,J,Jgm,basis);







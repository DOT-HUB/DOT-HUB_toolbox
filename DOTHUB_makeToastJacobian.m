function [jac, jacFileName] = DOTHUB_makeToastJacobian(rmap,basis)

% Calculates Jacobian using Toast++.
%
% #########################################################################
% INPUTS ##################################################################
%
% rmap          : The full path to the rmap file, or the rmap structure itself, containing variables:
%
%                    % SD3Dmesh          :   The SD structure containing registered 3D optode
%                    %                       positions on the mesh, and
%                                            (critically) SD.MeasList
%
%                    % headVolumeMesh     :   The multi-layer volume mesh structure, registered
%                    %                        to the relevant individual. Contains fields:
%                    %                        node, face, elem, labels
%
%                    % gmSurfaceMesh      :   The gm surface mesh structure, registered
%                    %                           to the relevant individual. Contains fields:
%                    %                           node, face.
%
%                    % scalpSurfaceMesh   :   The scalp surface mesh structure, registered
%                    %                           to the relevant individual. Contains fields:
%                    %                           node, face.
%
%                    % vol2gm                :   The sparse matrix mapping from head volume mesh
%                    %                           space to GM surface mesh space
%
% basis        :   (Optiona) 1x3 vector specifying basis dimensions if desired. 
%                  If you don't want to specify a basis don't parse or parse empty
%                  in which case calculation will be in volume unless
%                  nnodes>200k, in which case a basis of [50 50 50] is assigned.
%                  If you parse basis = 0, a basis will not be assigned no
%                  matter how large your volume is.
%            
% OUTPUTS #################################################################
%
% jac                      :  Structure containing all data inputs
%
% jacFileName              :  The full path of the resulting .jac file
%
% jacFileName.jac file     :  File containing all data inputs
%
% #########################################################################
% #########################################################################
% RJC, UCL, May 2019
% EVR, UCL, Modified on June 2019
% RJC, UCL, Modified April 2020 for Github first commit
% #########################################################################

fprintf('############### Running DOTHUB_makeToastJacobian #################\n');

% MANAGE VARIABLES
% #########################################################################

if ischar(rmap)
    rmapFileName = rmap;
    rmap = load(rmapFileName,'-mat');
else
    rmapFileName = rmap.fileName;
end

basisFlag = 0;
nNodeVol = size(rmap.headVolumeMesh.node,1);
nodeLimit = 2e5;%HARD CODE NODE LIMIT AT 200,000
if ~exist('basis','var') %Not parsed. Volume unless nodes exceed limit...
    if nNodeVol>nodeLimit 
        basisFlag = 1;
        basis = [50 50 50];
        fineBasis = basis.*2;
    else
        basis = []; %No basis
    end
elseif isempty(basis)  %Parsed empty. Volume unless nodes exceed limit...
    if nNodeVol>nodeLimit 
        basisFlag = 1;
        basis = [50 50 50];
        fineBasis = basis.*2;
    else
        basis = []; %No basis
    end    
elseif basis==0  %Parsed as zero = calc in volume no matter the size.
    if nNodeVol>nodeLimit 
        warning('Basis parsed as 0: Calculating Jacobian in volume space despite large node count...');
    end
    basis = []; %No basis
elseif ~isempty(basis) %Basis exists and is not empty, use.
    basisFlag = 1;
    fineBasis = basis.*2;
end
   
%testFlag==1 to run initial single-channel calculation test before running full
%calculation. This is useful because the test should only take a few
%minutes, so if it hangs for longer, you can be certain there is a problem,
%whereas the full calculations can take hours, so you can't be sure it is
%working.

% #########################################################################
% Load rmap ###############################################################
[rmapPath, rmapName,~] = fileparts(rmapFileName);
if isempty(rmapPath) %Ensure full path.
    rmapPath = pwd;
    rmapFileName = fullfile(rmapPath,[rmapName '.rmap']);
end
rmap = load(rmapFileName,'-mat'); %head, gm, scalp, SD3Dmesh, vol2gm, logData

% Unpack key fields #######################################################
headVolumeMesh = rmap.headVolumeMesh;
gmSurfaceMesh = rmap.gmSurfaceMesh;
SD3Dmesh = rmap.SD3Dmesh;
wavelengths = SD3Dmesh.Lambda;
vol2gm = rmap.vol2gm;

nElemVol = size(headVolumeMesh.elem,1);
nNodeGM = size(gmSurfaceMesh.node,1);
nChans = size(SD3Dmesh.MeasList,1);
nChansPerWav = sum(SD3Dmesh.MeasList(:,4)==1);
nWavs = length(wavelengths);

% #########################################################################
% MESH PREP ###############################################################
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
node_tmp = headVolumeMesh.node; %leave 4th column in place if exists;
included_nodes = ismember(1:nNodeVol,elem_tmp(:));
errnodes_ind = find(~included_nodes);
if ~isempty(errnodes_ind)
    warning('Erroneous nodes found, attempting to correct mesh...\n');
    %correct node list
    node_corr = node_tmp(included_nodes,:);
    %correct element list
    for wav = 1:length(errnodes_ind)
        elem_tmp(elem_tmp > (errnodes_ind(wav)-(wav-1))) = elem_tmp(elem_tmp > (errnodes_ind(wav)-(wav-1)))-1;
    end
    
    included_nodes = ismember(1:size(node_corr,1),elem_tmp(:));
    errnodes_ind = find(~included_nodes,1);
    if isempty(errnodes_ind)
        fprintf('Erroneous nodes removed...\n');
        headVolumeMesh.node = node_corr;
        headVolumeMesh.elem(:,1:4) = elem_tmp;
        nNodeVol = size(headVolumeMesh.node,1);
        clear elem_tmp node_tmp node_corr included_nodes errnodes_ind
        rewriteRMAP = 1;
    else
        error('Correction failed: erroneous nodes remain in rmap.headVolumeMesh');
    end
end

% Update nNodeVol

% Clear MEX
clear MEX

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
    DOTHUB_writeRMAP(rmapFileName,rmap.logData,rmap.SD3Dmesh,headVolumeMesh,rmap.gmSurfaceMesh,rmap.scalpSurfaceMesh,rmap.vol2gm);
end

% #########################################################################
% Set optical properties
% First assign
muaVec = zeros(nWavs,nNodeVol);
musPrimeVec = zeros(nWavs,nNodeVol);
refIndVec = zeros(nWavs,nNodeVol);

% Determine tissue optical properties and populate vectors
nTissues = length(headVolumeMesh.labels);
for tiss = 1:nTissues
    tmpInd = find(strcmpi({'scalp','skull','ECT','CSF','GM','WM'},headVolumeMesh.labels{tiss}), 1);
    if isempty(tmpInd)
        error('Unknown tissue label, please correct rmap.headVolumeMesh.labels and try again');
    end
    tissueNodeList = headVolumeMesh.node(:,4)==tiss;
    for wav = 1:nWavs
        [mua, musPrime, refInd] = DOTHUB_getTissueCoeffs(headVolumeMesh.labels{tiss},wavelengths(wav));%Potentially update this to be age-specific?
        muaVec(wav,tissueNodeList) = mua;
        musPrimeVec(wav,tissueNodeList) = musPrime;
        refIndVec(wav,tissueNodeList) = refInd;
        
        opticalPropertiesByTissue(tiss,wav,1) = mua;
        opticalPropertiesByTissue(tiss,wav,2) = musPrime;
        opticalPropertiesByTissue(tiss,wav,3) = refInd;
    end
end

if any([muaVec(:)==0; musPrimeVec(:)==0; refIndVec(:)==0])   
    error('Zeros remain in optical property vectors, aborting');
end

% #########################################################################
% Define linklist & QM file (remember to delete) ##########################
linklist = DOTHUB_SD2linklist(SD3Dmesh);
qmfilename = 'tmpfull.qm';
DOTHUB_writeToastQM(qmfilename,SD3Dmesh.SrcPos,SD3Dmesh.DetPos,linklist)

% Generate Jacobian #######################################################
% #########################################################################

% Jacobian inputs
c0 = 0.3; % speed of light in vacuum (m/ns ???)
jtype='bicgstab';
bicgstabtol=1e-12;

% Test how long this will take using single channel
% #########################################################################
% if testFlag == 1
%     try
%     disp('Running test calculation...');
%     DOTHUB_writeToastQM('tmp.qm',SD3Dmesh.SrcPos(1,:),SD3Dmesh.DetPos(1,:),1)
%     hMesh.ReadQM('tmp.qm');
%     qvec = hMesh.Qvec ('Neumann', 'Gaussian', 2);
%     mvec = hMesh.Mvec ('Gaussian', 2, refIndVec(1,:)');
%     
%     if basisFlag
%         hBasis = toastBasis(hMesh,basis,fineBasis);
%     else
%         hBasis = 0;
%     end
%     
%     if hMesh.isvalid
%         fprintf('Running toastJacobianCW test...\n');
%         tic;
%         Jtest = toastJacobianCW(hMesh, hBasis, qvec, mvec, muaVec(1,:)', musPrimeVec(1,:)', refIndVec(1,:)', jtype, bicgstabtol); %#ok<NASGU>
%         duration = toc;
%     end
%     Assume linear with optode number (very approximate)
%     fprintf('Test complete...\n');
%     fprintf(['Estimated processing time estimate for full Jacobian, per wavelength = ' num2str(duration*0.5*(SD3Dmesh.nDets + SD3Dmesh.nSrcs)/60,'%0.2f') ' mins\n']);
%     delete('tmp.qm');
%     catch
%         fprintf('Test of toastJacobianCW failed. Please check inputs...\n');
%     end
% end

% #########################################################################
% Calculate full Jacobian
%Re-set parameters
hMesh.ReadQM(qmfilename)
qvec = hMesh.Qvec ('Neumann', 'Gaussian', 2);
mvec = hMesh.Mvec ('Gaussian', 2, refIndVec(1,:)');

% Set basis
if basisFlag
    disp('Setting basis...');   
    hBasis = toastBasis(hMesh,basis,fineBasis);
    nNodeNat = hBasis.slen;
else
    hBasis = 0;
    nNodeNat = nNodeVol;
end

%pre-allocate
J = cell(nWavs,1);
c_medium = (c0./refIndVec(1,:));
for wav = 1:nWavs
    fprintf(['Producing Jacobian at wavelength ', num2str(wav), '...\n']);
    
    tic
    Jtmp = toastJacobianCW(hMesh, hBasis, qvec, mvec, muaVec(wav,:)', musPrimeVec(wav,:)', refIndVec(wav,:)', jtype, bicgstabtol);
    
    if basisFlag %Multiply by c, then map to volume to GM, delete volume
        J{wav}.basis = Jtmp.*repmat(hBasis.Map('M->S',c_medium),1,size(Jtmp,1))';
        for chan = 1:nChansPerWav
            J{wav}.vol(chan,:) = hBasis.Map('S->M',J{wav}.basis(chan,:)');
        end
        J{wav}.gm = (vol2gm*J{wav}.vol')'; 
        
        %Clear things, empty J.vol as we are in basis
        clear Jtmp
        J{wav}.vol = [];
        
    else %In volume nodes
        J{wav}.vol = Jtmp.*repmat(c_medium,nChansPerWav,1);
        J{wav}.gm = (vol2gm*J{wav}.vol')'; 
        J{wav}.basis = [];
    end
    duration = toc;
    fprintf(['Completed Jacobian at wavelength ', num2str(wav), ' in ' num2str(duration/60) ' minutes\n']);
end
delete(qmfilename);

% #########################################################################
%Process complete, write .jac file
% USE CODE SNIPPET FROM DOTHUB_writeJAC.m
ds = datestr(now,'yyyymmDDHHMMSS');
[pathstr, name, ~] = fileparts(rmapFileName);
jacFileName = fullfile(pathstr,[name '.jac']);
transportPackage = 'toast';
logData(1,:) = {'Created on: ', ds};
logData(2,:) = {'Derived from rmap file: ', rmapFileName};
logData(3,:) = {'Calculated using: ', transportPackage};
logData(4,:) = {'Optical properties (tissueInd, wavelength,[mua musPrime refInd]): ', opticalPropertiesByTissue};

% Write .jac file ########################################################
[jac, jacFileName] = DOTHUB_writeJAC(jacFileName,logData,J,basis);







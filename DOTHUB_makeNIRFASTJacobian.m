function [jac, jacFileName] = DOTHUB_makeNIRFASTJacobian(rmap)

% Calculates Jacobian using TNIRFAST
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
% RJC, UCL, April 2020 for Github first commit
% #########################################################################

fprintf('############## Running DOTHUB_makeNIRFASTJacobian ################\n');

% MANAGE VARIABLES
% #########################################################################

if ischar(rmap)
    rmapFileName = rmap;
    rmap = load(rmapFileName,'-mat');
else
    rmapFileName = rmap.fileName;
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
NIRFASTmesh.nodes = headVolumeMesh.node(:,1:3);
NIRFASTmesh.elem = headVolumeMesh.elem(:,1:4);
NIRFASTmesh.link = SD3Dmesh.MeasList(SD3Dmesh.MeasList(:,4)==1,[1 2 4]);

% #########################################################################
% Set optical properties
% First assign
muaVec = zeros(nWavs,nNodeVol);
musPrimeVec = zeros(nWavs,nNodeVol);
refIndVec = zeros(nWavs,nNodeVol);

% Determine tissue optical properties and populate vectors
nTissues = length(headVolumeMesh.labels);
for wav = 1:nTissues
    tmpInd = find(strcmpi({'scalp','skull','ECT','CSF','GM','WM'},headVolumeMesh.labels{wav}), 1);
    if isempty(tmpInd)
        error('Unknown tissue label, please correct rmap.headVolumeMesh.labels abd try again');
    end
    tissueNodeList = headVolumeMesh.node(:,4)==wav;
    for j = 1:nWavs
        [mua, musPrime, refInd] = DOTHUB_getTissueCoeffs(headVolumeMesh.labels{wav},wavelengths(j));%Potentially update this to be age-specific?
        muaVec(j,tissueNodeList) = mua;
        musPrimeVec(j,tissueNodeList) = musPrime;
        refIndVec(j,tissueNodeList) = refInd;
        
        opticalPropertiesByTissue(wav,j,1) = mua;
        opticalPropertiesByTissue(wav,j,2) = musPrime;
        opticalPropertiesByTissue(wav,j,3) = refInd;
    end
end

if any([muaVec(:)==0; musPrimeVec(:)==0; refIndVec(:)==0])   
    error('Zeros remain in optical property vectors, aborting');
end


% Generate Jacobian #######################################################
% #########################################################################

% Jacobian inputs
c0 = 0.3; % speed of light in vacuum (m/ns ???)
jtype='bicgstab';
bicgstabtol=1e-12;

% Test how long this will take using single channel
% #########################################################################
if testFlag == 1
    try
    disp('Running test calculation...');
    DOTHUB_writeToastQM('tmp.qm',SD3Dmesh.SrcPos(1,:),SD3Dmesh.DetPos(1,:),1)
    hMesh.ReadQM('tmp.qm');
    qvec = hMesh.Qvec ('Neumann', 'Gaussian', 2);
    mvec = hMesh.Mvec ('Gaussian', 2, refIndVec(1,:)');
    
    if basisFlag
        hBasis = toastBasis(hMesh,basis,fineBasis);
    else
        hBasis = 0;
    end
    
    if hMesh.isvalid
        fprintf('Running toastJacobianCW test...\n');
        tic;
        Jtest = toastJacobianCW(hMesh, hBasis, qvec, mvec, muaVec(1,:)', musPrimeVec(1,:)', refIndVec(1,:)', jtype, bicgstabtol); %#ok<NASGU>
        duration = toc;
    end
    %Assume linear with optode number (very approximate)
    fprintf('Test complete...\n');
    fprintf(['Estimated processing time estimate for full Jacobian, per wavelength = ' num2str(duration*0.5*(SD3Dmesh.nDets + SD3Dmesh.nSrcs)/60,'%0.2f') ' mins\n']);
    delete('tmp.qm');
    catch
        fprintf('Test of toastJacobianCW failed. Please check inputs...\n');
    end
end

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







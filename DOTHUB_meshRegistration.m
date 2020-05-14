function [rmap, rmapFileName] = DOTHUB_meshRegistration(SD3DFileName,origMeshFileName)

% This function takes the landmarks and source and detector locations for an individual
% and registers the selected mesh (atlas or subject specific) into the native space of the
% individual. It uses the landmarks to transform between the two spaces.
% We hope to expand registration options, but at present this is performed
% with a simple Affine transformation.
%
%The result is saved to the same location as the SD3DFileName
%
% ################################ INPUTS #################################
%
% SD3DFullFileName     :   The full path of the 3D .SD file (or .nirs file), which should
%                           itself be named after the posData.csv from which is
%                           was derived
%
% origMeshFullFileName  :   The full path to the original mesh .mshs file, which contains 
%                           the following variables:
%    
%                             % headVolumeMesh    :   The multi-layer volume mesh structure, Contains
%                             %                       fields node, face, elem, labels
%                                                   
%                             % gmSurfaceMesh     :   The gm surface mesh structure, Contains fields:
%                             %                       node, face.
%
%                             % scalpSurfaceMesh  :   The scalp surface mesh structure, contains fields:
%                             %                       node, face.
%
%                             % landmarks         :   5 x 3 matrix of the cranial landmark coordinates on
%                             %                       the mesh surface. Nz, Iz, Ar, Al, Cz.
%
%                             % vol2gm            :   (Optional) The sparse matrix mapping from head volume mesh
%                             %                       space to GM surface mesh space. Calculated if not parsed
%
% ####################### OUTPUTS #########################################
%
% rmapFileName      :   The path of the .rmap file that results from this
%                       registration.
%
% ####################### Dependencies ####################################
%
% #########################################################################
% UCL, EVR, 5th June 2019 & RJC, April 2020
%
% ############################# To Do #####################################
% Convert inputs to a single .mshs file, which would be the output of a
% mask-to-meshes function?
% #########################################################################
%
% ############################# Updates ###################################
% #########################################################################

fprintf('################ Running DOTHUB_meshRegistration #################\n');

% MANAGE VARIABLES
% #########################################################################

%Enforce full path names here
% Load SD3D
load(SD3DFileName,'SD3D','-mat');
% Load SD3D
load(origMeshFileName,'-mat');

if ~exist('vol2gm','var')
    radius = 3;
    vol2gm = DOTHUB_vol2gmMap(headVolumeMesh.node,gmSurfaceMesh.node,radius,1);
elseif isempty(vol2gm)
    radius = 3;
    vol2gm = DOTHUB_vol2gmMap(headVolumeMesh.node,gmSurfaceMesh.node,radius,1);
end



%% AFFINE METHOD (may add other options in future)
regMethod = 'Affine';
% Get affine transformation matrices
[A,B] = DOTHUB_affineMap(landmarks, SD3D.Landmarks);

% Transform meshnodes;
headVolumeMeshReg = headVolumeMesh;
headVolumeMeshReg.node(:,1:3) = DOTHUB_affineTrans(headVolumeMeshReg.node(:,1:3),A,B);

gmSurfaceMeshReg = gmSurfaceMesh;
gmSurfaceMeshReg.node(:,1:3) = DOTHUB_affineTrans(gmSurfaceMeshReg.node(:,1:3),A,B);

scalpSurfaceMeshReg = scalpSurfaceMesh;
scalpSurfaceMeshReg.node(:,1:3) = DOTHUB_affineTrans(scalpSurfaceMeshReg.node(:,1:3),A,B);

%Ensure sources are on the mesh - first force to scalp, then to volume.
SD3Dmesh = SD3D;
for i = 1:SD3D.nSrcs
    tmp = DOTHUB_nearestNode(SD3D.SrcPos(i,:),scalpSurfaceMeshReg.node(:,1:3));
    SD3Dmesh.SrcPos(i,:) = DOTHUB_nearestNode(tmp,headVolumeMeshReg.node(:,1:3));
end
for i = 1:SD3D.nDets
    tmp = DOTHUB_nearestNode(SD3D.DetPos(i,:),scalpSurfaceMeshReg.node(:,1:3));
    SD3Dmesh.DetPos(i,:) = DOTHUB_nearestNode(tmp,headVolumeMeshReg.node(:,1:3));
end

%% Write out .rmap file ####################################################
% USE CODE SNIPPET FROM DOTHUB_writeRMAP to define name and logData
ds = datestr(now,'yyyymmDDHHMMSS');
[~, origMeshFileName, ~] = fileparts(origMeshFileName);
[SD3DPath, SD3DFileName, ~] = fileparts(SD3DFileName);
rmapFileName = fullfile(SD3DPath,[SD3DFileName '.rmap']);
logData(1,:) = {'Created on: ',ds};
logData(2,:) = {'Positions derived from: ', SD3DFileName};
logData(3,:) = {'Meshes derived from: ', origMeshFileName};
logData(4,:) = {'Registration method: ', regMethod};

% Write
[rmap, rmapFileName] = DOTHUB_writeRMAP(rmapFileName,logData,SD3Dmesh,headVolumeMeshReg,gmSurfaceMeshReg,scalpSurfaceMeshReg,vol2gm);


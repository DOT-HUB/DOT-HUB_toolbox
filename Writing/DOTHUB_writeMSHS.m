function [rmap, rmapFileName] = DOTHUB_writeRMAP(rmapFileName,logData,SD3Dmesh,headVolumeMesh,gmSurfaceMesh,scalpSurfaceMesh,vol2gm)

% This script creates a registered mesh and positions (.rmap) file and associated struct 

% ####################### INPUTS ##########################################

% rmapFileName      :  The desired path &/ filename for the .rmap file.
%                      This can be anything, but we recommend this variable be defined with the
%                      following code snippet, where: origMeshFullFileName = full path and name
%                      of mesh being registered; SD3DFullFileName = full path and name of 
%                      SD3D file that is the basis of the registration; and regMethod = 'Affine' 
%                      or equivalent. This snippet also provides recommended input variable 'logData'. 

                       %ds = datestr(now,'yyyymmDDHHMMSS');
                       %[~, origMeshFileName, ~] = fileparts(origMeshFullFileName);
                       %[SD3DPath, SD3DFileName, ~] = fileparts(SD3DFullFileName);
                       %rmapFileName = fullfile(SD3DPath,[origMeshFileName '_Reg2_' SD3DFileName '.rmap']);
                       %logData(1,:) = {'Created on: ',ds};
                       %logData(2,:) = {'Positions derived from: ', SD3DFullFileName};
                       %logData(3,:) = {'Meshes derived from: ', origMeshFullFileName};
                       %logData(4,:) = {'Registration method: ', regMethod};

% logData           :  (Optional). logData is a cell array of strings containing useful
%                      info as per snippet above. Parse empty to ignore.

% SD3Dmesh          :  The SD structure containing registered 3D optode
%                      positions
%
% headVolumeMesh    :  The multi-layer volume mesh structure, registered
%                      to the relevant individual. Contains fields:
%                      node, face, elem, labels
%
% gmSurfaceMesh     :  The gm surface mesh structure, registered
%                      to the relevant individual. Contains fields:
%                      node, face.
%
% scalpSurfaceMesh  :  The scalp surface mesh structure, registered
%                      to the relevant individual. Contains fields:
%                      node, face.
%
% vol2gm            :  The sparse matrix mapping from head volume mesh
%                      space to GM surface mesh space
%
%
% ####################### OUTPUTS #########################################
%
% rmap              : A structure containing all fields for rmap
%
% rmapFileName      : The full path of the resulting .rmap file
%
% .rmap             :  A file containing a structure of:
%                      rmap.logData             - as defined above
%                      rmap.SD3Dmesh            - as defined above
%                      rmap.headVolumeMesh      - as defined above
%                      rmap.gmSurfaceMesh       - as defined above
%                      rmap.scalpSurfaceMesh    - as defined above
%                      rmap.vol2gm              - as defined above
%                      rmap.fileName            - the path of the saved rmap file


% ####################### Dependencies ####################################
% #########################################################################
% RJC, UCL, April 2020
%
% ############################# Updates ###################################
% #########################################################################

% MANAGE VARIABLES
% #########################################################################
% Check all optodes positions are on-mesh to allow save
for i = 1:SD3Dmesh.nSrcs
    tmp = sum(headVolumeMesh.node(:,1) == SD3Dmesh.SrcPos(i,1) & headVolumeMesh.node(:,2) == SD3Dmesh.SrcPos(i,2) & headVolumeMesh.node(:,3) == SD3Dmesh.SrcPos(i,3));
    if tmp==0
        error('Optodes found that are not on nodes of headVolumeMesh, please correct before saving rmap');
    end
end
for i = 1:SD3Dmesh.nDets
    tmp = sum(headVolumeMesh.node(:,1) == SD3Dmesh.DetPos(i,1) & headVolumeMesh.node(:,2) == SD3Dmesh.DetPos(i,2) & headVolumeMesh.node(:,3) == SD3Dmesh.DetPos(i,3));
    if tmp==0
        error('Optodes found that are not on nodes of headVolumeMesh, please correct before saving rmap');
    end
end

if isempty(logData)
    logData = {};
    warning('logData is empty: this might make it harder to keep track of your data...');
end

%Create rmap struct #######################################################
rmap = struct('headVolumeMesh',headVolumeMesh,'gmSurfaceMesh',gmSurfaceMesh,'scalpSurfaceMesh',scalpSurfaceMesh,...
    'vol2gm',vol2gm,'SD3Dmesh',SD3Dmesh);
rmap.logData = logData;

%Create filename ##########################################################
[pathstr, name, ext] = fileparts(rmapFileName);
if isempty(ext) || ~strcmpi(ext,'.rmap')
    ext = '.rmap';
end
if isempty(pathstr)
    pathstr = pwd;
end

rmapFileName = fullfile(pathstr,[name ext]);
rmap.fileName = rmapFileName; %including the fileName within the structure is very useful 
%for tracking and naming things derived further downstream.

if exist(rmapFileName,'file')
    warning([name ext ' will be overwritten...']);
end

%Save .rmap file ###########################################################
save(rmapFileName,'-struct','rmap');
fprintf('###################### Writing .rmap file ########################\n');
fprintf(['.rmap data file saved as ' rmapFileName '\n']);
fprintf('\n');

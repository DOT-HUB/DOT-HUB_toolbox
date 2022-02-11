function [mshs, mshsFileName] = DOTHUB_writeMSHS(mshsFileName,logData,headVolumeMesh,gmSurfaceMesh,scalpSurfaceMesh,vol2gm,landmarks,tenFive)

% This script creates a mesh file created from a given head model. Could be from an atlas or from an individual 
%
% ####################### INPUTS ##########################################
%
% mshsFileName      :  The desired path &/ filename for the .mshs file.
%                      This can be anything, but we recommend this variable be defined with the
%                      following code snippet, where: origMeshFullFileName = full path and name
%                      of mesh being registered. This snippet also provides recommended input variable 'logData'. 

                       %ds = datestr(now,'yyyymmDDHHMMSS');
                       % @Rob not sure if there are other things you want
                       % added in here
                       %logData(1,:) = {'Created on: ',ds};

% logData           :  (Optional). logData is a cell array of strings containing useful
%                      info as per snippet above. Parse empty to ignore.
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
% landmarks         :  (Optional) A matrix of cranial landmark positions on this mesh.
%                      (5x3 Nz, Iz, Ar, Al, Cz)
%
% tenFive           :  (Optional) The 10-5 locations for this mesh: .positions (n x 3)
%                      'labels {n}
%
% ####################### OUTPUTS #########################################
%
% mshs              : A structure containing all fields for mshs
%
% mshsFileName      : The full path of the resulting .mshs file
%
% .mshs             :  A file containing a structure of:
%                      mshs.logData             - as defined above
%                      mshs.headVolumeMesh      - as defined above
%                      mshs.gmSurfaceMesh       - as defined above
%                      mshs.scalpSurfaceMesh    - as defined above
%                      mshs.vol2gm              - as defined above
%                      mshs.Landmarks           - as defined above
%                      mshs.tenFive             - as defined above
%                      mshs.fileName            - the path of the saved mshs file


% ####################### Dependencies ####################################
% #########################################################################
% RJC, UCL, April 2020
%
% ############################# Updates ###################################
% #########################################################################

% MANAGE VARIABLES
% #########################################################################


if isempty(logData)
    logData = {};
    warning('logData is empty: this might make it harder to keep track of your data...');
end
if ~exist('logData','var')
    mshs.logData = [];
else
    mshs.logData = logData;
end

%Create mshs struct #######################################################
mshs = struct('headVolumeMesh',headVolumeMesh,'gmSurfaceMesh',gmSurfaceMesh,'scalpSurfaceMesh',scalpSurfaceMesh,'vol2gm',vol2gm);
if ~exist('tenFive','var')
    mshs.tenFive = [];
else
    mshs.tenFive = tenFive;
end
if ~exist('landmarks','var')
    mshs.landmarks = [];
else
    mshs.landmarks = landmarks;
end

mshs.logData = logData;


%Create filename ##########################################################
[pathstr, name, ext] = fileparts(mshsFileName);
if isempty(ext) || ~strcmpi(ext,'.mshs')
    ext = '.mshs';
end
if isempty(pathstr)
    pathstr = pwd;
end

mshsFileName = fullfile(pathstr,[name ext]);
mshs.fileName = mshsFileName; %including the fileName within the structure is very useful 
%for tracking and naming things derived further downstream.

if exist(mshsFileName,'file')
    warning([name ext ' will be overwritten...']);
end

%Save .mshs file ###########################################################
save(mshsFileName,'-struct','mshs');
fprintf('###################### Writing .mshs file ########################\n');
fprintf(['.mshs data file saved as ' mshsFileName '\n']);
fprintf('\n');
end

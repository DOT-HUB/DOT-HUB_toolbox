function [allMeshes,allMeshesFileName] = DOTHUB_makeMesh(allMeshesFileName,logData,mask,tissueLabels,landmarks,opt)

% Creates volumetric and surface meshes from tissue mask and creates the associated struct.
%
% #########################################################################
% INPUTS ##################################################################
%
% allMeshesFileName       :  The desired path &/ filename for the .mshs file.
%                            This can be anything, but we recommend this variable be defined with the
%                            following code snippet, where: origMaskFullFileName = full path and name
%                            of mask being used; landmarkFullFileName = full path and name of 
%                            landmarks file; This snippet also provides recommended input variable 'logData'. 
%        
                             %ds = datestr(now,'yyyymmDDHHMMSS');
                             %[maskPath, origMaskFileName, ~] = fileparts(origMaskFullFileName);
                             %[~, landmarkFileName, ~] = fileparts(landmarkFullFileName);
                             %allMeshesFileName = fullfile(maskPath,[origMaskFileName '_landmark_' landmarkFileName '.mshs']);
                             %logData(1,:) = {'Created on: ',ds};
                             %logData(2,:) = {'Landmarks derived from: ', landmarkFullFileName};
                             %logData(3,:) = {'Mask derived from: ', origMaskFullFileName};
%
% logData                 :  (Optional). logData is a cell array of strings containing useful
%                            info as per snippet above. Parse empty to ignore.
%
% mask                    :  The full path to the tissue mask file or the tissue
%                            mask. The mask should be saved in mask.img while in
%                            mask.voxelSize we should have the size of the voxel
%
% tissueLabels            :  The labels of the tissues contained in the tissue mask.
%                            It should be organized as a cell array. The allowed entries are: 
%                            scalp, skull, ect, csf, gm, wm and air. There are two
%                            compulsory entries: scalp (or ect) and gm.
%
% landmarks               :  A matrix with cranial landmark positions on
%                            the mesh (5x3: Nz, Iz, Ar, Al, Cz)
%
% opt                     :  struct containing options to build the mesh.
%                            Possible fields are: 
%                            
%                           - maxnodeVolMesh/maxnodeScalpMesh/maxnodeGmMesh:
%                             max number of nodes for each mesh type
%                             (default 800000 for headVolumeMesh and 50000
%                             for surface meshes)
%
%                           - radboundVolMesh/radboundScalpMesh/radboundGmMesh:
%                             max radius of the Delaunay sphere (defaul = 1)
%
%                           - side: 'upper' or 'lower', defines if threshold is at 
%                             upper or lower interface (default: lower)
%
%                           - maxVol: target maximum tetrahedral elem
%                             volume. It can be either a scalar or a different value    
%                             for each tissue type (array of values) (default = 1) 
%
%                           - downSamplingScalpMesh/downSamplingGmMesh:
%                             downsampling factor for the low density
%                             surface meshes (default = 0.3)
%
%
% OUTPUTS #################################################################
%
% allMeshes               :  A structure containing all fields for allMeshes
%
% allMeshesFileName       :  The full path of the resulting .mshs file
%
% .mshs                   :  A file containing a structure of:
%
%                           % headVolumeMesh     :   The multi-layer volume mesh structure. 
%                                                    Contains fields: node, face, elem, labels
%
%                           % gmSurfaceMesh      :   The gm surface mesh structure. 
%                                                    Contains fields: node, face.
%
%                           % scalpSurfaceMesh   :   The scalp surface mesh structure.
%                                                    Contains fields: node, face.
%
%                           % vol2gm             :   The sparse matrix mapping from head volume mesh
%                                                    space to GM surface mesh space
%
%                           % landmarks          :   A matrix containing the
%                                                    landmarks coordinate
%
%                           % tenFive            :   ten-five locations for
%                                                    the mesh (.positions
%                                                    (nx3) and .labels {n})
%
%                           % logData            :   As defined above
%
%                           % fileName           :   The path of the saved mshs file
%
%
%
% ####################### Dependencies ####################################
% #########################################################################
%
% iso2mesh, interparc
%
% ############################# Updates ###################################
% #########################################################################
%
% Sabrina Brigadoi, 21/05/2020
%
% #########################################################################

%%%%% Checking all inputs are correct %%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    error('Not enough input arguments, please check the function help');
end

% Check that scalp/ect and gm labels are present in the mask
nTissueScalp = [];
nTissueGm = [];
nTissueWm = [];
for iT = 1:length(tissueLabels)
    if strcmpi(tissueLabels{iT},'scalp')
        nTissueScalp = iT;
    elseif strcmpi(tissueLabels{iT},'ect')
        nTissueScalp = iT;
    elseif strcmpi(tissueLabels{iT},'gm')
        nTissueGm = iT;
    elseif strcmpi(tissueLabels{iT},'wm')
        nTissueWm = iT;
    end    
end

if isempty(nTissueScalp) 
    error('No scalp or ECT tissue found in the mask, please correct before creating the mesh');
end

if isempty(nTissueGm) 
    error('No GM tissue found in the mask, please correct before creating the mesh');
end

if isempty(logData)
    logData = {};
    warning('logData is empty: this might make it harder to keep track of your data...');
end

% Check that 5 landmarks are provided
if size(landmarks,1) ~= 5
    error('5 landmarks should be provided, please correct before creating the mesh');
end

if nargin < 6 % no opt
    opt = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% Create volumetric mesh %%%%%%%%%%%%%%%%%%%%%%%%%%

% Create allMeshes struct with given fields
allMeshes.logData = logData;

% If the mask is passed as fileName, load it
if ischar(mask)
    load(mask);
end

% Create volumetric mesh 
if ~isfield(opt,'maxnodeVolMesh')
    optVolMesh.maxnode = 800000;
else
    optVolMesh.maxnode = opt.maxnodeVolMesh;
end
if ~isfield(opt,'radboundVolMesh')
    optVolMesh.radbound = 1;
else
    optVolMesh.radbound = opt.radboundVolMesh;
end
if ~isfield(opt,'side')
    optVolMesh.side = 'lower';
else
    optVolMesh.side = opt.side;
end

% max volume of the tetrahedron
if ~isfield(opt,'maxVol')
    maxVol = 1;
else
    if length(opt.maxVol) == 1
        maxVol = opt.maxVol;
    else % one size for each tissue label
        maxVol = [];
        for iT = 1:length(opt.maxVol)
            maxVol = [maxVol num2str(iT) '=' num2str(opt.maxVol(iT)) ':'];
        end   
        maxVol = maxVol(1:end-1);
    end
end

% Size of the mask
dim = size(mask.img);

% Create tetrahedral mesh
[node,elem,face] = vol2mesh(uint8(mask.img),1:dim(1),1:dim(2),1:dim(3),optVolMesh,maxVol,1,'cgalmesh');

% Reorient the element of the mesh
elem_or = meshreorient(node(:,1:3),elem(:,1:4)) ;
elem = [elem_or,elem(:,5)];

% The output of iso2mesh is in grid/voxel_size: in order to make the mesh having the same dimension
% of the real MRI we need to multiply for the real voxel size
if length(mask.voxelSize) > 1
    mask.voxelSize = repmat(mask.voxelSize,size(node,1),1);
end

node(:,1:3) = node(:,1:3).*mask.voxelSize;     

% The 4th column of node should contain the label of the tissue type
tiss = unique(elem(:,5));

% Use histograms to compute how frequently a node is shared by elements in
% the same tissue and then assign to that node the tissue having the
% highest frequency
tmp = zeros(size(node,1),length(tiss));
for i_t = 1:length(tiss)
    tiss_nodes = elem(elem(:,5) == i_t,1:4);
    [n,~] = hist(tiss_nodes(:),1:size(node,1));
    tmp(unique(tiss_nodes(:)),i_t) = n(unique(tiss_nodes(:)));
end
[~,ind] = max(tmp,[],2);
node(:,4) = ind;

% Save all the volumetric mesh output in allMeshes                                           
allMeshes.headVolumeMesh.node = node;
allMeshes.headVolumeMesh.elem = elem;
allMeshes.headVolumeMesh.face = face;
allMeshes.headVolumeMesh.labels = tissueLabels;

% Scale the landmarks coordinates by voxel size and save in the struct
allMeshes.landmarks = landmarks.*mask.voxelSize;

%%%%%%%%%%%%%% Create scalp surface mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select the scalp tissue from the mask
mask_tissue = zeros(dim);
mask_tissue(mask.img==nTissueScalp) = 1;

% Create low density surface mesh
if ~isfield(opt,'maxnodeScalpMesh')
    optScalpMesh.maxnode = 500000;
else
    optScalpMesh.maxnode = opt.maxnodeScalpMesh;
end
if ~isfield(opt,'radboundScalpMesh')
    optScalpMesh.radbound = 1;
else
    optScalpMesh.radbound = opt.radboundScalpMesh;
end
if ~isfield(opt,'downSamplingScalpMesh')
    optScalpMesh.downSampling = 0.3;
else
    optScalpMesh.downSampling = opt.downSamplingScalpMesh;
end

[nodeS,faceS] = createSurfMesh(mask_tissue,optScalpMesh);

% Save all the scalp surface mesh output in allMeshes   
allMeshes.scalpSurfaceMesh.node = nodeS.*mask.voxelSize;
allMeshes.scalpSurfaceMesh.face = faceS;

%%%%%%%%%%%%%% Create GM surface mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select the GM tissue from the mask
if exist('nTissueWm','var')
    nTissueGm = [nTissueGm nTissueWm];
end
mask_tissue = zeros(dim);
for i_n = 1:length(nTissueGm)
    mask_tissue(mask.img==nTissueGm(i_n)) = 1;
end

% Create low density surface mesh
if ~isfield(opt,'maxnodeGmMesh')
    optGmMesh.maxnode = 500000;
else
    optGmMesh.maxnode = opt.maxnodeGmMesh;
end
if ~isfield(opt,'radboundGmMesh')
    optGmMesh.radbound = 1;
else
    optGmMesh.radbound = opt.radboundGmMesh;
end
if ~isfield(opt,'downSamplingGmMesh')
    optGmMesh.downSampling = 0.3;
else
    optGmMesh.downSampling = opt.downSamplingGmMesh;
end

[nodeGM,faceGM] = createSurfMesh(mask_tissue,optGmMesh);

% Save all the scalp surface mesh output in allMeshes   
allMeshes.gmSurfaceMesh.node = nodeGM.*mask.voxelSize;
allMeshes.gmSurfaceMesh.face = faceGM;

%%%%%%%%%%%%%%% Create vol2gm matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
radius = 3;
allMeshes.vol2gm = Vol2GM_Transform(allMeshes.headVolumeMesh.node,allMeshes.gmSurfaceMesh.node,radius,0);

%%%%%%%%%%%%%%% Get 10-5 positions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Esctract the node of the external surface of the scalp
face_idx_ext_surf = faceneighbors(allMeshes.headVolumeMesh.elem(:,1:4),'surface');
node_idx_ext_surf = unique(face_idx_ext_surf(:));
surf = allMeshes.headVolumeMesh.node(node_idx_ext_surf,1:3); % creates the surface
    
% Bring landmarks' coordinates to the surface mesh
pts = bringPtsToSurf(surf,allMeshes.landmarks);

% Compute the 10-5 position
[refpts_10_5,refpts_10_5_label] = DOTHUB_getTen5points(surf,refpts);

allMeshes.tenFive.positions = refpts_10_5;
allMeshes.tenFive.labels = refpts_10_5_label;
    
%%%%%%%%%%%%%%% Create filename %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[pathstr,name,ext] = fileparts(allMeshesFileName);
if isempty(ext) || ~strcmpi(ext,'.mshs')
    ext = '.mshs';
end
if isempty(pathstr)
    pathstr = pwd;
end

allMeshesFileName = fullfile(pathstr,[name ext]);
allMeshes.fileName = allMeshesFileName; % including the fileName within the structure is very useful 
%for tracking and naming things derived further downstream.

if exist(allMeshesFileName,'file')
    warning([name ext ' will be overwritten...']);
end

%%%%%%%%%%%%%% Save .mshs file $%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(allMeshesFileName,'-struct','allMeshes');
fprintf('###################### Writing .mshs file ########################\n');
fprintf(['.mshs data file saved as ' allMeshesFileName '\n']);
fprintf('\n');

end

%%%%%%% functions %%%%%%%%%%%%%%%%%%%%%%

function [node_s,face_r] = createSurfMesh(mask_tissue,opt)

dim = size(mask_tissue);

% Fill the obtained tissue mask; in this way the surface mesh
% created will have only the external surface of the tissue and not
% the internal boundary one
mask_filled = zeros(dim);
for i_slice = 1:dim(3)
    mask_filled(:,:,i_slice) = imfill(mask_tissue(:,:,i_slice),'holes');
end

% Create surface mesh
[node_surf,face_surf] = vol2surf(mask_filled,1:dim(1),1:dim(2),1:dim(3),opt,1);

% Repairing the mesh, downsampling it and smoothing it
[node_r,face_r] = meshcheckrepair(node_surf,face_surf);
[newnode,newface] = meshresample(node_r,face_r,opt.downSampling); % downsample
[node_r,face_r] = meshcheckrepair(newnode,newface);
conn = meshconn(face_r,size(node_r,1));
node_s = smoothsurf(node_r,[],conn,10,0.7,'lowpass'); % smoothing
[node_s,face_r]=removeisolatednode(node_s,face_r);
[node_s,face_r] = meshcheckrepair(node_s,face_r);
[node_s,face_r] = surfreorient(node_s,face_r);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pts = bringPtsToSurf(surf,pts)

n = size(pts, 1);
for i=1:n
    pts(i,:) = nearestPoint(surf,pts(i,:));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p2_closest,ip2_closest] = nearestPoint(p2,p1)

% AUTHOR: Jay Dubb (jdubb@nmr.mgh.harvard.edu)

m=size(p1,1);

if(~isempty(p2) && ~isempty(p1))
    p2_closest=zeros(m,3);
    ip2_closest=zeros(m,1);
    dmin=zeros(m,1);
    for k=1:m
        d=sqrt((p2(:,1)-p1(k,1)).^2+(p2(:,2)-p1(k,2)).^2+(p2(:,3)-p1(k,3)).^2);
        [dmin(k),ip2_closest(k)]=min(d);
        p2_closest(k,:)=p2(ip2_closest(k),:);
    end
else
    p2_closest=[];
    ip2_closest = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function vol2gm = Vol2GM_Transform(VolNodes,GMNodes,radius,saveFlag)

%This function performs the mapping from a tetrahedral volume mesh to the
%associated GM mesh.  Output is in the form of a sparse transformation matrix with
%dimensions NxM where M is the number of tetrahedral mesh nodes and N is
%the number of GM surface nodes;

%Defualt radius = 3 mm
if ~exist('radius','var')
    radius = 3;
end

count = 1;
h = waitbar(0,'Calculating vol2gm transform...');
for n = 1:length(GMNodes)
    waitbar(n/length(GMNodes));
    p = GMNodes(n,:);
    ind = find(VolNodes(:,1) < (p(1)+radius) & VolNodes(:,2) < (p(2)+radius) & VolNodes(:,3) < (p(3)+radius) ...
        & VolNodes(:,1) > (p(1)-radius) & VolNodes(:,2) > (p(2)-radius) & VolNodes(:,3) > (p(3)-radius));
    nind = length(ind);
    if nind == 0
        err_message = sprintf('GM node %d has zero Head nodes within radius.  Check meshes are aligned or increase radius...',n);
        error(err_message)
        return
    end
    i(count:count+nind-1) = n;
    j(count:count+nind-1) = ind;
    s(count:count+nind-1) = 1/nind;
    count = count+nind;
end
delete(h)

vol2gm = sparse(i,j,s,length(GMNodes),length(VolNodes));

if ~exist('saveFlag','var')
    saveFlag = 1;
end

if saveFlag==1
    [mapping_filename, mapping_pathname] = uiputfile('*.mat','Save vol2gm mapping as...');
    out_full = [mapping_pathname, mapping_filename ];
    save(out_full,'vol2gm', 'radius');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


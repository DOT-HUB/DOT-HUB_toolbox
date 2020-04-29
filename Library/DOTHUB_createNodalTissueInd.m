function headVolumeMesh = DOTHUB_createNodalTissueInd(headVolumeMesh)

%This function takes the element wise tissue specification (mesh.elem(:,5))
%and uses it and the mesh to calculate the node-wise tissue specification,
%which is output as mesh.node(:,4);

% ####################### INPUTS ##########################################

% headVolumeMesh     :   A multi-layer volume mesh structure. Contains fields:
%                        node(nnodex4), face(nfacex3), elem(nelemx5),
%                        tissue_labels{1xntiss}

% ####################### OUTPUTS #########################################

% headVolumeMesh     :   A multi-layer volume mesh structure. Contains fields:
%                        node([nnodex4), face(nfacex3), elem(nelemx5),
%                        tissue_labels{1xntiss}, with node(:,4) updated.

% ####################### Dependencies ####################################

% #########################################################################
% Written by S. Brigadoi
% Tidied by RJC, April 2020
%
% ############################# Updates ###################################
% #########################################################################


%Can this be done without looping around lots of nodes? Yes, we can use
%histograms!!!! 

%[fname,pathname] = uigetfile('*.mat','Select Mesh');
%load([pathname '/' fname]);

tiss = unique(headVolumeMesh.elem(:,5));

% Use histograms to compute how frequently a node is shared by elements in
% the same tissue and then assign to that node the tissue having the
% highest frequency
tmp = zeros(size(headVolumeMesh.node,1),length(tiss));
for i_t = 1:length(tiss)
    tiss_nodes = headVolumeMesh.elem(headVolumeMesh.elem(:,5) == i_t,1:4);
    [n,~] = hist(tiss_nodes(:),1:size(headVolumeMesh.node,1));
    tmp(unique(tiss_nodes(:)),i_t) = n(unique(tiss_nodes(:)));
end
[~,ind] = max(tmp,[],2);
headVolumeMesh.node(:,4) = ind;

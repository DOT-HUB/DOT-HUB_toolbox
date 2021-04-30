
function [Fig1, Fig2]=selectseedCorrMap(dotimgFileName,gmSurfaceMesh,GMSensitivityMask);

% Input:
% dotimgFileName            % full path to the dotimg file name
% gmSurfaceMesh             % gmSurfaceMesh as a .mat structure (fields
%                               nodes and face)
% GMSensitivityMask         % Sensitivity mask calculated from the group,
%                               size 1 x # nodes x # files/subjects

% HbFlag (optional)         % set flag = 1 for HbR, 0 for HbO. Default is 0 
%                               (no flag).


if ~exist('infantFlag','var')
    HbFlag=0;
else 
    HbFlag=1
end 
    
% Output:
% Figure that plots the correlation between the averaged seed signal and 
% all other nodes in the gmSurfaceMesh  

% Set the figure
image=zeros(length(gmSurfaceMesh.node),1);
Fig1=figure; DOTHUB_plotSurfaceImage(gmSurfaceMesh,image);
  
% Activate datacursormode on the figure to allow the user to select a seed
 
disp('Select point on the surface that will serve as the seed location');
disp('Press any key when finished');

d=datacursormode;
 pause

% Collect the user selected point
 point=getfield(getCursorInfo(d), 'Position');

% Find the neared gmSurfaceMesh node to the user selected point
 [nnode, ind, offset]=DOTHUB_nearestNode(point,gmSurfaceMesh.node);

% Calculate the distance between that node and all other nodes
  nodes=gmSurfaceMesh.node;
  dists = sqrt(sum((nodes - repmat(nnode,size(nodes,1),1)).^2,2));
 
% Find which nodes fall within the user specified radius,
% create a binary mask (1 if node is within seed rgion, 0 if not)
  seedMask = dists<=10;

% Convert mask from logical to double to use for intensity
  seed=double(seedMask');

% Sum the rows of the GMSensitivityMask mask
  summedMask=sum(squeeze(GMSensitivityMask(1,:,:))');

% Load the dotimg file
  dotimg=load(dotimgFileName, '-mat');

% Multiply group sensitivity mask by the dotimg to apply the mask
if HbFlag==1
    masked=summedMask.*dotimg.hbr.gm;
else
    masked=summedMask.*dotimg.hbo.gm;
end

% Take the average values of the seed from the masked image
  seedavg=mean(masked(:,(seed'==1)),2);

% Calculate the correlation between the averaged seed value
% seed and all other nodes in the masked image. Remove any NaNs.
  corrcoefs=corr(seedavg, masked);
  corrcoefs(isnan(corrcoefs))=0;

% Plot the correlation map
 Fig2=figure; DOTHUB_plotSurfaceImage(gmSurfaceMesh,corrcoefs);

end
  


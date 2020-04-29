function [nnode, ind, offset] = DOTHUB_nearestNode(point,nodes)

%This function outputs the node (1x3) present in the list nodes (Mx3) that is closest
%to the 3D location specified by point (1x3) and the offset in whatever
%dimensions are provided.

%Inputs (point,nodes) : the single point and the list of nodes
%Outputs [nnode, offset]: the index of the specific node in the list
%'nodes' which is nearest to the input point, and the offset (euclidean
%error).

% #########################################################################
% RJC, UCL, April 2020

if size(point,2)~=3
    point = point';
    if size(point,2)~=3
        error('Dimensions of input point not 1x3');
    end
end

if size(nodes,2)~=3
    nodes = nodes';
    if size(nodes,2)~=3
        error('Dimensions of input nodes not Mx3');
    end
end

M = size(nodes,1);

if size(point,1)>1
    for i = 1:size(point,1)
        tmpPoint = point(i,:);
        dists = sqrt(sum((nodes - repmat(tmpPoint,M,1)).^2,2));
        [offset(i),ind(i)] = min(dists);
    end
else
    dists = sqrt(sum((nodes - repmat(point,M,1)).^2,2));
    [offset,ind] = min(dists);
end

nnode = nodes(ind,:);

function [source_pos,detector_pos] = SrcDetMapping(landmarks_atlas,SD,headMesh)

% SrcDetMapping.m
% Mapping of the source and detector locations recorded during the expriment
% into an atlas or subject specific mesh. It uses the landmarks to
% transform between the two spaces.
%
% The function produces the transformation matrices A and b that satisfy 
% the equation: 
%                  A*landmark_Polhemus + b = landmarks_atlas

% Then A and b are used to map sources/detector into the FEM mesh

% ************************************************************************
% UCL, EVR, 5th June 2019.

% Get affine transformation matrices
    [A,B] = affinemap(SD.Landmarks,landmarks_atlas); % Function from AtlasViewer
    
% Transformed source and detector locations    
    mappedSrc = affine_trans(SD.SrcPos,A,B);
    mappedDet = affine_trans(SD.DetPos,A,B);
    
% Find sources in the mesh
    K = dsearchn(headMesh.node(:,1:3),headMesh.elem(:,1:4),mappedSrc);
    source_pos = headMesh.node(K,1:3);
    [~,Ps] = tsearchn(headMesh.node(:,1:3),headMesh.elem(:,1:4),source_pos);
    i = find(isnan(Ps(:,1)));
    source_pos(:,2)=source_pos(:,2)+1;
    source_pos(i,2)=source_pos(i,2)+1;
    
% Find detectors in the mesh    
    K = dsearchn(headMesh.node(:,1:3),headMesh.elem(:,1:4),mappedDet);
    detector_pos = headMesh.node(K,1:3);
    [~,Pd] = tsearchn(headMesh.node(:,1:3),headMesh.elem(:,1:4),detector_pos);
    i = find(isnan(Pd(:,1)));
    detector_pos(:,2)=detector_pos(:,2)+1;
    detector_pos(i,2)=detector_pos(i,2)+1;

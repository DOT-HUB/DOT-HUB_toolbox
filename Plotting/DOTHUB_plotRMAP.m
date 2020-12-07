function DOTHUB_plotRMAP(rmap,viewAng)

%This function uses dotimg and rmap to plot a volumetric node-wise image overlayed onto a volume 
%It creates a figure per image distribution (e.g. 2 for HbO and HbR)
%
% INPUTS ##################################################################
%
% rmap          : rmap or mshs structure or path. must contain
%                 gmSurfaceMesh variable
%
% viewAng       : view angle vector 1 x 2 (optional), defaults [-37.5 30];
%
% OUTPUTS #################################################################
%
% RJC UCL, Dec 2020 #####################################################

% MANAGE VARIABLES
% #########################################################################
if ischar(rmap)
    rmapFileName = rmap;
    rmap = load(rmapFileName,'-mat');
else
    rmapFileName = rmap.fileName;
end

if ~exist('viewAng','var')
    viewAng = [-37.5 30];
elseif isempty(viewAng)
    viewAng = [-37.5 30];
end

plotmesh(rmap.SD3Dmesh.SrcPos,'r.','MarkerSize',30); hold on;
plotmesh(rmap.SD3Dmesh.DetPos,'b.','MarkerSize',30);
plotmesh(rmap.SD3Dmesh.Landmarks,'c.','MarkerSize',30);
legend('Source','Detector','Landmark','AutoUpdate','off')
plotmesh(rmap.scalpSurfaceMesh.node,rmap.scalpSurfaceMesh.face,'FaceColor',[0.9 0.85 0.75],'EdgeColor',[0.3 0.3 0.3]);
set(gcf,'Color','w');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
hold on;




function [hAxis, hPatch, hColorbar] = DOTHUB_plotSurfaceImage(mesh,intensity,viewAng,shadingtype,cmap)

% Displays node-wise distribution on mesh. Heavily influenced by Jay Dubb's 
% Homer2 function (see https://www.nitrc.org/projects/homer2/)
%
% INPUTS ##################################################################
% mesh          : mesh structure with .node and .face fields
%    
% intensity     : a distribution of length = length(mesh.node)
%
% viewAng       : view angle vector 1 x 2 (optional), defaults [-37.5 30];
%
% shadingtype   : 'interp', 'flat', 'faceted'. (Optional). Defaults to interp.
%
% cmap          : colormap. (Optional). Defaults to matlab pref.
%
% OUTPUTS #################################################################
% h         : trisurf handle
%
% RJC UCL, Jan 2020 #######################################################

if ~exist('viewAng','var')
    viewAng = [-37.5 30];
elseif isempty(viewAng)
    viewAng = [-37.5 30];
end

if ~exist('shadingtype','var')
    shadingtype = 'interp';
elseif isempty(shadingtype)
    shadingtype = 'interp';
end

nodes = mesh.node;
face  = mesh.face;
hAxis = gca;
hPatch = trisurf(face(:,1:3), nodes(:,1), nodes(:,2), nodes(:,3),intensity,'EdgeColor',[0.8 0.8 0.8],'EdgeAlpha',1);
shading(shadingtype);
set(hPatch,'diffusestrength',.7,'specularstrength',.2,'ambientstrength',.2);
set(hPatch,'Facelighting','phong');
view(viewAng);
camlight(viewAng(1),viewAng(2));
camlight(viewAng(1)+90,0);
camlight(viewAng(1)+180,0);
camlight(viewAng(1)+270,0);
axis equal;axis off;

%Set balanced colorbar;
if exist('cmap','var')
    if ischar(cmap)
        eval(['cmap = ' cmap ';']); %In case e.g. 'jet' is parsed
    end
    hAxis.Colormap = cmap;
end
hColorbar = colorbar;
climits = caxis;
limit = max(abs(climits));
caxis([-limit limit]);




function DOTHUB_movieSurfaceDOTIMG(dotimg,rmap,varargin)

% Takes dotimg and rmap/mshs files or structures and displays node-wise distribution on surface mesh
%
% INPUTS ##################################################################
%
% dotimg        : The dotimg structure or path to dotimg file. Contains
%                 hbo.gm, hbr.gm and/or mua{1}.gm, stat.gm etc.
%
% rmap          : rmap or mshs structure or path. must contain
%                 gmSurfaceMesh variable
%
%
% varargin      : input argument pairs, with options:
%
%                 'condition'   : integer specifying condition to plot
%                 'shading'     : 'interp', 'flat', 'faceted'. (Optional). Defaults to interp.
%                 'imageType'   : 'haem', 'mua', default 'haem'
%                 'colormap'    : preferred colormap array
%                 'view'        : view angle, defaults to [-37.5 30]
%
% OUTPUTS #################################################################
%
% RJC UCL, April 2020 #####################################################

% Manage Variables ########################################################
varInputs = inputParser;
varInputs.CaseSensitive = false;
validateShading = @(x) assert(any(strcmpi({'flat','interp','faceted'},x)));
validateImageType = @(x) assert(any(strcmpi({'haem','mua'},x)));
addParameter(varInputs,'shading','interp',validateShading);
addParameter(varInputs,'imageType','haem',validateImageType);
addParameter(varInputs,'colormap','greyJet');
addParameter(varInputs,'condition',1,@isnumeric);
addParameter(varInputs,'view',[-37.5 30],@isnumeric);
parse(varInputs,varargin{:});
varInputs = varInputs.Results;

viewAng = varInputs.view;
shadingtype = varInputs.shading;
condition = varInputs.condition;
hrfExplorer = false;

if ischar(dotimg)
    dotimgFileName = dotimg;
    dotimg = load(dotimgFileName,'-mat');
end

%Work out how many frames there are in dotimg.hbo.gm
%Work out max scale of any famr to fix colorbar (add scale as varargin to
%plot call.
%Run loop around frames
f1 = figure('color','w');
for frame = 1:nFrames
    DOTHUB_plotSurfaceDOTIMG(dotimg,rmap,frames,varargin)
    drawnow;
    MovToWrite(i) = getframe(f1);
end
v = VideoWriter(tit,'MPEG-4');
v.FrameRate = 10;
open(v);
writeVideo(v,MovToWrite);
close(v);
    

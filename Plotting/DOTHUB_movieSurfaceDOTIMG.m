function DOTHUB_movieSurfaceDOTIMG(dotimg,rmap,saveName,varargin)

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
% saveName      : output movie name, default is 'SurfaceDOTIMGMovie.mp4';
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

% TO DO               #####################################################
% Add input that allows sub-selection of frames within DOTIMG to moviefy

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
addParameter(varInputs,'cbScaleFactor',1,@isnumeric);
addParameter(varInputs,'title',[]);

parse(varInputs,varargin{:});
varInputs = varInputs.Results;

viewAng = varInputs.view;
shadingtype = varInputs.shading;
condition = varInputs.condition;
cbScaleFactor = varInputs.cbScaleFactor;

if ischar(dotimg)
    dotimgFileName = dotimg;
    dotimg = load(dotimgFileName,'-mat');
end

outname = 'SurfaceDOTIMGMovie.mp4';
if exist('saveName','var')
    if ~isempty(saveName)
        outname = saveName;
    end
end

%Work out how many frames we need to videoify, and work out max scale of any frame to fix colorbar (add scale as varargin to
if strcmpi(varInputs.imageType,'haem')
    nFrames = size(dotimg.hbo.gm,1);
    scalMax = max(abs([dotimg.hbo.gm(:); dotimg.hbr.gm(:)]));
elseif strcmpi(varInputs.imageType,'mua')
    nFrames = size(dotimg.mua{1}.gm,1); % is this indexing correct?
    scalMax = max(abs([dotimg.mua{1}.gm(:); dotimg.mua{2}.gm(:)]));
end

%plot call.
%Run loop around frames
MovToWrite(nFrames) = struct('cdata',[],'colormap',[]);
f1 = figure('color','w');
for frame = 1:nFrames
    DOTHUB_plotSurfaceDOTIMG(dotimg,rmap,frame,varargin{:});
    ax = findall(gcf, 'type', 'axes');
    for i = 1:length(ax)
        caxis(ax(i),[-scalMax scalMax]*cbScaleFactor);
    end
    drawnow;
    MovToWrite(frame) = getframe(f1);
end
v = VideoWriter(outname,'MPEG-4');
v.FrameRate = 10;
open(v);
writeVideo(v,MovToWrite);
close(v);
    

function DOTHUB_plotSurfaceDOTIMG(dotimg,rmap,frames,varargin)

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
% frames        : a vector specifying the frame or frames (which are then
%                 averaged) to display. If not parsed, defaults to 1
%
% varargin      : input argument pairs, with options:
%
%                 'shading'   : 'interp', 'flat', 'faceted'. (Optional). Defaults to interp.
%                 'imageType' : 'haem', 'mua', default 'haem'
%                 'colormap'  : preferred colormap array
%                 'condition' : integer specifying condition to plot
%                 'view'      : view angle, defaults to [-37.5 30]
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
cond = varInputs.condition;

if ischar(dotimg)
    dotimgFileName = dotimg;
    dotimg = load(dotimgFileName,'-mat');
end
if ischar(rmap)
    rmapFileName = rmap;
    rmap = load(rmapFileName,'-mat');
end
if ~exist('frames','var')
    frames = 1;
end
if ischar(varInputs.colormap)
    load('greyJet.mat');
    varInputs.colormap = greyJet;
end

% Define image to display and run plotting routines #######################
if strcmpi(varInputs.imageType,'haem')
    if ndims(dotimg.hbo) == 3 %Conditions exist
        img(1,:) = squeeze(mean(dotimg.hbo.gm(cond,frames,:),2));
        img(2,:) = squeeze(mean(dotimg.hbr.gm(cond,frames,:),2));
    else
        img(1,:) = squeeze(mean(dotimg.hbo.gm(frames,:),1));
        img(2,:) = squeeze(mean(dotimg.hbr.gm(frames,:),1));
    end
    nSubplot = 2;
    subplotLabels = {'HbO, \muM','HbR, \muM'};
    
else strcmpi(varInputs.imageType,'mua')
    nWavs = length(dotimg.mua);
    for i = 1:nWavs
        if ndims(dotimg.mua{i}.gm) == 3 %Conditions exist
            img(i,:) = squeeze(mean(dotimg.mua{i}.gm(cond,frames,:),2));
        else
            img(i,:) = squeeze(mean(dotimg.mua{i}.gm(frames,:),1));
        end
        subplotLabels{1,i} = [' \Delta\muA at Wav. ' num2str(i) ' mm^-^1'];
    end
    nSubplot = nWavs;
end

hFig = gcf;
set(gcf,'Color','w','Units','Normalized');
for i = 1:nSubplot
    subplot(1,nSubplot,i);
    [hAxis, hPatch, hColorbar] = DOTHUB_plotSurfaceImage(rmap.gmSurfaceMesh,img(i,:),viewAng,shadingtype,varInputs.colormap);
    set(hAxis,'FontSize',16);
    hColorbar.Location = 'South';
    tmp = hColorbar.Position;
    hColorbar.Position = [tmp(1) tmp(2) tmp(3) tmp(4)];
    hColorbar.AxisLocation = 'out';
    ylabel(hColorbar,subplotLabels{1,i})
end
[~,fname,~] = fileparts(dotimg.fileName);
sgtitle(fname,'FontSize',16,'Interpreter','none');

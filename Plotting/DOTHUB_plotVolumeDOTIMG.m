function DOTHUB_plotVolumeDOTIMG(dotimg,rmap,frames,varargin)

%This function uses dotimg and rmap to plot a volumetric node-wise image overlayed onto a volume 
%It creates a figure per image distribution (e.g. 2 for HbO and HbR)
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
%                 'condition'   : integer specifying condition to plot
%                                 (defaults to 1)
%                 'imageType'   : 'haem', 'mua', default 'haem'
%
%                 'sliceDim'    : The dimension to take the slice. If a vector, 
%                                 create a sublot for each up to N = 3;
%                                 Defaults to [1 2 3];
%                 'slicePos'    : The dist (in same units as mesh) of the slice to take for 
%                                 each dimension. Must have the same size
%                                 as sliceDim. Defaults to xyz of max(abs(image))
%                 'imageThresh' : The value below which the values of the image are
%                                 rendered transparent (e.g. a significance threshold)   
%                                 -Defaults to 0.1*max(abs(image));
%                 'cRange'      : Range applied to caxis. Defaults to
%                                 max/min
%                 'colormap'    : preferred colormap array
%
%
% OUTPUTS #################################################################
%
% RJC UCL, April 2020 #####################################################

% Manage Variables ########################################################
varInputs = inputParser;
varInputs.CaseSensitive = false;
validateImageType = @(x) assert(any(strcmpi({'haem','mua'},x)));
addParameter(varInputs,'imageType','haem',validateImageType);
addParameter(varInputs,'sliceDim',[1 2 3],@isnumeric);
addParameter(varInputs,'slicePos',[],@isnumeric);
addParameter(varInputs,'imageThresh',[],@isnumeric);
addParameter(varInputs,'cRange',[],@isnumeric);
addParameter(varInputs,'colormap','greyJet');
addParameter(varInputs,'condition',1,@isnumeric);
parse(varInputs,varargin{:});
varInputs = varInputs.Results;

sliceDim = varInputs.sliceDim;
condition = varInputs.condition;
slicePos = varInputs.slicePos;

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
cmap = varInputs.colormap;

% Define image to display and run plotting routines #######################
if strcmpi(varInputs.imageType,'haem')
    if ndims(dotimg.hbo) == 3 %Conditions exist
        img(1,:) = squeeze(mean(dotimg.hbo.vol(frames,:,condition),2));
        img(2,:) = squeeze(mean(dotimg.hbr.vol(frames,:,condition),2));
    else
        img(1,:) = squeeze(mean(dotimg.hbo.vol(frames,:,condition),1));
        img(2,:) = squeeze(mean(dotimg.hbr.vol(frames,:,condition),1));
    end
    nSubplot = 2;
    imageLabels = {'HbO, \muM','HbR, \muM'};
    
else strcmpi(varInputs.imageType,'mua')
    nWavs = length(dotimg.mua);
    for i = 1:nWavs
        if ndims(dotimg.mua{i}.vol) == 3 %Conditions exist
            img(i,:) = squeeze(mean(dotimg.mua{i}.vol(frames,:,condition),2));
        else
            img(i,:) = squeeze(mean(dotimg.mua{i}.vol(frames,:),1));
        end
        imageLabels{1,i} = [' \Delta\muA at Wav. ' num2str(i) ' mm^-^1'];
    end
    nSubplot = nWavs;
end

% Define input variables dependent on image  ##############################
if isempty(varInputs.slicePos)
    for i = 1:nSubplot
        [~,tmp] = max(img(i,:));
        slicePosAll(i,:) = rmap.headVolumeMesh.node(tmp,1:3);
    end
else
    slicePosAll = repmat(slicePos,nSubplot,1);
end

if isempty(varInputs.imageThresh)
    for i = 1:nSubplot
        imageThreshAll(i) = 0.1*max(abs(img(i,:)));
    end
end
if isempty(varInputs.cRange)
    for i = 1:nSubplot
        cRangeAll(i,:) = [-max(abs(img(i,:))) max(abs(img(i,:)))];
    end
end

for i = 1:nSubplot
    figure;
    set(gcf,'Color','w','Units','Normalized');
    slicePos = slicePosAll(1,:);
    imageThresh = imageThreshAll(i);
    cRange = cRangeAll(i,:);
    imageLabel = imageLabels{1,i};
    
    DOTHUB_plotVolumeImage(rmap.headVolumeMesh,img(i,:),sliceDim,slicePos,imageThresh,cmap,cRange,imageLabel)
end

% set(hAxis,'FontSize',16);
% hColorbar.Location = 'South';
% tmp = hColorbar.Position;
% hColorbar.Position = [tmp(1) tmp(2) tmp(3) tmp(4)];
% hColorbar.AxisLocation = 'out';
% ylabel(hColorbar,subplotLabels{1,i})
% [~,fname,~] = fileparts(dotimg.fileName);
% sgtitle(fname,'FontSize',16,'Interpreter','none');

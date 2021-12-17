function [hAxis, hPatch, hColorbar] = DOTHUB_plotSurfaceDOTIMG(dotimg,rmap,frames,varargin)

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
%                 'condition'   : integer specifying condition to plot
%                 'shading'     : 'interp', 'flat', 'faceted'. (Optional). Defaults to interp.
%                 'imageType'   : 'haem', 'mua', default 'haem'
%                 'colormap'    : preferred colormap array
%                 'view'        : view angle, defaults to [-37.5 30]
%                 'hrfExplorer' : 'on' or 'off' to add subplot showing selected node pseudochannel
%                                 default 'off';
%                 'title'       : Force title of figure
%
% OUTPUTS #################################################################
%
% RJC UCL, April 2020 #####################################################

% Manage Variables ########################################################
varInputs = inputParser;
varInputs.CaseSensitive = false;
validateShading = @(x) assert(any(strcmpi({'flat','interp','faceted'},x)));
validateImageType = @(x) assert(any(strcmpi({'haem','mua'},x)));
validatehrfExplorer = @(x) assert(any(strcmpi({'on','off'},x)));
addParameter(varInputs,'shading','interp',validateShading);
addParameter(varInputs,'imageType','haem',validateImageType);
addParameter(varInputs,'hrfExplorer',false,validatehrfExplorer);
addParameter(varInputs,'colormap','greyJet');
addParameter(varInputs,'title',[]);
addParameter(varInputs,'condition',1,@isnumeric);
addParameter(varInputs,'view',[-37.5 30],@isnumeric);
parse(varInputs,varargin{:});
varInputs = varInputs.Results;

viewAng = varInputs.view;
shadingtype = varInputs.shading;
condition = varInputs.condition;
hrfExplorer = varInputs.hrfExplorer;
title = varInputs.title;

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

conditionFlag = 0;
% Define image to display and run plotting routines #######################
if strcmpi(varInputs.imageType,'haem')
    if ndims(dotimg.hbo.gm) == 3 %Conditions exist
        conditionFlag = 1;
        img(1,:) = squeeze(mean(dotimg.hbo.gm(frames,:,condition),1));
        img(2,:) = squeeze(mean(dotimg.hbr.gm(frames,:,condition),1));
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
            conditionFlag = 1;
            img(i,:) = squeeze(mean(dotimg.mua{i}.gm(frames,:,condition),2));
        else
            img(i,:) = squeeze(mean(dotimg.mua{i}.gm(frames,:,condition),1));
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

if isempty(title)
    if conditionFlag
        sgtitle([fname ', Condition ' num2str(condition)],'FontSize',16,'Interpreter','none');
    else
        sgtitle(fname,'FontSize',16,'Interpreter','none');
    end
else
    sgtitle(title)
end

if hrfExplorer
    dcmObj = datacursormode(hFig);
    spPosition = [0.4 0.7 0.2 0.2];
    aHRF = axes('Position',spPosition);
    set(dcmObj,'UpdateFcn',{@dcmUpdateHRF,dotimg,rmap,condition,varInputs.imageType,aHRF},'Enable', 'on');
end



function out = dcmUpdateHRF(hDataTip, event_obj, dotimg, rmap, condition, imageType, aHRF)

cla(aHRF);
axes(aHRF);
out = '';

%Determine which subplot we are clicked on;
tmp = event_obj.Target.Parent.Position;
plt = round( (tmp(1)+(tmp(3)/2)) / (tmp(1)+tmp(3)));

[~,ind] = DOTHUB_nearestNode(event_obj.Position,rmap.gmSurfaceMesh.node(:,1:3));
if ndims(dotimg.hbo.gm) == 3 %Conditions exist
    if strcmpi(imageType,'haem')
        hrf(:,1) = dotimg.hbo.gm(:,ind,condition);
        hrf(:,2) = dotimg.hbr.gm(:,ind,condition);
        hrf(:,3) = sum(hrf(:,1:2),2);
        plot(dotimg.tImg,hrf(:,1),'r');hold on;
        plot(dotimg.tImg,hrf(:,2),'b');
        plot(dotimg.tImg,hrf(:,3),'g');hold off;
        xlabel('Time (s)');
        ylabel('\muM');
        legend('\DeltaHbO','\DeltaHbR','\DeltaHbT','location','best')
    else
        hrf(:,1) = dotimg.mua{1}.gm(:,ind,condition);
        hrf(:,2) = dotimg.mua{2}.gm(:,ind,condition);        
        plot(dotimg.tImg,hrf);
    end
else
    if strcmpi(imageType,'haem')
        hrf(:,1) = dotimg.hbo.gm(:,ind);
        hrf(:,2) = dotimg.hbr.gm(:,ind);
        hrf(:,3) = sum(hrf(:,1:2),2);
        plot(dotimg.tImg,hrf(:,1),'r');hold on;
        plot(dotimg.tImg,hrf(:,2),'b');
        plot(dotimg.tImg,hrf(:,3),'g');
        xlabel('Time (s)');
        ylabel('\muM');
        legend('\DeltaHbO','\DeltaHbR','\DeltaHbT','location','best')
    else
        hrf(:,1) = dotimg.mua{1}.gm(:,ind);
        hrf(:,2) = dotimg.mua{2}.gm(:,ind);        
        plot(dotimg.tImg,hrf);
    end
end

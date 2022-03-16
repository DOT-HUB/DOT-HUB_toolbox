function hFig = DOTHUB_LUMOplotArray(y,tHRF,SD,yLimits,distRange,overSizeRatios,titString,fontSize,lineWidth)

%This function is designed to plot HRFs derived from LUMO data in the 2D 
%arrangement dictated by the SD variable, with indications of the tile 
%positions. Heavily inspired by HOMER2 plotProbe.
%
%Note that this plotting routine will resize the current figure window (or
%create a new one) such that the aspect ratio of array is correct.
%Manually re-sizing the window will likely create an incorrect aspect
%ratio.  If you want a different size figure window, create a figure of the
%size you want and make it active before running this function: the width
%will be maintained and the height adjusted to get the correct proportions.
%
% INPUTS: ################################################################
% y                 = HRF data (timex3xchannels) in micromolar
% tHRF              = timebase in seconds
% SD                = SD (2D representation)
% ylimits           = The limits of the y axis [lower upper]
% distRange         = The distance range of channels to plot. Defaults to [0 60]
% overSizeRatios    = A factor by which to extend the display size of each
%                   axis [x y]
% titString         = A string to title the figure;
% fontSize          = Figure font size;
% lineWidth         = Thickness of HRF plot lines.
%
% Outputs: ###############################################################
% hFig              = Figure handle
%
%#########################################################################
% RJC, UCL, April 2020
%#########################################################################
%
% UPDATES ################################################################
% TO DO LIST #############################################################

if ~exist('distRange','var')
    distRange = [0 60];
elseif isempty(distRange)
    distRange = [0 60];
end

if ~exist('overSizeRatios','var')
    overSizeRatios = [1 1];
elseif isempty(overSizeRatios)
    overSizeRatios = [1 1];
end

if ~exist('fontSize','var')
    fontSize = 16;
elseif isempty(fontSize)
    fontSize = 16;
end

if ~exist('lineWidth','var')
    lineWidth = 1;
elseif isempty(lineWidth)
    lineWidth = 1;
end

tmpy = y(:,:,SD.MeasListAct(1:end/2)==1);
mnVal = min(tmpy(:));
mxVal = max(tmpy(:));
fprintf(['Maximum and minimum data values = [' num2str(mnVal,'%3f') ' ' num2str(mxVal,'%3f') ']  µM\n']);

if ~exist('yLimits','var')
    yl(1) = mnVal;
    yl(2) = mxVal;
elseif isempty(yLimits)
    yl(1) = mnVal;
    yl(2) = mxVal;
else
    yl(1) = ylimits(1);
    yl(2) = ylimits(2);
end

if yl(2)<1e-4
    warning('The data values seem too small to be in micro-molar units.');
end

%Define chanpos and scaling;
nchan = size(SD.MeasList,1)/2;
for i = 1:nchan
    pos(i,:) = mean([SD.SrcPos(SD.MeasList(i,1),1:2);SD.DetPos(SD.MeasList(i,2),1:2)]);
end
centre = mean(pos);

%Determine hex corner positions
nTiles = SD.nSrcs/3;
ang = 0;
hexCorners(1,:) = [0 14.5 0];
for i = 2:6
    hexCorners(i,:) = rotz(hexCorners(i-1,:),60,[0 0 0]);
end

for i = 1:nTiles
    u = SD.SrcPos((i-1)*3+1,:) - SD.DetPos((i-1)*4+3,:); %Src 1 of each tile - Det 3 of each tile - the center
    v = [0 1 0];
    CosTheta = dot(u,v)/(norm(u)*norm(v));
    theta = acosd(CosTheta);
    for j = 1:6
        tmp(j,:) = rotz(hexCorners(j,:),theta,[0 0 0]) + SD.DetPos((i-1)*4+3,:);
    end
    hexCornersArray((i-1)*6+1:i*6,:) = tmp;
end
hexCornersArray = hexCornersArray(:,1:2);

%Define maxima etc.
allposOrig = [pos;hexCornersArray];
Xupper = max(pos(:,1));
Xlower = min(pos(:,1));
Yupper = max(pos(:,2));
Ylower = min(pos(:,2));

nPlotsX = length(unique(pos(:,1)));
nPlotsY = length(unique(pos(:,2)));
%nPlotsX_border = nPlotsX + 4;
%nPlotsY_border = nPlotsY + 4;

plotWidth = 8*overSizeRatios(1)/nPlotsX; % There is no correct way of doing this because it is array dependent. 
plotHeight = 5*overSizeRatios(2)/nPlotsY; % These values provide reasonable proportions for LUMO12 array

%Normalise positioning values for plotting
pos_cent = pos - repmat(centre,size(pos,1),1);
hexCornersArray = hexCornersArray - repmat(centre,size(hexCornersArray,1),1);

allpos = [pos_cent;hexCornersArray];
pos_cent(:,1) = pos_cent(:,1) + abs(min(allpos(:,1)));
pos_cent(:,2) = pos_cent(:,2) + abs(min(allpos(:,2)));
hexCornersArray(:,1) = hexCornersArray(:,1) + abs(min(allpos(:,1)));
hexCornersArray(:,2) = hexCornersArray(:,2) + abs(min(allpos(:,2)));

allpos = [pos_cent;hexCornersArray];
pos_cent(:,1) = pos_cent(:,1)./max(allpos(:,1));
pos_cent(:,2) = pos_cent(:,2)./max(allpos(:,2));
hexCornersArray(:,1) = hexCornersArray(:,1)./max(allpos(:,1));
hexCornersArray(:,2) = hexCornersArray(:,2)./max(allpos(:,2));

borderWidth = 0.1; %Screen ratio of border width each edge
borderScal = 1 - 2*borderWidth;
%Scale
pos_cent(:,1) = pos_cent(:,1)*(borderScal);
pos_cent(:,2) = pos_cent(:,2)*(borderScal);
hexCornersArray(:,1) = hexCornersArray(:,1)*(borderScal);
hexCornersArray(:,2) = hexCornersArray(:,2)*(borderScal);

allpos = [pos_cent;hexCornersArray];
%Re-center x
pos_cent(:,1) = pos_cent(:,1) - (DOTHUB_range(allpos(:,1))/2) +  min(allpos(:,1)) + 0.5;
hexCornersArray(:,1) = hexCornersArray(:,1) - (DOTHUB_range(allpos(:,1))/2) +  min(allpos(:,1)) + 0.5;
%Re-center y with offset for gnomons
pos_cent(:,2) = pos_cent(:,2) - (DOTHUB_range(allpos(:,2))/2) +  min(allpos(:,2)) + 0.55;
hexCornersArray(:,2) = hexCornersArray(:,2) - (DOTHUB_range(allpos(:,2))/2) +  min(allpos(:,2)) + 0.55;

%Pos_cent is now the centre position of each axis, need to convert this to
%the position from left and bottom
pos_cent(:,1) = pos_cent(:,1) - (plotWidth/2);
pos_cent(:,2) = pos_cent(:,2) - (plotHeight/2);

hFig = gcf;
clf;

set(gcf,'Color','w');
%First draw hexagons;
axes('Position',[0 0 1 1]);
for i = 1:nTiles
    for j = 2:6
        line([hexCornersArray((i-1)*6+j-1,1) hexCornersArray((i-1)*6+j,1)],[hexCornersArray((i-1)*6+j-1,2) hexCornersArray((i-1)*6+j,2)],'color',[0.6 0.6 0.6]);
    end
    line([hexCornersArray((i-1)*6+j,1) hexCornersArray((i-1)*6+1,1)],[hexCornersArray((i-1)*6+j,2) hexCornersArray((i-1)*6+1,2)],'color',[0.6 0.6 0.6]);
end
xlim([0 1]);
ylim([0 1]);
set(gca,'Position',[0 0 1 1]);
axis off

if exist('titString','var')
    tHand = title(titString,'FontSize',16);axis off;
    tHand.Position = [0.5 0.95 0];
end

dists = DOTHUB_getSDdists(SD);
for i = 1:nchan
    
    if SD.MeasListAct(i)==1 & dists(i)>= distRange(1) & dists(i)<=distRange(2)
        
        axHand(i) = axes;
        plot(tHRF, squeeze(y(:,1,i)), 'r','LineWidth',lineWidth);hold on;
        plot(tHRF, squeeze(y(:,2,i)), 'b','LineWidth',lineWidth);
        plot(tHRF, squeeze(y(:,3,i)), 'g','LineWidth',lineWidth);
        ylim(yl);
        xlim([tHRF(1) tHRF(end)]);
        set(axHand(i),'Position',[pos_cent(i,1) pos_cent(i,2) plotWidth plotHeight],'Color','none');
        axis off
        if i==1
            yT = get(gca,'YTick');
        end
        %Add manual lines/ticks
        %line([tHRF(1) tHRF(end)],[0 0],'Color','k');
        %line([0 0],[yT(1) yT(end)],'Color','k');
    end
end

%Add gnomons - try to keep these the same as the plot axes, but warn if
%there is no room to fit them in the border.
plotAR = plotWidth/plotHeight;
if plotHeight > borderWidth
    warning('Gnomon height exceeds border');
    gnomonHeight = 0.8*borderWidth;
    gnomonWidth = plotAR*gnomonHeight;
else
    gnomonHeight = plotHeight;
    gnomonWidth = plotWidth;
end

axGnomon = axes;
ylim(yl);
xlim([tHRF(1) tHRF(end)]);
set(axGnomon,'Position',[0.4 0.05 gnomonWidth gnomonHeight],'LineWidth',3,'FontSize',10,'YTick',[],'XTick',[]);
xlabel([num2str(DOTHUB_range(tHRF),'%2.1f') ' s']);
ylabel([num2str(sum(abs(yl)),'%2.1f') ' µM']);

%Add Legend
axLegend = axes;
set(axLegend,'Position',[0.6-gnomonWidth 0.05 gnomonWidth gnomonHeight],'FontSize',fontSize);
plot(0,0, 'r','LineWidth',lineWidth*2);hold on;
plot(0,0, 'b','LineWidth',lineWidth*2);
plot(0,0, 'g','LineWidth',lineWidth*2);
ylim([-0.5 0.5]);
xlim([-0.5 0.5]);
legHand = legend({'HbO','HbR','HbT'});
set(legHand,'Position',[0.6-gnomonWidth 0.05 gnomonWidth gnomonHeight]);
axis off

%Size figure to ensure hexagon aspect ratio is correct.
h = DOTHUB_range(allposOrig(:,2))/DOTHUB_range(allposOrig(:,1));
tmp = get(hFig,'Position');
tmp(4) = tmp(3)*h;
set(hFig,'Position',tmp);

function hFig = DOTHUB_LUMOplotArray(y,tHRF,SD,yLimits,distRange,overSizeRatios,titString,fontSize,lineWidth)

%This function is designed to plot HRFs derived from LUMO data in the 2D 
%arrangement dictated by the SD variable, with indications of the tile 
%positions. Heavily inspired by HOMER2 plotProbe.
%
%Note that the 2D SD, while never a perfect representation of a real world
%array, should be scalled approximately correctly (tiles sizes correct,
%etc). Otherwise this routine will likely fail
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
    lineWidth = 1.5;
elseif isempty(lineWidth)
    lineWidth = 1.5;
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
dists = DOTHUB_getSDdists(SD);
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

%nPlotsX = length(unique(xPosSort));
%nPlotsY = length(unique(yPosSort));
%nPlotsX_border = nPlotsX + 4;
%nPlotsY_border = nPlotsY + 4;

% First Centre around (0,0)
pos_cent = pos - repmat(centre,size(pos,1),1); 
hexCornersArray = hexCornersArray - repmat(centre,size(hexCornersArray,1),1);

%Put all in positive quadrant
allpos = [pos_cent;hexCornersArray];
pos_cent(:,1) = pos_cent(:,1) + abs(min(allpos(:,1)));
pos_cent(:,2) = pos_cent(:,2) + abs(min(allpos(:,2)));
hexCornersArray(:,1) = hexCornersArray(:,1) + abs(min(allpos(:,1)));
hexCornersArray(:,2) = hexCornersArray(:,2) + abs(min(allpos(:,2)));

%Estimate the required width and height of the plots on the basis of
% a) define occupied space in x and y as total range of x and y channel position values minus
% any gaps larger than 60 mm
% b) finding the number of unique channel positions (rounded to nearest 5
% mm) in both x and y
% c) number of plots in X and Y = occupued length/nPlots.
gapThresh = 60;
posCropRound = (round(pos(dists<gapThresh,:)./5))*5;

xPosSort = sort(unique(posCropRound(:,1)));
dfX = diff(xPosSort);
xOccupied = range(xPosSort) - dfX(dfX>gapThresh);

yPosSort = sort(unique(posCropRound(:,1)));
dfY = diff(yPosSort);
yOccupied = range(yPosSort) - dfY(dfY>gapThresh);

%find maximum number of channels at any given y or x
uniY = unique(yPosSort);
nPlotsX = 0;
for i = 1:length(uniY)
    tmp = length(unique(posCropRound(posCropRound(:,2)==uniY(i),1)));
    if tmp>nPlotsX
        nPlotsX = tmp;
    end
end

uniX = unique(xPosSort);
nPlotsY = 0;
for i = 1:length(uniX)
    tmp = length(unique(posCropRound(posCropRound(:,1)==uniX(i),2)));
    if tmp>nPlotsY
        nPlotsY = tmp;
    end
end

%Normalize units
allpos = [pos_cent;hexCornersArray];
pos_cent(:,1) = pos_cent(:,1)./max(allpos(:,1)); %In positive quadrant so max = range
pos_cent(:,2) = pos_cent(:,2)./max(allpos(:,2));
hexCornersArray(:,1) = hexCornersArray(:,1)./max(allpos(:,1)); %This is where things get tricky.
hexCornersArray(:,2) = hexCornersArray(:,2)./max(allpos(:,2));
aspRot = max(allpos(:,1))/max(allpos(:,2)); %width/height

plotWidthNative = (xOccupied/max(allpos(:,1)))/nPlotsX; %Normalized units
plotHeightNative = (yOccupied/max(allpos(:,2)))/nPlotsY;
plotWidth = overSizeRatios(1)*plotWidthNative; 
plotHeight = overSizeRatios(2)*plotHeightNative;  

%Add border
borderWidth = 0.1; %Screen ratio of border width at each edge
borderScal = 1 - 2*borderWidth;

%Scale
pos_cent(:,1) = pos_cent(:,1)*(borderScal);
pos_cent(:,2) = pos_cent(:,2)*(borderScal);
hexCornersArray(:,1) = hexCornersArray(:,1)*(borderScal);
hexCornersArray(:,2) = hexCornersArray(:,2)*(borderScal);

%Re-center x
pos_cent(:,1) = pos_cent(:,1) + borderWidth;
hexCornersArray(:,1) = hexCornersArray(:,1) + borderWidth;

%Re-center y with offset of borderWidth/2 for gnomon
pos_cent(:,2) = pos_cent(:,2) + borderWidth + borderWidth/2;
hexCornersArray(:,2) = hexCornersArray(:,2) + borderWidth + borderWidth/2;
allpos = [pos_cent;hexCornersArray];

%Pos_cent is now the centre position of each channel in axis units. Need to convert this to
%the position from left and bottom corner to set axis position?
%Perhaps make the display such that the mid point of the x axis (y=0)
%aligns with the channel location, in which case we need to consider yl
pos_cent(:,1) = pos_cent(:,1) - (plotWidth/2);
pos_cent(:,2) = pos_cent(:,2) - (abs(yl(1))/(yl(2)-yl(1)))*plotHeight;

%Sort out figure and aspect ratio
hFig = gcf;
clf(hFig);
set(hFig,'Color','w');
%Need to figure out aspect ratio
set(hFig,'Units','Normalized','Position',[0 0 1 1]);
set(hFig,'Units','Pixels');
screenWidth = hFig.Position(3);
screenHeight = hFig.Position(4);
if aspRot > 1 %target width>height. Make width the smaller of screenHeight or screenwidth, and set height accordingly
    set(hFig,'Position',[hFig.Position(1) hFig.Position(2) min([screenWidth screenHeight]) min([screenWidth screenHeight])/aspRot]);
else %target width < height. Make height the smaller of screenHeight or screenwidth , and set width accordingly
    set(hFig,'Position',[hFig.Position(1) hFig.Position(2) min([screenWidth screenHeight])*aspRot min([screenWidth screenHeight])]);
end
set(hFig,'Units','Normalized');

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
    
    if SD.MeasListAct(i)==1 && dists(i)>= distRange(1) && dists(i)<=distRange(2)
        
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

set(hFig,'Units','Normalized');
%Resize figure to ensure hexagon aspect ratio is correct? 
% Haven't figured this out yet
% h = DOTHUB_range(allposOrig(:,2))/DOTHUB_range(allposOrig(:,1));
% tmp = get(hFig,'Position');
% tmp(4) = tmp(3)*h;
% set(hFig,'Position',tmp);


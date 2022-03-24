function DOTHUB_plotIntVDist(d,SD,labelFlag,xAxisUpperLim,varargin)

% This function creates an intensity versus SD-distance drop-off curve
% figure.
%
% ######################## INPUTS ##########################################
%
% d             : The intensity data matrix (timepoints x channels) (Homer2 style)
% SD            : The source-detector structure (Homer2 style). Should be 3D.
% labelFlag     : (Optional). If true, the datacursor is activated to allow point
%                 labelling (optional, default false);
% xAxisUpperLim : (Optional). The upper limit on the x axis. Defaults to
%                 max distance in array, but cropping useful for papers.
%
% varargin  =  optional input pairs:
%              'hAxes' - optional axis handle
%              'noiseFloor' - flag to estimate and plot noise floor. Default true.
%
% ######################## OUTPUTS #########################################
%
% Outputs are figures...
%
% ######################## Dependencies ####################################
% This script requires other functions in the DOTHUB function library
%
% #########################################################################
% RJC, UCL, April 2020
%
% ############################# Updates ###################################
% #########################################################################

if ~exist('labelFlag','var')
    labelFlag = 0;
end

varInputs = inputParser;
addParameter(varInputs,'hAxes','',@ishandle);
validateNFFlag = @(x) assert(any(strcmpi({'on','off'},x)));
addParameter(varInputs,'noiseFloor','on',validateNFFlag);
parse(varInputs,varargin{:});
varInputs = varInputs.Results;
if isempty(varInputs.hAxes)
    varInputs.hAxes = gca;
end
if isempty(varInputs.noiseFloor)
    noiseFloorFlag = true;
else
    if strcmpi(varInputs.noiseFloor,'on')
        noiseFloorFlag = true;
    elseif strcmpi(varInputs.noiseFloor,'off')
        noiseFloorFlag = false;
    end
end
      
dists = DOTHUB_getSDdists(SD);
dists = [dists dists];
nWavs = length(SD.Lambda);

if ~exist('xAxisUpperLim','var')
    xAxisUpperLim = max(dists)+1;
end

if size(d,1)==1 %if only one sample
    mnD = d;
else
    mnD = mean(d);
end

markersize = 40;
Markers = {'o','s','*','x','d'};
count = 1;
for i = 1:nWavs
    scatter(varInputs.hAxes,dists(SD.MeasList(:,4)==i),mnD(SD.MeasList(:,4)==i),markersize,[0 0 1],Markers{i},'LineWidth',1);
    hold(varInputs.hAxes,'on');
    legText{count} = ['All chan. ' num2str(SD.Lambda(i)) ' nm'];
    count = count+1;
    scatter(varInputs.hAxes,dists(SD.MeasListAct==1 & (SD.MeasList(:,4)==i)),mnD(SD.MeasListAct==1 & (SD.MeasList(:,4)==i)),markersize,[1 0 0],Markers{i},'LineWidth',1);
    legText{count} = ['Good chan. ' num2str(SD.Lambda(i)) ' nm'];
    count = count+1;
end

xlim(varInputs.hAxes,[0 xAxisUpperLim]);
set(varInputs.hAxes,'YScale','log','XGrid','on','YGrid','on','box','on','FontSize',16);
xlabel(varInputs.hAxes,'S-D Distance (mm)');
ylabel(varInputs.hAxes,'Intensity (arb.)');

if noiseFloorFlag
    if max(dists)>70
        noisefloorest = mean(mnD(dists>70));
    else
        noisefloorest = min(mnD(:));
    end
    line(varInputs.hAxes,[0 xAxisUpperLim],[noisefloorest noisefloorest],'LineWidth',2,'LineStyle','-.','Color','k');
    text(varInputs.hAxes,1,noisefloorest*0.75,['Noise floor ~ ' num2str(noisefloorest,'%0.2e')]);
end

hold(varInputs.hAxes,'off');
legend(varInputs.hAxes,legText);

% Determine if data is LUMO data (to allow module labelling) 
% ##########################################################################
% ##### (Temporary solution? Need to define module elements in .SD3D?) #####
if SD.nSrcs == (3/4)*SD.nDets && rem(SD.nSrcs,3)==0 && rem(SD.nDets,4)==0 &&min(dists)<12
    lumoFlag = 1;
else
    lumoFlag = 0;
end
% ##########################################################################

if labelFlag
    dcm = datacursormode(gcf);
    datacursormode on
    if lumoFlag
        set(dcm,'updatefcn',{@DOTHUB_LUMOintvdistplot_tiplabels,SD})
    else
        set(dcm,'updatefcn',{@DOTHUB_intvdistplot_tiplabels,SD})
    end
    dcm.updateDataCursors
end

function labels = DOTHUB_LUMOintvdistplot_tiplabels(obj,event_obj,SD)

index = event_obj.Target.Children.DataIndex; 
%this could be channel index or index of channel(goodchan)
%determine if it is a point from all chan, good chan, and define
%wavelength;
if strcmpi(event_obj.Target.Marker,'square')
    wav = '850nm';
else
    wav = '735nm';
end
if event_obj.Target.CData(1) == 1 %Is red, so index is relative to good data
    tmp = 1:size(SD.MeasList,1);
    tmp2 = tmp(SD.MeasListAct==1);
    index = tmp2(index);
end
   
tmpS = mod(SD.MeasList(:,1),3); tmpS(tmpS==0) = 3;
srcLabs = {'A','B','C'};
tmpD = mod(SD.MeasList(:,2),4); tmpD(tmpD==0) = 4;
si = SD.SrcPos(SD.MeasList(index,1),:);
di = SD.DetPos(SD.MeasList(index,2),:);
labels = {['Channel ',num2str(index),', dist = ',num2str(sqrt(sum((si-di).^2)),4),' mm'],...
          ['Source = ' num2str(SD.MeasList(index,1)) ' (Tile ' num2str(ceil(SD.MeasList(index,1)/3)) ', Src = ' srcLabs{tmpS(index)} ', Wav = ' wav ')'],...
          ['Detector = ' num2str(SD.MeasList(index,2)) ' (Tile ' num2str(ceil(SD.MeasList(index,2)/4)) ', Detector = ' num2str(tmpD(index)) ')']};

      
function labels = DOTHUB_intvdistplot_tiplabels(obj,event_obj,SD)

index = event_obj.Target.Children.DataIndex; 
%this could be channel index or index of channel(goodchan)
%determine if it is a point from all chan, good chan, and define
%wavelength;
if strcmpi(event_obj.Target.Marker,'square')
    wav = 2;
else
    wav = 1;
end
if event_obj.Target.CData(1) == 1 %Is red, so index is relative to good data
    tmp = 1:size(SD.MeasList,1);
    tmp2 = tmp(SD.MeasListAct==1);
    index = tmp2(index);
end
   
si = SD.SrcPos(SD.MeasList(index,1),:);
di = SD.DetPos(SD.MeasList(index,2),:);
labels = {['Channel ',num2str(index),', dist = ',num2str(sqrt(sum((si-di).^2)),4),' mm'],...
          ['Source = ' num2str(SD.MeasList(index,1))],...
          ['Detector = ' num2str(SD.MeasList(index,2))]};



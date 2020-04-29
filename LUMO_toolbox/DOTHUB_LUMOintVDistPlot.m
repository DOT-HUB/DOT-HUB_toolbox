function DOTHUB_LUMOintVDistPlot(d,SD,labelFlag,xAxisUpperLim)

% This function creates an intensity versus SD-distance drop-off curve
% figure.

%######################## INPUTS ##########################################

% d             : The intensity data matrix (timepoints x channels) (Homer2 style)
% SD            : The source-detector structure (Homer2 style). Should be 3D.
% labelFlag     : (Optional). If true, the datacursor is activated to allow point
%                 labelling (optional, default false);
% xAxisUpperLim : (Optional). The upper limit on the x axis. Defaults to
%                 max distance in array, but cropping useful for papers.

%######################## OUTPUTS #########################################

%Outputs are figures...

%######################## Dependencies ####################################
%This script requires other functions in the DOTHUB function library

% #########################################################################
% RJC, UCL, April 2020
%
% ############################# Updates ###################################
% #########################################################################

if ~exist('labelFlag','var')
    labelFlag = 0;
end

dists = DOTHUB_getSDdists(SD);
dists = [dists dists];

if ~exist('xAxisUpperLim','var')
    xAxisUpperLim = max(dists)+1;
end

if size(d,1)==1 %if only one sample
    mnD = d;
else
    mnD = mean(d);
end

markersize = 40;
scatter(dists(SD.MeasList(:,4)==1),mnD(SD.MeasList(:,4)==1),markersize,[0 0 1],'LineWidth',1);hold on;
scatter(dists(SD.MeasList(:,4)==2),mnD(SD.MeasList(:,4)==2),markersize,[0 0 1],'Marker','s','LineWidth',1); %Second wav must be square for datatips to work.
scatter(dists(SD.MeasListAct==1 & (SD.MeasList(:,4)==1)),mnD(SD.MeasListAct==1 & (SD.MeasList(:,4)==1)),markersize,[1 0 0],'LineWidth',1);
scatter(dists(SD.MeasListAct==1 & (SD.MeasList(:,4)==2)),mnD(SD.MeasListAct==1 & (SD.MeasList(:,4)==2)),markersize,[1 0 0],'Marker','s','LineWidth',1);

ylim([1e-6 5]);
xlim([0 xAxisUpperLim]);
set(gca,'YScale','log','XGrid','on','YGrid','on','FontSize',16);
xlabel('S-D Distance (mm)');
ylabel('Intensity (arb.)');
legend(['All chan. ' num2str(SD.Lambda(1)) ' nm'],['Good chan. ' num2str(SD.Lambda(1)) ' nm'],['All chan. ' num2str(SD.Lambda(2)) ' nm'],['Good chan. ' num2str(SD.Lambda(2)) ' nm']);
box on

if max(dists)>70
    noisefloorest = mean(mnD(dists>70));
    line([0 xAxisUpperLim],[noisefloorest noisefloorest],'LineWidth',2,'LineStyle','-.','Color','k');
    text(1,noisefloorest*0.75,['Noise floor ~ ' num2str(noisefloorest,'%0.2e')]);
end

hold off;

if labelFlag
    dcm = datacursormode(gcf);
    datacursormode on
    set(dcm,'updatefcn',{@DOTHUB_LUMOintvdistplot_tiplabels,SD})
    dcm.updateDataCursors
end

function labels = DOTHUB_LUMOintvdistplot_tiplabels(obj,event_obj,SD)

index = event_obj.Target.Children.DataIndex; 
%this could be channel index or index of channel(goodchan)
%determine if it is a point from all chan, good chan, and define
%wavelength;
if strcmpi(event_obj.Target.Marker,'square');
    wav = 2;
else
    wav = 1;
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
          ['Source = ' num2str(SD.MeasList(index,1)) ' (Tile ' num2str(ceil(SD.MeasList(index,1)/3)) ', Src = ' srcLabs{tmpS(index)} ')'],...
          ['Detector = ' num2str(SD.MeasList(index,2)) ' (Tile ' num2str(ceil(SD.MeasList(index,2)/4)) ', Detector = ' num2str(tmpD(index)) ')']};



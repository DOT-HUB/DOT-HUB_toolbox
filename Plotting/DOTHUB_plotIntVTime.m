function DOTHUB_plotIntVTime(t,d,SD,labelFlag,varargin)

% This function creates an intensity versus time figure with active labels
% for LUMO
%
%######################## INPUTS ##########################################
%
% t             : time vector (Homer2 style)
% d             : The intensity data matrix (timepoints x channels) (Homer2 style)
% SD            : The source-detector structure (Homer2 style). Note that
%                 this plotting function will exclude channels listed as zero in
%                 SD.MeasListAct. If you want all data displayed, parse a
%                 clearn SD.MeasListAct
% labelFlag     : If true, the datacursor is activated to allow point
%                 labelling (optional, default off);
%
% varargin  =  optional input pairs:
%              'hAxes' - optional axis handle
%              'noiseFloor' - flag to estimate and plot noise floor. Default true.
%
%######################## OUTPUTS #########################################
%
%Outputs are figures...
%
%######################## Dependencies ####################################
%This script requires other functions in the DOTHUB function library
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

plot(varInputs.hAxes,t,d(:,SD.MeasListAct==1));
xlim([min(t) max(t)]);
set(varInputs.hAxes,'YScale','log','XGrid','on','YGrid','on','FontSize',16);
xlabel('Time (s)');
ylabel('Intensity (arb.)');
box on

if noiseFloorFlag
    if max(dists)>70
        mnD = mean(d);
        noisefloorest = mean(mean(d(:,dists>70)));
    else
        noisefloorest = min(mnD(:));
    end
    line([min(t) max(t)],[noisefloorest noisefloorest],'LineWidth',2,'LineStyle','-.','Color','k');
    text(1,noisefloorest*0.75,['Noise floor ~ ' num2str(noisefloorest,'%0.2e')]);
end

hold off;

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
        set(dcm,'updatefcn',{@DOTHUB_LUMOintvtimeplot_tiplabels,SD,d})
    else
        set(dcm,'updatefcn',{@DOTHUB_intvtimeplot_tiplabels,SD,d})    
    end
    dcm.updateDataCursors
end

function labels = DOTHUB_LUMOintvtimeplot_tiplabels(~,event_obj,SD,d)

tindex = event_obj.DataIndex;
values = event_obj.Position; %This will be the values of the clicked point (t,d)
index = find(d(tindex,:)==values(2) & SD.MeasListAct'==1,1,'first');

tmpS = mod(SD.MeasList(:,1),3); tmpS(tmpS==0) = 3;
srcLabs = {'A','B','C'};
tmpD = mod(SD.MeasList(:,2),4); tmpD(tmpD==0) = 4;
si = SD.SrcPos(SD.MeasList(index,1),:);
di = SD.DetPos(SD.MeasList(index,2),:);
if SD.MeasList(index,4)== 1
    lambda = SD.Lambda(1);
elseif SD.MeasList(index,4)== 2
    lambda = SD.Lambda(2);
end
labels = {['Channel ',num2str(index),', dist = ',num2str(sqrt(sum((si-di).^2)),4),' mm'],...
          ['Source = ' num2str(SD.MeasList(index,1)) ' (Tile ' num2str(ceil(SD.MeasList(index,1)/3)) ', Src = ' srcLabs{tmpS(index)} ...
          ', Wavelength = ',num2str(lambda),' nm)'],...
          ['Detector = ' num2str(SD.MeasList(index,2)) ' (Tile ' num2str(ceil(SD.MeasList(index,2)/4)) ', Detector = ' num2str(tmpD(index)) ')']};



function labels = DOTHUB_intvtimeplot_tiplabels(~,event_obj,SD,d)

tindex = event_obj.DataIndex;
values = event_obj.Position; %This will be the values of the clicked point (t,d)
index = find(d(tindex,:)==values(2) & SD.MeasListAct'==1,1,'first');

si = SD.SrcPos(SD.MeasList(index,1),:);
di = SD.DetPos(SD.MeasList(index,2),:);
labels = {['Channel ',num2str(index),', dist = ',num2str(sqrt(sum((si-di).^2)),4),' mm'],...
          ['Source = ' num2str(SD.MeasList(index,1))],...
          ['Detector = ' num2str(SD.MeasList(index,2))]};


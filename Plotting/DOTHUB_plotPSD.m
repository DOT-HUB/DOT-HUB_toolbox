function DOTHUB_plotPSD(t,d,SD,labelFlag)

% This function creates a PSD plot with active labels
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

dists = DOTHUB_getSDdists(SD);
dists = [dists dists];

fs = length(t)/DOTHUB_range(t);
[psd,freq] = periodogram(d,rectwin(size(d,1)),size(d,1),fs,'psd');
psddb = 10*log10(psd);
plot(freq,psddb(:,SD.MeasListAct==1));
xlim([min(freq) max(freq)]);
set(gca,'YScale','linear','XGrid','on','YGrid','on','FontSize',16);
xlabel('Freq (Hz))');
ylabel('Power/freq. (dB/Hz)');
box on
hold off;

% Determine if data is LUMO data (to allow module labelling) 
% ##########################################################################
% ##### (Temporary solution? Need to define module elements in .SD3D?) #####
if SD.nSrcs == (3/4)*SD.nDets && rem(SD.nSrcs,3)==0 && rem(SD.nDets,4)==0 && min(dists)<12
    lumoFlag = 1;
else
    lumoFlag = 0;
end
% ##########################################################################

if labelFlag
    dcm = datacursormode(gcf);
    datacursormode on
    if lumoFlag
        set(dcm,'updatefcn',{@DOTHUB_LUMOpsdplot_tiplabels,SD,psddb})
    else
        set(dcm,'updatefcn',{@DOTHUB_psdplot_tiplabels,SD,psddb})
    end
    dcm.updateDataCursors
end

function labels = DOTHUB_LUMOpsdplot_tiplabels(~,event_obj,SD,psddb)

findex = event_obj.DataIndex;
values = event_obj.Position; %This will be the values of the clicked point (t,d)
index = find(psddb(findex,:)==values(2) & SD.MeasListAct'==1,1,'first');

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


function labels = DOTHUB_psdplot_tiplabels(~,event_obj,SD,psddb)

findex = event_obj.DataIndex;
values = event_obj.Position; %This will be the values of the clicked point (t,d)
index = find(psddb(findex,:)==values(2) & SD.MeasListAct'==1,1,'first');

si = SD.SrcPos(SD.MeasList(index,1),:);
di = SD.DetPos(SD.MeasList(index,2),:);
labels = {['Channel ',num2str(index),', dist = ',num2str(sqrt(sum((si-di).^2)),4),' mm'],...
          ['Source = ' num2str(SD.MeasList(index,1))],...
          ['Detector = ' num2str(SD.MeasList(index,2))]};

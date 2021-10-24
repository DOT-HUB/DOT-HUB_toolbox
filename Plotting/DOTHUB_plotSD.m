function DOTHUB_plotSD(SD,SDRange,labelFlag)
 
%This function plots an SD layout in 3D, with source, detector, channels
%and landmarks (if present).  All channels in SD.MeasListAct are plotted,
%unless filtered via SDrange
%
%############################### INPUTS ###################################
%
% SD      :   Homer2 style SD structure (2D or 3D) or path thereto.
%
% SDrange :   (Optional) Two-element vector of specified distange range of channels to
%             include in plot.  Default to [0 inf];
%
%############################# Dependencies ###############################
%
% #########################################################################
% RJC, UCL, May 2020
% 
% ############################## Updates ##################################
% #########################################################################
%
% ############################### TO DO ###################################
% #########################################################################

if ~exist('SDRange','var')
    SDRange = [0 inf];
elseif isempty(SDRange)
    SDRange = [0 inf];
end

if ~exist('labelFlag','var')
    labelFlag = false;
end

if ischar(SD)
    tmp = load(SD,'-mat');
    fn = fieldnames(tmp);
    eval(['SD = tmp.' fn{1}]);
end

if ~isfield(SD,'MeasListAct')
    SD.MeasListAct = true(size(SD.MeasList,1),1);
end

plotmesh(SD.SrcPos,'r.','MarkerSize',30); hold on;
plotmesh(SD.DetPos,'b.','MarkerSize',30);


dists = DOTHUB_getSDdists(SD);
dists = [dists dists];

flag = 1;
for i = 1:size(SD.MeasList,1)
    if dists(i) > SDRange(1) && dists(i)<=SDRange(2) && SD.MeasListAct(i)==1
        lh = line([SD.SrcPos(SD.MeasList(i,1),1) SD.DetPos(SD.MeasList(i,2),1)], ...
            [SD.SrcPos(SD.MeasList(i,1),2) SD.DetPos(SD.MeasList(i,2),2)],...
            [SD.SrcPos(SD.MeasList(i,1),3) SD.DetPos(SD.MeasList(i,2),3)],'LineWidth',2.5);
        lh.Color=[0,1,0,0.2];
        if flag
            legend('Source','Detector','Channel','AutoUpdate','off')
            flag = 0;
        end  
    end
end

if isfield(SD,'Landmarks')
    plotmesh(SD.Landmarks,'c.','MarkerSize',30);
end

if labelFlag
    for i = 1:SD.nSrcs
        plotmesh(SD.SrcPos(i,:),'r.','MarkerSize',30); hold on;
        text(SD.SrcPos(i,1)+2,SD.SrcPos(i,2)+2,SD.SrcPos(i,3)+2,['S' num2str(i)],'Color','r');
    end
    for i = 1:SD.nDets
        plotmesh(SD.DetPos(i,:),'b.','MarkerSize',30); hold on;
        text(SD.DetPos(i,1)+2,SD.DetPos(i,2)+2,SD.DetPos(i,3)+2,['D' num2str(i)],'Color','b');
    end
    if isfield(SD,'Landmarks')
        plotmesh(SD.Landmarks,'c.','MarkerSize',30);
    end
else
    plotmesh(SD.SrcPos,'r.','MarkerSize',30); hold on;
    plotmesh(SD.DetPos,'b.','MarkerSize',30);
    if isfield(SD,'Landmarks')
        plotmesh(SD.Landmarks,'c.','MarkerSize',30);
    end
end

hold off;
axis equal;
set(gca,'FontSize',16);

        



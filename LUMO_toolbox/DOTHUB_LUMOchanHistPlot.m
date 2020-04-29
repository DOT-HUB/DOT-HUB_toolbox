function DOTHUB_LUMOchanHistPlot(SD,xAxisUpperLim)

% This function creates a histogram of good channels (defined by SD.MeasListAct) 
% overlayed on all possible channels;
%
%######################## INPUTS ##########################################
%
% SD            : The source-detector structure (Homer2 style). Should be 3D.
% xAxisUpperLim : (Optional). The upper limit on the x axis. Defaults to
%                 62.5 mm, but cropping can be useful for papers.
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

dists = DOTHUB_getSDdists(SD);
dists = [dists dists];
SD.MeasListAct(end/2+1:end) = SD.MeasListAct(1:end/2);
goodChan = SD.MeasListAct==1;

if ~exist('xAxisUpperLim','var')
    xAxisUpperLim = 62.5;
end

%Plotting
hg = histogram(dists);
hg.BinEdges = [7.5:5:xAxisUpperLim];
hg.FaceColor = 'b';
hg.FaceAlpha = 1;
hold on;
histogram(dists(goodChan),hg.BinEdges,'FaceColor','r','FaceAlpha',1);
xlabel('Source-detector distance (mm)');
ylabel('Channel Count');
legend('All Channels','Good Channels');
set(gca,'FontSize',16,'Xgrid','on','Ygrid','on');
xlim([0 xAxisUpperLim]);
box on

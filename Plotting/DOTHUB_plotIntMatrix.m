function DOTHUB_plotIntMatrix(d,SD)

% This function creates a matrix plot of log intensity
%
%######################## INPUTS ##########################################
%
% d             : The intensity data matrix (timepoints x channels) (Homer2 style)
% SD            : The source-detector structure (Homer2 style).
%
%######################## OUTPUTS #########################################
%
%Outputs are figures...
%
%######################## Dependencies ####################################
%This script requires other functions in the DOTHUB function library
%
% #########################################################################
% RJC, UCL, Nov 2020
%
% ############################# Updates ###################################
% #########################################################################

dmn = log10(mean(d,1));
dMat1 = ones(SD.nSrcs,SD.nDets)*-10;
dMat2 = ones(SD.nSrcs,SD.nDets)*-10;

for i = 1:SD.nSrcs
    dInd = SD.MeasList(SD.MeasList(:,1)==i & SD.MeasList(:,4)==1,2);
    dMat1(i,dInd) = dmn(SD.MeasList(:,1)==i & SD.MeasList(:,4)==1);
    dInd = SD.MeasList(SD.MeasList(:,1)==i & SD.MeasList(:,4)==2,2);
    dMat2(i,dInd) = dmn(SD.MeasList(:,1)==i & SD.MeasList(:,4)==2);
end

ax1 = subplot(1,2,1);
set(gcf,'color','w');
imagesc(ax1,dMat1);axis square
caxis(ax1,[-8 0]);
cb1 = colorbar(ax1,'northoutside');
ylabel(cb1,'-log10(intensity) 735nm','FontSize',12)
colormap(ax1,'hot')
set(ax1,'Box','on','FontSize',12,'XTick',1:SD.nDets,'YTick',1:SD.nSrcs);
axis(ax1, 'tight');
xlabel(ax1,'Detector','FontSize',12);
ylabel(ax1,'Source','FontSize',12);

ax2 = subplot(1,2,2);
imagesc(ax2,dMat2);axis square
caxis(ax2,[-8 0]);
cb2 = colorbar(ax2,'northoutside');
ylabel(cb2,'-log10(intensity) 850nm','FontSize',12)
colormap(ax2,'hot')
set(ax2,'Box','on','FontSize',12,'XTick',1:SD.nDets,'YTick',1:SD.nSrcs);
axis(ax2, 'tight');
xlabel(ax2,'Detector','FontSize',12);
ylabel(ax2,'Source','FontSize',12);

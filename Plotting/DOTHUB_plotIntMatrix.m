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

dMat = log10(mean(d,1));
dMat1 = reshape(dMat(1:end/2),SD.nDets,SD.nSrcs)';
dMat2 = reshape(dMat(end/2+1:end),SD.nDets,SD.nSrcs)';

ax1 = subplot(1,2,1);
set(gcf,'color','w');
imagesc(ax1,dMat1);axis square
caxis(ax1,[-6 0]);
cb1 = colorbar(ax1,'northoutside');
ylabel(cb1,'-log10(intensity) 735nm','FontSize',12)
colormap(ax1,'gray')
set(ax1,'Box','on','FontSize',12,'XTick',1:SD.nDets,'YTick',1:SD.nSrcs);
axis(ax1, 'tight');
xlabel(ax1,'Detector','FontSize',12);
ylabel(ax1,'Source','FontSize',12);

ax2 = subplot(1,2,2);
imagesc(ax2,dMat2);axis square
caxis(ax2,[-6 0]);
cb2 = colorbar(ax2,'northoutside');
ylabel(cb2,'-log10(intensity) 850nm','FontSize',12)
colormap(ax2,'gray')
set(ax2,'Box','on','FontSize',12,'XTick',1:SD.nDets,'YTick',1:SD.nSrcs);
axis(ax2, 'tight');
xlabel(ax2,'Detector','FontSize',12);
ylabel(ax2,'Source','FontSize',12);

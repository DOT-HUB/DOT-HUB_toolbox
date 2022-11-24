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

dMatOrig = log10(mean(d,1));

nWav = max(SD.MeasList(:,4));
nRowsPerWav = size(d,2)/nWav;
for i = 1:nWav
    dMat = reshape(dMatOrig((i-1)*nRowsPerWav+1:i*nRowsPerWav),SD.nDets,SD.nSrcs)';

    ax(i) = subplot(1,nWav,i);
    set(gcf,'color','w');
    imagesc(ax(i),dMat);axis square
    caxis(ax(i),[-6 0]);
    cb1 = colorbar(ax(i),'northoutside');
    ylabel(cb1,['-log10(intensity)' num2str(SD.Lambda(i)) 'nm'],'FontSize',12)
    colormap(ax(i),'gray')
    set(ax(i),'Box','on','FontSize',12,'XTick',1:SD.nDets,'YTick',1:SD.nSrcs);
    axis(ax(i), 'tight');
    xlabel(ax(i),'Detector','FontSize',12);
    ylabel(ax(i),'Source','FontSize',12);

end

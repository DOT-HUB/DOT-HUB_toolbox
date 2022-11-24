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
% Updated by GCVL 2022 Cambridge to handle n wavelengths
% #########################################################################




%How Rob would do it ###################
% dMatOrig = log10(mean(d,1));
% 
% nWav = max(SD.MeasList(:,4));
% nRowsPerWav = size(d,2)/nWav;
% for i = 1:nWav
%     dMat = reshape(dMatOrig((i-1)*nRowsPerWav+1:i*nRowsPerWav),SD.nDets,SD.nSrcs)';
% 
%     ax(i) = subplot(1,nWav,i);
%     set(gcf,'color','w');
%     imagesc(ax(i),dMat);axis square
%     caxis(ax(i),[-6 0]);
%     cb1 = colorbar(ax(i),'northoutside');
%     ylabel(cb1,['-log10(intensity)' num2str(SD.Lambda(i)) 'nm'],'FontSize',12)
%     colormap(ax(i),'gray')
%     set(ax(i),'Box','on','FontSize',12,'XTick',1:SD.nDets,'YTick',1:SD.nSrcs);
%     axis(ax(i), 'tight');
%     xlabel(ax(i),'Detector','FontSize',12);
%     ylabel(ax(i),'Source','FontSize',12);
% 
% end
%###################



dMat = log10(mean(d,1));

if length(SD.Lambda) == 2
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

else
    % Assume dMats are being split by wavelength
    % Do some messy splitting work to get five wavelengths (according to order
    % in spreadsheet)
    n_sourcetiles = SD.nSrcs;
    n_detectortiles = SD.nSrcs;
    n_detectorspertile = 4;
    
    n_channels = n_detectorspertile * n_detectortiles * n_sourcetiles;
    n_wavs = length(SD.Lambda);
    n_dataptstotal = n_wavs * n_channels;
    
    
    r = n_sourcetiles*n_detectorspertile;
    
    % Sorting order: column = wavelength, row = channel number
    % (mathematical proof for this covering n wavelengths available)
    for p = 1:n_wavs
        data(:,p) = [((p-1)*r)+1:((p-1)*r)+r ((p-1)*r)+(n_dataptstotal/2)+1:((p-1)*r)+(n_dataptstotal/2)+r];
    end
    % Assign the indexing to the actual data points
    for j = 1:n_wavs
        for k = 1:length(data)
            d_new(k,j) = dMat(data(k,j));
        end
    end
    % Split into cells where each cell is one wavelength and contains an
    % array of Y x X where Y is the number of sources and X the number of
    % detectors
    for l = 1:n_wavs
        dMat_multi{l} = d_new(:,l);
        dMat_multi{l} = reshape(dMat_multi{l},[r,n_sourcetiles]);
    end
    
    % Plotting
    for h = 1:n_wavs
        ax1 = subplot(2,3,h);
        set(gcf,'color','w');
        imagesc(ax1,dMat_multi{h});%axis square
        caxis(ax1,[-6 0]);
        cb1 = colorbar(ax1,'northoutside');
        ylabel(cb1,'-log10(intensity)','FontSize',12)
        colormap(ax1,'gray')
        %set(ax1,'Box','on','FontSize',12,'XTick',1:SD.nDets,'YTick',1:SD.nSrcs);
        
        % I think potentially this has been the wrong way around? And x
        % axis is sources while y axis is detectors? Certainly works when
        % I change round the xtick and ytick - will tentatively change this
        % going forwards
        set(ax1,'Box','on','FontSize',12,'XTick',1:SD.nSrcs,'YTick',1:SD.nDets);
        axis(ax1, 'tight');
        xlabel(ax1,'Source','FontSize',12);
        ylabel(ax1,'Detector','FontSize',12);
        title(['Source ',num2str(SD.Lambda(h)) 'nm'])
        hold on
    end
    hold off
end

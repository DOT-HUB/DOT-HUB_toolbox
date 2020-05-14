function DOTHUB_plotHRF(timebase,data,errdata,blockendtime,fontsize)
 
%This function plots a single-channel HRF, with shaded error bars etc and
%other embellishments 

%############################### INPUTS ###################################

% timebase      =   a vector of time points for the HRF

% data          =   hrf data of dimensions time x chromophore [oxy deoxy tot]
%                   in micromolar units!

% errdata       =   error on the hrf data with same dimensions as data (optional)

% blockendtime  =   timebase value at which stimualtion trial/block ends
%                   (optional). Assued to start at timebase = 0. If parsed, block period is shaded.

% fontsize      =   axis fontsize (optional). Defaults to 16.

%############################# Dependencies ###############################
%shadedErrorBar.m by Rob Campbell https://github.com/raacampbell/shadedErrorBar

% #########################################################################
% RJC, UCL, April 2020
% 
% ############################## Updates ##################################
% #########################################################################

% ############################### TO DO ###################################
% #########################################################################

if ~exist('errdata','var')
    errFlag = 0;
else
    errFlag = 1;
end

if ~exist('blockendtime','var')
    blockFlag = 0;
else
    blockFlag = 1;
end

if ~exist('fontsize','var')
    fontsize = 16;
end

if size(data,1) ~= size(timebase,1) %Correct timebase orientation if needed
    timebase = timebase';
end

if max(data(:))<1e-4;
    warning('Are you sure you parsed data in micromolar units?');
end

fs = round(1/mean(diff(timebase)));

if errFlag  
    plot(timebase, data(:,1), 'r','LineWidth',1.5);hold on;
    plot(timebase, data(:,2), 'b','LineWidth',1.5);
    plot(timebase, data(:,3), 'g','LineWidth',1.5);
    legend('HbO', 'HbR', 'HbT');
    
    shadedErrorBar(timebase, data(:,3),errdata(:,3),{'g','LineWidth',1.5},0.6);    
    shadedErrorBar(timebase, data(:,2),errdata(:,2),{'b','LineWidth',1.5},0.6);
    shadedErrorBar(timebase, data(:,1),errdata(:,1),{'r','LineWidth',1.5},0.6);
    
    xlabel('Time (s)','fontsize',14);
    ylabel('Concentration change (\muM)','fontsize',fontsize);
    set(gca,'fontsize',fontsize);
    hold off;   
else
    plot(timebase, data(:,1), 'r','LineWidth',1.5);hold on;
    plot(timebase, data(:,2), 'b','LineWidth',1.5);
    plot(timebase, data(:,3), 'g','LineWidth',1.5);
    xlabel('Time (s)','fontsize',fontsize);
    ylabel('Concentration change (\muM)','fontsize',fontsize);
    legend('HbO', 'HbR', 'HbT');
    set(gca,'fontsize',fontsize);
    hold off;
end

xlim([min(timebase) max(timebase)]);
yl = ylim;

if blockFlag
    hold on
    pt = patch([0 0 blockendtime blockendtime],[yl(1)-1 yl(2)+1 yl(2)+1 yl(1)-1],[0.9 0.9 0.9]);
    pt.EdgeColor = 'None';
    pt.FaceAlpha = 0.5;
    set(gca,'children',flipud(get(gca,'children')))
    hold off
    ylim(yl);
end
%line([0 0],[yl(1)-1 yl(2)+1],'LineWidth',1.5,'Color','k');
%ylim(yl);


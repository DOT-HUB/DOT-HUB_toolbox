function DOTHUB_dataQualityCheck(nirsFileName,printFigFlag,timeSelectFlag)

%This script reads a .nirs file and provides a series of figures
%demonstrating the quality of the data recorded in that file. This is to
%allow rapid checking of the quality of a recording. The results are NOT
%saved into the .nirs file.  Pruning must be repeated in the downstream
%pre-processing steps.
%
%######################## INPUTS ##########################################
%
% nirsFilename:     The path of the .nirs data file
%
% printFigFlag:     (Optional). If true, the figures are printed to the directory
%                   of the nirs file. Default is true.

% timeSelectFlag:   (Optional). If true, user will be prompted to manually
%                   select a time period from the recording for analysis.
%                   Default is false
%
%######################## OUTPUTS #########################################
%
% Outputs are figures and a _dataQualityCheck.txt file
%
%######################## Dependencies ####################################
% This script requires the DOTHUB function library, export_fig, and homer2.
%
% #########################################################################
% RJC, UCL, April 2020
%
% ############################# Updates ###################################
% #########################################################################
% nirsFileName = '/Users/RCooper/Dropbox/Projects/LUMO/Lumo_AdultVisual_Pilot_150419/Rob_150419/Rob_150419_Visual2.nirs';
% #########################################################################
% #########################################################################

% MANAGE VARIABLES  ##################################################
if ~exist('nirsFileName','var')
    disp('Select .nirs file...');
    [file,path] = uigetfile('*.nirs','Select .nirs file');
    nirsFileName = fullfile(path,file);
end

[fpath,fname,fext] = fileparts(nirsFileName);
if strcmpi(fext,'.LUMO')
    error('Please convert .LUMO to .nirs first...');
end

if ~exist('printFigFlag','var')
    printFigFlag = 0;
end

if ~exist('timeSelectFlag','var')
    timeSelectFlag = 0;
end

%Load data ############################################################
load(nirsFileName,'-mat');

%Determine SD distances  ##############################################
if ~exist('SD3D','var')
    warning('SD3D not found in .nirs file, taking distance estimates from SD')
    SDtmp = SD;
else
    SDtmp = SD3D;
end
dists = DOTHUB_getSDdists(SDtmp);

% MOTION CHECK  #######################################################
%timeSelectFlag allows the user to select a clean period (motion free) for analysis.
%A subset of channels are plotted to decrease computer load (this can
%always be changed to plot all channels). % For manual time period selection, click once at the beginning of the segment to
% select the start point. Drag the cursor and click again to select the end
% point. Hit enter when finished. After, the figure will reappear with the
% selected time period highlighted in blue. This is to make sure that the correct period(s) was(were)
%selected beforre proceeding. If the period is correct, type 'y' and hit enter.
if (timeSelectFlag==1)
    my_ans = 'n';
    while my_ans == 'n'
        fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 0.8 0.8], 'color', [1 1 1])
        if size(d,2) >= 500
            plot(d(:,randi(size(d, 2),[500, 1])));
        else
            plot(d);
        end
        xlim([0 size(d,1)]); ylim([0 1.5]); title('Select segment to include (start-end pair) - Press enter when finished', 'fontsize', 20)
        xlabel ('Time (samples)','Fontsize', 15); ylabel('Intensity (raw)', 'Fontsize', 15)
        set(gca, 'Yscale', 'log')
        
        hold on
        [x, ~] = ginput;
        % If you select segments
        if  ~isempty(x)
            % Repeat if x data is not an even number
            if  mod(size(x,1),2) ~= 0
                close all; clc
                warning('Select start-end of each segments, should be pairs of numbers (even)')
                clear x
                continue
            end
            % If first selected point is before 0 assume 1
            if x(1)<0
                x(1) = 1;
            end
            % If last selected point is bigger than size data, assume is the last point
            if x(end)>size(d,1)
                x(end) = size(d,1);
            end
            
            mrk_data = zeros(size(x,1)/2,2);
            mrk_data(:,1) = ceil(x(1:2:size(x,1)));
            mrk_data(:,2) = ceil(x(2:2:size(x,1)));
            
            % Then Check outcome
            close all
            
            fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 0.8 0.8], 'color', [1 1 1])
            if size(d,2) >= 500
                plot(d(:,randi(size(d, 2),[500, 1])));
            else
                plot(d);
            end
            xlim([0 size(d,1)]); title('Select segment to include (start-end pair) - Press enter when finished', 'fontsize', 14)
            
            hold on
            for i = 1: size(mrk_data,1)
                area([mrk_data(i,1) mrk_data(i,2)], [1.5 1.5], 'FaceColor', 'c',...
                    'LineStyle', 'none'); alpha 0.2
                hold on
            end
            xlabel ('Time (samples)','Fontsize', 15);ylabel('Optical density', 'Fontsize', 15)
            set(gca,'fontsize',15); set(gca, 'Yscale', 'log')
            
            % Ask if ok and continue, or not ok and repeat selection
            prompt = {'Selection correct? Yes = y , No = n '};
            dlgtitle = 'Input';
            dims = 1;
            definput = {'y','hsv'};
            my_ans = inputdlg(prompt,dlgtitle,dims,definput);
            my_ans = my_ans{1};
            
            if my_ans == 'n'
                clear x mrk_data
                close all
                clc
                warning('Repeat selection')
                continue
            end
            
            clean_idx = zeros(size(d,1),1);
            for i = 1: size(mrk_data,1)
                clean_idx(mrk_data(i,1):mrk_data(i,2)) = 1;
            end
            % if noisy segments are not selected
        else
            clean_idx = ones(size(d,1),1);
            my_ans = 'y';
        end
        dcrop=d(clean_idx==1, :);
        
        close all
    end
end

% - this script below does not focus on motion content, but providing
%an estimate is very useful to allow the other metrics to be evaluated on
%non-motion-corrupted data. Problem is, this is very system specific, so
if timeSelectFlag==0
    fs = length(t)/DOTHUB_range(t);
    SDtmp.MeasListAct = ones(size(SDtmp.MeasList,1),1);
    SDtmp.MeasListAct = ([dists dists] <20)'; %look for motion in short channels only
    tMotionArtifact = hmrMotionArtifact(d,fs,SDtmp,ones(length(t),1),1,1,10,0.5);
    %plot(t,d);
    %hold on
    %plot(t,tMotionArtifact);
    
    %Find longest motion-free segment;
    fntMot = find(tMotionArtifact==0);
    ind = diff(fntMot);
    [lngth,loc] = max(ind);
    goodRange = [fntMot(loc)+1 fntMot(loc)+lngth];
    if isempty(goodRange) || (goodRange(2)-goodRange(1) < fs*5) %if there is no goodRange, or goodRange is less than 5 seconds in length
        dcrop = d;
    else
        dcrop = d(goodRange(1):goodRange(2),:);
    end
end

% SNR PRUNE  ##############################################################
dRange = [0 2.5];
SNRthresh = 12;
SDrange = [0 1e6];
SDtmp.MeasListAct = ones(size(SDtmp.MeasList,1),1);
SDtmp = enPruneChannels(dcrop,SDtmp,ones(size(dcrop,1),1),dRange,SNRthresh,SDrange,0);
SDtmp.MeasListAct(end/2+1:end) = SDtmp.MeasListAct(1:end/2);

%Create all-channels SDclean
SDclean = SDtmp;
SDclean.MeasListAct = ones(size(SD.MeasListAct));

% PLOTTING  ###############################################################
% First plot intvtime plot for all data (better for fault finding)
f1 = figure('Units','Normalized','Position',[0 0 0.8 0.8],'Color','w');
f1.Color = 'w';
DOTHUB_plotIntVTime(t,d,SDclean,1)
title(fname,'Interpreter','none','FontSize',16,'FontWeight','Bold');

% Plot PSD of user-selected data if timeSelectFlag=1 or automatically
% detected clean data if timeSelectFlag=0
f2 = figure('Units','Normalized','Position',[0 0 0.8 0.8],'Color','w');
DOTHUB_plotPSD(t,dcrop,SDclean,1)
title(fname,'Interpreter','none','FontSize',16,'FontWeight','Bold');

% Plot Intensity versus Distance scatter plot and histogram
f3 = figure('Units','Normalized','Position',[0 0 0.8 0.8],'Color','w');
subplot(1,2,1);
DOTHUB_plotIntVDist(dcrop,SDtmp,1);
subplot(1,2,2);
DOTHUB_plotChanHist(SDtmp);
sgtitle(fname,'Interpreter','none','FontSize',16,'FontWeight','Bold');

% Plot greyscale intensity matrix
f4 = figure('Units','Normalized','Position',[0 0 0.8 0.8],'Color','w');
DOTHUB_plotIntMatrix(d,SD)

% PRINT SUMMARY ###########################################################
diary(fullfile(fpath,[fname '_dataQualityCheck.txt']));
noisefloorest = mean(mean(dcrop(:,(dists>70))));
if timeSelectFlag==0
    fprintf(['Estimated motion burden = ' num2str(100*sum(~tMotionArtifact)*(1/fs)/DOTHUB_range(t),'%0.1f') ' %% of recording\n']);
end
fprintf(['Estimated noise floor = ' num2str(noisefloorest,'%0.2e') '\n']);
rng = [0 20];
tmp = length(find(SDtmp.MeasListAct(1:end/2)==1 & dists'>rng(1) & dists'<rng(2))) / length(find(dists'>rng(1) & dists'<rng(2)));
fprintf(['Percentage of good channels in range [' num2str(rng(1)) '  ' num2str(rng(2)) '] = ' num2str(tmp*100,'%0.1f') '%% \n']);
rng = [20 27.5];
tmp = length(find(SDtmp.MeasListAct(1:end/2)==1 & dists'>rng(1) & dists'<rng(2))) / length(find(dists'>rng(1) & dists'<rng(2)));
fprintf(['Percentage of good channels in range [' num2str(rng(1)) ' ' num2str(rng(2)) '] = ' num2str(tmp*100,'%0.1f') '%% \n']);
rng = [27.5 32.5];
tmp = length(find(SDtmp.MeasListAct(1:end/2)==1 & dists'>rng(1) & dists'<rng(2))) / length(find(dists'>rng(1) & dists'<rng(2)));
fprintf(['Percentage of good channels in range [' num2str(rng(1)) ' ' num2str(rng(2)) '] = ' num2str(tmp*100,'%0.1f') '%% \n']);
rng = [32.5 37.5];
tmp = length(find(SDtmp.MeasListAct(1:end/2)==1 & dists'>rng(1) & dists'<rng(2))) / length(find(dists'>rng(1) & dists'<rng(2)));
fprintf(['Percentage of good channels in range [' num2str(rng(1)) ' ' num2str(rng(2)) '] = ' num2str(tmp*100,'%0.1f') '%% \n']);
rng = [37.5 42.5];
tmp = length(find(SDtmp.MeasListAct(1:end/2)==1 & dists'>rng(1) & dists'<rng(2))) / length(find(dists'>rng(1) & dists'<rng(2)));
fprintf(['Percentage of good channels in range [' num2str(rng(1)) ' ' num2str(rng(2)) '] = ' num2str(tmp*100,'%0.1f') '%% \n']);
rng = [42.5 inf];
tmp = length(find(SDtmp.MeasListAct(1:end/2)==1 & dists'>rng(1) & dists'<rng(2))) / length(find(dists'>rng(1) & dists'<rng(2)));
fprintf(['Percentage of good channels in range [' num2str(rng(1)) ' ' num2str(rng(2)) '] = ' num2str(tmp*100,'%0.1f') '%% \n']);
fprintf(['Total number of good fNIRS channels = ' num2str(sum(SDtmp.MeasListAct(1:end/2))) '/' num2str(length(SDtmp.MeasListAct(1:end/2))) '\n']);
diary off;

% SAVE FIGURES ######
if printFigFlag
    if isempty(fpath)
        fpath = pwd;
    end
    fprintf(['Printing figures to ' fpath '\n']);
    
    filename = [fpath '/' fname '_intvtimeplot.jpg'];
    print(f1,filename,'-djpeg','-r400');
    filename = [fpath '/' fname '_psdplot.jpg'];
    print(f2,filename,'-djpeg','-r400');
    filename = [fpath '/' fname '_intvdistplot.jpg'];
    print(f3,filename,'-djpeg','-r400');
    filename = [fpath '/' fname '_intMatrix.jpg'];
    print(f4,filename,'-djpeg','-r400');
end

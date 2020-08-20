function DOTHUB_dataQualityCheck(nirsFileName,dRange,SNRthresh,printFigFlag)

%This script reads a .nirs file and provides a series of figures
%demonstrating the quality of the data recorded in that file. This is to
%allow rapid checking of the quality of a recording. The results are NOT
%saved into the .nirs file.  Pruning must be repeated in the downstream
%pre-processing steps.
%
%######################## INPUTS ##########################################
%
% nirsFilename: The path of the .nirs data file
%
% printFigFlag: (Optional). If true, the figures are printed to the directory 
%               of the nirs file. Default is true.
%               
%
%######################## OUTPUTS #########################################
%
%Outputs are figures and a _dataQualityCheck.txt file
%
%######################## Dependencies ####################################
%This script requires the DOTHUB function library, export_fig, and homer2.
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
    [file,path] = uigetfile('*.nirs','Select .nirs file');
    nirsFileName = fullfile(path,file);
end

[fpath,fname,fext] = fileparts(nirsFileName);
if strcmpi(fext,'.LUMO')
    error('Please convert .LUMO to .nirs first...');
end

if ~exist('dRange','var')
    dRange = [0 1e15];
end

if ~exist('SNRthresh','var')
    SNRthresh = 12;
end

if ~exist('printFigFlag','var')
    printFigFlag = 1;
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
% - this script does not focus on motion content. Providing
% an estimate is very useful to allow the other metrics to be evaluated on
% non-motion-corrupted data, but this is very hard to achieve in a way that
% is hardware/system agnostic.

% SNR PRUNE  ##############################################################
SDrange = [0 inf];
SDtmp = enPruneChannels(d,SDtmp,ones(size(d,1),1),dRange,SNRthresh,SDrange,0);
SDtmp = DOTHUB_balanceMeasListAct(SDtmp);

%Create all-channels SDclean
SDclean = SDtmp;
SDclean.MeasListAct = ones(size(SD.MeasListAct));

% PLOTTING  ###############################################################
% First plot intvtime plot for all data (better for fault finding)
f1 = figure('Units','Normalized','Position',[0 0 0.8 0.8],'Color','w');
f1.Color = 'w';
DOTHUB_plotIntVTime(t,d,SDclean,1)
title(fname,'Interpreter','none','FontSize',16,'FontWeight','Bold');

% Plot PSD of all data
f2 = figure('Units','Normalized','Position',[0 0 0.8 0.8],'Color','w');
subplot(1,2,1);
DOTHUB_plotPSD(t,d,SDclean,1);
subplot(1,2,2);
DOTHUB_plotPSD(t,d,SDtmp,1);
title(fname,'Interpreter','none','FontSize',16,'FontWeight','Bold');

% Plot channel-wise intvdist plots
f3 = figure('Units','Normalized','Position',[0 0 0.8 0.8],'Color','w');
subplot(1,2,1);
DOTHUB_plotIntVDist(d,SDtmp,1)
subplot(1,2,2);
DOTHUB_plotChanHist(SDtmp);
sgtitle(fname,'Interpreter','none','FontSize',16,'FontWeight','Bold');

% PRINT SUMMARY ###########################################################
diary(fullfile(fpath,[fname '_dataQualityCheck.txt']));
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
end

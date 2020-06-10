function DOTHUB_dataQualityCheck(nirsFileName,printFigFlag)

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
% - this script does not focus on motion content, but providing
%an estimate is very useful to allow the other metrics to be evaluated on
%non-motion-corrupted data. Problem is, this is very system specific, so
fs = length(t)/range(t);
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
dcrop = d(goodRange(1):goodRange(2),:);

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

% Plot PSD of all data
f2 = figure('Units','Normalized','Position',[0 0 0.8 0.8],'Color','w');
DOTHUB_plotPSD(t,dcrop,SDclean,1)
title(fname,'Interpreter','none','FontSize',16,'FontWeight','Bold');

% Plot channel-wise intvdist plots
f3 = figure('Units','Normalized','Position',[0 0 0.8 0.8],'Color','w');
subplot(1,2,1);
DOTHUB_plotIntVDist(dcrop,SDtmp,1);
subplot(1,2,2);
DOTHUB_plotChanHist(SDtmp);
sgtitle(fname,'Interpreter','none','FontSize',16,'FontWeight','Bold');

% PRINT SUMMARY ###########################################################
diary(fullfile(fpath,[fname '_dataQualityCheck.txt']));
noisefloorest = mean(mean(dcrop(:,(dists>70))));
fprintf(['Estimated motion burden = ' num2str(100*sum(~tMotionArtifact)*(1/fs)/range(t),'%0.1f') ' %% of recording\n']);
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
end

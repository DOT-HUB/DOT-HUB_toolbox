function [SD2D,SD3D] = DOTHUB_LUMOlayout2SD(layoutFilename,plotFlag,saveFlag)

%This function converts the LUMO .json layout file to SD. It assumes all
%docks are populated (as it must), and uses the 3D data to determine the
%channel list on the basis of all channels <60 mm existing.

%############################# INPUTS #####################################

% jsoNFilename :    The path of the .json layout file to convert to SD.
%                   If layoutFile is not parsed (or parsed empty)
%                   it will be requested.

%############################# OUTPUTS ####################################

%Outputs are SD structures for the 2D and 3D positioning information. These
%are also saved as .SD files with the same name as the .json.

%############################ Dependencies ################################
% #########################################################################

% RJC, UCL, Summer 2019

% ############################# Updates ###################################
% #########################################################################

% ############################## TO DO: ###################################
% #########################################################################

% MANAGE VARIABLES  #######################################################
if ~exist('layoutFilename','var')
    [filename, pathname, ~] = uigetfile('*.json','Select .json layout file');
    layoutFilename = [pathname '/' filename];
elseif isempty(layoutFilename)
    [filename, pathname, ~] = uigetfile('*.json','Select .json layout file');
    layoutFilename = [pathname '/' filename];
end

if ~exist('plotFlag','var')
    plotFlag = 1;
end
if ~exist('saveFlag','var')
    saveFlag = 1;
end

%Load data
layoutData = jsondecode(fileread(layoutFilename));

%Target SD output names
outname2D = [layoutFilename(1:end-5) '_2D.SD'];
outname3D = [layoutFilename(1:end-5) '_3D.SD'];

%Wavelengths
wavelength1 = 735;
wavelength2 = 850;

%distLimit
distLimit = 1e6; %mm

%#####
nNodes = length(layoutData.docks);
SD.nSrcs = nNodes*3;
SD.nDets = nNodes*4;

% Now determine optode positions from 2D information in layout JSON file
SD.SrcPos = zeros(SD.nSrcs,3);
SD.DetPos = zeros(SD.nDets,3);
for n = 1:nNodes
    nid = n;
    for det = 1:4
        SD.DetPos(det+(n-1)*4,1) = layoutData.docks(nid).optodes(det).coordinates_2d.x;
        SD.DetPos(det+(n-1)*4,2) = layoutData.docks(nid).optodes(det).coordinates_2d.y;
    end
    for src = 1:3
        SD.SrcPos(src+(n-1)*3,1) = layoutData.docks(nid).optodes(src+4).coordinates_2d.x;
        SD.SrcPos(src+(n-1)*3,2) = layoutData.docks(nid).optodes(src+4).coordinates_2d.y;
    end
end

SD.Lambda = [wavelength1 wavelength2];
SD.SpatialUnit = 'mm';

%#########
count = 1;
for i = 1:SD.nSrcs
    for j = 1:SD.nDets
        %Insert channel distance limit filter?
        distTmp = sqrt(sum((SD.SrcPos(i,:) - SD.DetPos(j,:)).^2));
        if distTmp <= distLimit
            dists(count) = distTmp; %might be useful later
            SD.MeasList(count,1) = i;
            SD.MeasList(count,2) = j;
            SD.MeasList(count,3) = 0;
            SD.MeasList(count,4) = 1;
            count = count+1;
        end
    end
end
SD.MeasList = [SD.MeasList; SD.MeasList];
SD.MeasList(end/2 + 1:end,4) = 2;

SD2D = SD;

%####################################### Create 3D SD too ####################################
SD3D = SD;
for n = 1:nNodes
    nid = n;
    for det = 1:4
        SD3D.DetPos(det+(n-1)*4,1) = layoutData.docks(nid).optodes(det).coordinates_3d.x;
        SD3D.DetPos(det+(n-1)*4,2) = layoutData.docks(nid).optodes(det).coordinates_3d.y;
        SD3D.DetPos(det+(n-1)*4,3) = layoutData.docks(nid).optodes(det).coordinates_3d.z;
    end
    for src = 1:3
        SD3D.SrcPos(src+(n-1)*3,1) = layoutData.docks(nid).optodes(src+4).coordinates_3d.x;
        SD3D.SrcPos(src+(n-1)*3,2) = layoutData.docks(nid).optodes(src+4).coordinates_3d.y;
        SD3D.SrcPos(src+(n-1)*3,3) = layoutData.docks(nid).optodes(src+4).coordinates_3d.z;
    end
end

% Display #######################################
if plotFlag;
    figure;
    subplot(1,2,1);
    for i = 1:length(SD.MeasList)
        line([SD.SrcPos(SD.MeasList(i,1),1) SD.DetPos(SD.MeasList(i,2),1)],[SD.SrcPos(SD.MeasList(i,1),2) SD.DetPos(SD.MeasList(i,2),2)],'LineWidth',2,'Color','c');
        hold on;
    end
    plot(SD.SrcPos(:,1),SD.SrcPos(:,2),'r.','MarkerSize',20);
    plot(SD.DetPos(:,1),SD.DetPos(:,2),'b.','MarkerSize',20);
    axis equal;
    
    subplot(1,2,2);
    for i = 1:length(SD.MeasList)
        line([SD3D.SrcPos(SD3D.MeasList(i,1),1) SD3D.DetPos(SD3D.MeasList(i,2),1)],[SD3D.SrcPos(SD3D.MeasList(i,1),2) SD3D.DetPos(SD3D.MeasList(i,2),2)],[SD3D.SrcPos(SD3D.MeasList(i,1),3) SD3D.DetPos(SD3D.MeasList(i,2),3)],'LineWidth',2,'Color','c');
        hold on;
    end
    plot3(SD3D.SrcPos(:,1),SD3D.SrcPos(:,2),SD3D.SrcPos(:,3),'r.','MarkerSize',20);hold on;
    plot3(SD3D.DetPos(:,1),SD3D.DetPos(:,2),SD3D.DetPos(:,3),'b.','MarkerSize',20);
    axis equal;
end

% Save #######################################
if saveFlag
    save(outname2D,'SD'); SD = SD3D;
    save(outname3D,'SD');
end

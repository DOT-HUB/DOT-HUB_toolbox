function [nirs, nirsFileName, SD3DFileName] = DOTHUB_LUMO2nirs(lumoDIR,layoutFileName,posCSVFileName,eventsCSVFileName)

%This script reads a .LUMO file (which is actually just a directory)
%and (if it is not contained within the .LUMO - which it should be)
%also loads the associated .JSON layout from a seperate file.

%If subject-specific registration data (from polhemus or other method) is available
%in .cvs format, this too can be parsd. If it is parsed, the function
%DOTHUB_LUMOPolhemus2SD is called to determine optode positions from the
%polhemus data. The resulting 3D optode positions overwrite the default 3D
%data from the .JSON file.

%This function then outputs the .nirs variables and (if saveFlag == 1 or isempty)
%produces a .nirs file with the same name as the .LUMO. Note re-running this function
%will overwrite previous .nirs derived from the parse LUMO file. Note that the .nirs
%file contains both 2D (SD) and 3D (SD_3D) info.

%######################## INPUTS ##########################################

% lumoDIR :     The path of the .LUMO data directory (if not parsed, requested)

% layoutFileName    :   (Optional) The path of the .json layout file (only required if not
%                       present inside .LUMO directory). If layoutFile variable is
%                       parsed it is assued to take precedence over any .json within
%                       the .LUMO directory. If layoutFile is not parsed (or parsed empty)
%                       and it is not present in the .LUMO direcroty, it will be requested.

% posCSVFileName    :   (Optional) The path of a .csv file containing the 3D positioning data associated
%                       with this LUMO recording. This data is a four-column CSV of
%                       position label, then x, y, z coordinate (assumed in cm).
%                       The first five rows should be Nz, Iz, Ar, Al, Cz. The
%                       following rows should be SrcA, SrcB, SrcC of each of the N
%                       tiles in turn. If polhemusCSV is not parsed, the 3D
%                       information in the .json layout file is assumed.

% eventsCSVFileName :   (Optional) The path of a .csv file containing two columns: time of event
%                       (in seconds) and event label. If the eventsCSV is parsed,
%                       it takes precedence over the contents of the .lumo file.

%######################## OUTPUTS #########################################

% nirs                      :   the contents of the resulting .nirs file as a structure.

% nirsFileName              :   Full path and name of resulting .nirs file

% nirsFileName.nirs         :   A .nirs file saved in the same location as the 
%                               .LUMO variable. Has extra variable: SD_3D, which
%                               is the same as SD, but with 3D source and detector 
%                               locations.  SD_3D also contains a 'Landmarks' field
%                               that contains the five cranial landmark locations.

%######################## Dependencies ####################################
%This script requires the TOML library for Matlab, which is available here:
%https://github.com/g-s-k/matlab-toml

% #########################################################################
% RJC, UCL, Summer 2019
%
% ############################# Updates ###################################
% RJC 20190815 - Updated to correctly re-arrange chan_list and data to match
% homer2 arrangement (sort by wavelength, source, detector)
% RJC 20190829 - Updated to work when events field is empty
% RJC 20191128 - Updated to handle multiple binary files (long datasets)
% RJC 20191202 - Updated to parse source powers into SD.SrcPowers (organized a la Homer2)
% Also note that the plan is to cap the max channel length in the .json layout files to
% 60 mm by default, which will trickle down to this conversion.
% RJC 20200217 - Tidied for GITHUB first commit.
% EGJ 20200923 - Updated to handle zero intensity values
% #########################################################################

% #########################################################################
% TO DO:
% Update event marking options to allow both .csv and automated/serial
% #########################################################################

% MANAGE VARIABLES  #######################################################
if ~exist('lumoDIR','var')
    disp('Select .LUMO directory...');
    lumoDIR = uigetdir(pwd,'Select .LUMO directory');
elseif isempty(lumoDIR)
    disp('Select .LUMO directory...');
    lumoDIR = uigetdir(pwd,'Select .LUMO directory');
elseif ~exist(lumoDIR,'dir')
    disp('Specified directory not found, please select .LUMO directory...');
    lumoDIR = uigetdir(pwd,'Select .LUMO directory');
end

%Ensure full path
[lumoPath,lumoName,ext] = fileparts(lumoDIR);
if ~strcmpi(ext,'.lumo')
    error('The lumoDIR input must point to a .LUMO directory');
end
if isempty(lumoPath)
    lumoDIR = fullfile(pwd,lumoDIR);
    [lumoPath,lumoName,~] = fileparts(lumoDIR);
end

if ~exist('layoutFileName','var') %Not parsed so check exists in .LUMO
    jsonTmp = dir([lumoDIR '/layout.json']);
    if isempty(jsonTmp) %Not contained in .LUMO, so load
        disp('Layout file not parsed or found, please select .json layout file...');
        [filename, pathname, ~] = uigetfile({'*.json';'*.JSON'},'Select .json layout file');
        layoutFileName = [pathname '/' filename];
    else
        disp('Using default layout.json within .LUMO directory');
        layoutFileName = [lumoDIR '/' jsonTmp(1).name];
    end
elseif isempty(layoutFileName)
 jsonTmp = dir([lumoDIR '/layout.json']);
    if isempty(jsonTmp) %Not contained in .LUMO, so load
        disp('Layout file not parsed or found, please select .json layout file...');
        [filename, pathname, ~] = uigetfile({'*.json';'*.JSON'},'Select .json layout file');
        layoutFileName = [pathname '/' filename];
    else
        disp('Using layout.json within .LUMO directory');
        layoutFileName = [lumoDIR '/' jsonTmp(1).name];
    end
end

polhemusFlag = 0;
if exist('posCSVFileName','var')
    if ~isempty(posCSVFileName)
        polhemusFlag = 1;
    end
end

eventsFileFlag = 0;
if exist('eventsCSVFileName','var')
    if ~isempty(eventsCSVFileName)
        eventsFileFlag = 1;
    end
end

% LOAD DATA  ###############################################################
%Load toml contents
disp('Loading .toml contents...');
metadata = toml.read([lumoDIR '/metadata.toml']);
recordingdata = toml.read([lumoDIR '/' metadata.file_names.recordingdata_file]); %To be updated to recording.toml
events = toml.read([lumoDIR '/' metadata.file_names.event_file]);
if isfield(metadata.file_names,'hardware_file') %Updated file naming structure
    hardware = toml.read([lumoDIR '/' metadata.file_names.hardware_file]);
elseif ~contains(metadata.file_names.layout_file,'json') %Old file naming structure (layout_file is hardware.toml!)
    hardware = toml.read([lumoDIR '/' metadata.file_names.layout_file]);
else
    error('No hardware.toml or layout.toml found in LUMO directory');
end

%Load layout JSON
layoutData = jsondecode(fileread(layoutFileName));

%Load log file
logFileID = fopen([lumoDIR '/' metadata.file_names.log_file]);
logtxt = textscan(logFileID,'%s');

%Load intensity binaries
disp('Loading intensity data');
nIntFiles = length(metadata.intensity_files);
for i = 1:nIntFiles
    disp(['Loading intensity data file ' num2str(i) ' of ' num2str(nIntFiles) '...']);
    fname = metadata.intensity_files{1,i}.file_name;
    %time_range = metadata.intensity_files{1,i}.time_range;
    fileID = fopen([lumoDIR '/' fname]);
    int(i).magic = fread(fileID,1);
    int(i).major_version = fread(fileID,1);
    int(i).minor_version = fread(fileID,1);
    int(i).patch = fread(fileID,1,'uint8',4);
    int(i).recordings_per_frame = fread(fileID,1,'uint64');
    int(i).nframe = fread(fileID,1,'uint64',24); %skip checksum
    int(i).data = fread(fileID,[int(i).recordings_per_frame int(i).nframe],'float32')';
    fclose(fileID);
end

% #########################################################################
% COMBINE MULTIPLE INTENSITY FILES and yield .nirs-style 'd' variable
intDataAll = [];
for i = 1:length(int)
    intDataAll = [intDataAll; int(i).data];
end
% #########################################################################

% Extract useful things for conversion ####################################
nChan = recordingdata.variables.n_chans;
chansList = reshape(recordingdata.variables.chans_list,3,nChan)';
fs = recordingdata.variables.framerate;
nFrames = recordingdata.variables.number_of_frames;
nDocks = size(layoutData.docks,1);
nodes = recordingdata.variables.nodes;
nNodes = length(nodes);

% Define t ################################################################
t = (0:1/fs:nFrames/fs - 1/fs)';

% Define SD ###############################################################
SD.nSrcs = recordingdata.variables.n_srcs;
SD.nDets = recordingdata.variables.n_dets;
SD.Lambda = recordingdata.variables.wavelength;
SD.SpatialUnit = 'mm';
MLAtmp = recordingdata.variables.chans_list_act'; %temporary saturation list

% Now determine optode positions from 2D information in layout JSON file
% Also use same loop to extract source power information
SD.SrcPos = zeros(SD.nSrcs,3);
SD.DetPos = zeros(SD.nDets,3);
SD.SrcPowers = zeros(SD.nSrcs*2,2);
SD.SrcPowers(1:SD.nSrcs,2) = 1;
SD.SrcPowers(SD.nSrcs+1:end,2) = 2;
for n = 1:nNodes
    nid = nodes(n);
    for det = 1:4
        SD.DetPos(det+(n-1)*4,1) = layoutData.docks(nid).optodes(det).coordinates_2d.x;
        SD.DetPos(det+(n-1)*4,2) = layoutData.docks(nid).optodes(det).coordinates_2d.y;
    end
    for src = 1:3
        SD.SrcPos(src+(n-1)*3,1) = layoutData.docks(nid).optodes(src+4).coordinates_2d.x;
        SD.SrcPos(src+(n-1)*3,2) = layoutData.docks(nid).optodes(src+4).coordinates_2d.y;
        
        SD.SrcPowers(src+(n-1)*3,1) = hardware.Hub.Group.Node{1,n}.Source{1,src*2-1}.Source_power;
        SD.SrcPowers(src+(n-1)*3 + SD.nSrcs,1) = hardware.Hub.Group.Node{1,n}.Source{1,src*2}.Source_power;
    end
end

%Make SD.MeasList
SD.MeasList = ones(size(chansList,1),4);
SD.MeasList(:,1:2) = chansList(:,1:2);
SD.MeasList(:,4) = chansList(:,3);
[SD.MeasList, tmpInd] = sortrows(SD.MeasList,[4,1,2]);
SD.MeasListAct = ones(size(SD.MeasList,1),1);

% ########## Use sorted ML order to correctly sort data and measlistact ##
d = intDataAll(:,tmpInd);
MLAtmp = MLAtmp(tmpInd);

% ########## Consider channels flagged as saturated ######################
% The default LUMO chan_list_act marks channels that saturate at any point
% in the recording. This can be unfairly conservative as transient motion
% can cause a transient saturation flag. For now, the saturation list is
% saved into SD.MeasListActSat, but the default is set to all active.
SD.MeasListActSat = MLAtmp;
% ########################################################################

% Include 3D information for future
if polhemusFlag %If polhemus information is parsed, calculate S-D positions from that file
    fprintf(['Using ' posCSVFileName ' to define SD3D...\n']);
    [SD_POL, SD3DFileName] = DOTHUB_LUMOpolhemus2SD3D(posCSVFileName); %This line saves the .SD3D
    SD_POL.MeasListAct = SD.MeasListAct;
    
    if nNodes<nDocks %Crop SD file accordingly if the number of docks populated in the datafile differs from the number in the SD_3D data
        SD3D = SD; %Define based on SD which contains correct measlist
        for n = 1:nNodes
            nid = nodes(n);
            for det = 1:4
                SD3D.DetPos(det+(n-1)*4,:) = SD_POL.DetPos(det+(nid-1)*4,:);
            end
            for src = 1:3
                SD3D.SrcPos(src+(n-1)*3,:) = SD_POL.SrcPos(src+(nid-1)*3,:);
            end
        end
        SD3D.Landmarks = SD_POL.Landmarks;
        
        %Plot cropped array
        f2 = figure;
        set(f2,'Name','Final (SUBSET) 3D Array Layout');
        for i = 1:size(SD3D.SrcPos,1)
            plot3(SD3D.SrcPos(i,1),SD3D.SrcPos(i,2),SD3D.SrcPos(i,3),'r.','MarkerSize',30);hold on;
            text(SD3D.SrcPos(i,1),SD3D.SrcPos(i,2)+3,SD3D.SrcPos(i,3),['S' num2str(i)],'Color','r');
        end
        for i = 1:size(SD3D.DetPos,1)
            plot3(SD3D.DetPos(i,1),SD3D.DetPos(i,2),SD3D.DetPos(i,3),'b.','MarkerSize',30);hold on;
            text(SD3D.DetPos(i,1),SD3D.DetPos(i,2)+3,SD3D.DetPos(i,3),['D' num2str(i)],'Color','b');
        end
        plotmesh(SD3D.Landmarks,'g.','MarkerSize', 30);hold on;
        landmarkLabels = {'Nz','Iz','Ar','Al','Cz'};
        for i = 1:size(SD3D.Landmarks,1)
            text(SD3D.Landmarks(i,1),SD3D.Landmarks(i,2)+3,SD3D.Landmarks(i,3)+3,landmarkLabels{i});
        end
        axis equal
        xlabel('X (mm)');ylabel('Y (mm)');zlabel('Z (mm)');
        title('Final Subset Array') 
    else
        SD3D = SD_POL;
    end
    
else %Assume 3D contents of layout file remains and save as SD3D
    
    SD3D = SD;
    for n = 1:nNodes
        nid = nodes(n);
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
    if isfield(layoutData,'Landmarks') %Nasion, Inion, Ar, Al, Cz
        for i = 1:size(layoutData.Landmarks,1)
            SD3D.Landmarks(i,1) = layoutData.Landmarks(i).x;
            SD3D.Landmarks(i,2) = layoutData.Landmarks(i).y;
            SD3D.Landmarks(i,3) = layoutData.Landmarks(i).z;
        end
    end
    SD3DFileName = fullfile(lumoPath, [lumoName '_default.SD3D']);
    fprintf(['Saving SD3D to ' SD3DFileName ' ...\n']);
    save(SD3DFileName,'SD3D');
end

% Define stim vector, s ###################################################
if eventsFileFlag %Read events from simple 2-column CSV file.
    delimiter = ',';
    startRow = 1;
    formatSpec = '%f%s%[^\n\r]';
    fileID = fopen(eventsCSVFileName,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    if isempty(dataArray)
        warning('Specified eventsCSV file contains no events!');
        s = zeros(size(t,1),1);
        CondNames = {''};
    else
        
        timeStamp = dataArray{1};
        eventStr = dataArray{2};
        [tmp,~,occuranceInd] = unique(eventStr,'stable'); %find unique events
        CondNames = arrayfun(@cellstr,tmp);
        nCond = size(CondNames,2); %number of conditions
        s = zeros(size(t,1),nCond); %preassign s matrix
        for i = 1:nCond
            timeStampTmp = timeStamp(occuranceInd==i); %Can have multiple entries
            [~,indTmp] = min(abs(repmat(t,1,length(timeStampTmp))-repmat(timeStampTmp,size(t,1),1)));
            s(indTmp,i) = 1;
        end
    end
    
else
    
    %Now check if there are any serial events
    if isempty(fieldnames(events))
        s = zeros(size(t,1),1); %return empty s vector
        CondNames = {''};
    else
        for i = 1:size(events.events,2)
            timeStamp(i) = str2num(events.events{1,i}.Timestamp)*1e-3; %to seconds
            eventStr{i} = events.events{1,i}.name;
        end
        [tmp,~,occuranceInd] = unique(eventStr,'stable'); %find unique events, maintain order in which they occured
        CondNames = arrayfun(@cellstr,tmp);
        nCond = size(CondNames,2); %number of conditions
        s = zeros(size(t,1),nCond); %preassign s matrix
        for i = 1:nCond
            timeStampTmp = timeStamp(occuranceInd==i); %Can have multiple entries
            [~,indTmp] = min(abs(repmat(t,1,length(timeStampTmp))-repmat(timeStampTmp,size(t,1),1)));
            s(indTmp,i) = 1;
        end
    end
end

% Define AUX matrix #######################################################
aux = zeros(length(t),8);

% Convert zero intensity values to noise floor estimate ###################
if size(d,1) == 1 % if only one sample
    mnD = d;
else
    mnD = mean(d);
end
dists_3D = DOTHUB_getSDdists(SD3D);
SDS_noise = 70; % SD channels greater than this (mm) are deemed to be in the the noise floor
n_zeros = sum(d==0,'all');
if n_zeros > 0
    if max(dists_3D) >= SDS_noise
        noisefloorest = mean(mnD(dists_3D>SDS_noise));
        d(d == 0) = noisefloorest;
        disp(['Warning - zero intensity values have been converted to noise floor estimate'])
    else
        d(d == 0) = 1e-6;
        disp(['Warning - zero intensity values have been converted to 1e-6'])
    end
end

% OUTPUT ##################################################################
nirs.SD = SD;
nirs.SD3D = SD3D;
nirs.aux = aux;
nirs.d = d;
nirs.s = s;
nirs.t = t;
nirs.CondNames = CondNames;

nirsFileName = fullfile(lumoPath, [lumoName '.nirs']);
fprintf(['Saving file to ' nirsFileName ' ...\n']);
save(nirsFileName,'-struct','nirs','-v7.3');



